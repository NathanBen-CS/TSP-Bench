"""Shared helpers for the TSP benchmark scripts (Windows only, stdlib only).

Measurement works like Linux `/usr/bin/time -v` + `taskset -c 0`:

  * pin_to_core0()  restricts this Python process to CPU 0.  Windows child
                    processes inherit the affinity mask, so every solver we
                    launch also runs on CPU 0 only.
  * run_measured()  launches a program and, after it exits but before its
                    process handle is closed, asks Windows for
                      - user CPU time    (GetProcessTimes)
                      - peak memory      (GetProcessMemoryInfo -> PeakWorkingSetSize,
                                          the Windows equivalent of max RSS)
                    for the whole program, file reading included.
"""

import ctypes
import json
import os
import re
import subprocess
import sys
import time
from ctypes import wintypes
from pathlib import Path

if sys.platform != "win32":
    sys.exit("These scripts measure via the Windows API; run them on Windows.")

ROOT = Path(__file__).resolve().parent.parent

# ---------------------------------------------------------------- Windows API

_kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)


class _PROCESS_MEMORY_COUNTERS(ctypes.Structure):
    _fields_ = [
        ("cb", wintypes.DWORD),
        ("PageFaultCount", wintypes.DWORD),
        ("PeakWorkingSetSize", ctypes.c_size_t),
        ("WorkingSetSize", ctypes.c_size_t),
        ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
        ("QuotaPagedPoolUsage", ctypes.c_size_t),
        ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
        ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
        ("PagefileUsage", ctypes.c_size_t),
        ("PeakPagefileUsage", ctypes.c_size_t),
    ]


_kernel32.GetCurrentProcess.restype = wintypes.HANDLE
_kernel32.SetProcessAffinityMask.argtypes = [wintypes.HANDLE, ctypes.c_size_t]
_kernel32.SetProcessAffinityMask.restype = wintypes.BOOL
_kernel32.GetProcessTimes.argtypes = [wintypes.HANDLE] + [ctypes.POINTER(wintypes.FILETIME)] * 4
_kernel32.GetProcessTimes.restype = wintypes.BOOL
_kernel32.K32GetProcessMemoryInfo.argtypes = [
    wintypes.HANDLE, ctypes.POINTER(_PROCESS_MEMORY_COUNTERS), wintypes.DWORD]
_kernel32.K32GetProcessMemoryInfo.restype = wintypes.BOOL


def _filetime_seconds(ft):
    return ((ft.dwHighDateTime << 32) | ft.dwLowDateTime) / 1e7   # 100 ns units


def pin_to_core0(core=0):
    if not _kernel32.SetProcessAffinityMask(_kernel32.GetCurrentProcess(), 1 << core):
        raise ctypes.WinError(ctypes.get_last_error())


def _process_stats(handle):
    creation, exit_, kernel, user = (wintypes.FILETIME() for _ in range(4))
    if not _kernel32.GetProcessTimes(handle, ctypes.byref(creation), ctypes.byref(exit_),
                                     ctypes.byref(kernel), ctypes.byref(user)):
        raise ctypes.WinError(ctypes.get_last_error())
    mem = _PROCESS_MEMORY_COUNTERS()
    mem.cb = ctypes.sizeof(mem)
    if not _kernel32.K32GetProcessMemoryInfo(handle, ctypes.byref(mem), mem.cb):
        raise ctypes.WinError(ctypes.get_last_error())
    return _filetime_seconds(user), _filetime_seconds(kernel), mem.PeakWorkingSetSize / 2**20


def run_measured(cmd, timeout=None, heartbeat=None, every=1.0):
    """Run cmd; return user/kernel/wall seconds, peak MB, exit code, output.

    heartbeat(elapsed_seconds) is called about every `every` seconds while the
    program runs, so callers can show progress.  Ctrl+C kills the program
    before re-raising, so no solver is left running in the background."""
    env = dict(os.environ, OMP_NUM_THREADS="1")
    start = time.perf_counter()
    proc = subprocess.Popen([str(c) for c in cmd], stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, text=True, env=env)
    timed_out = False
    try:
        while True:
            elapsed = time.perf_counter() - start
            wait = every if heartbeat else None              # how long to block
            if timeout is not None:
                left = timeout - elapsed
                if left <= 0:                                # out of time: stop it
                    proc.kill()
                    out, err = proc.communicate()
                    timed_out = True
                    break
                wait = left if wait is None else min(wait, left)
            try:
                out, err = proc.communicate(timeout=wait)   # finished: done
                break
            except subprocess.TimeoutExpired:               # still running:
                if heartbeat:                               #   report, keep waiting
                    heartbeat(time.perf_counter() - start)
    except KeyboardInterrupt:
        proc.kill()
        proc.communicate()
        raise
    wall = time.perf_counter() - start
    # proc._handle stays open until the Popen object is destroyed
    user, kernel, peak_mb = _process_stats(int(proc._handle))
    return {"user_s": user, "kernel_s": kernel, "wall_s": wall, "peak_mb": peak_mb,
            "exit_code": proc.returncode, "timed_out": timed_out,
            "stdout": out, "stderr": err}


# ---------------------------------------------------------------- config / files

def load_config(path=None):
    path = Path(path) if path else ROOT / "bench" / "config.json"
    with open(path) as f:
        return json.load(f)


def resolve(p):
    """Config paths are relative to the project root unless absolute."""
    p = Path(p)
    return p if p.is_absolute() else ROOT / p


def instance_name(path):
    """E1k.0 -> E1k.0, pr1002.tsp -> pr1002 (matches names in the bounds file)."""
    name = Path(path).name
    return name[:-4] if name.endswith(".tsp") else name


def instance_is_geometric(path):
    """True if the instance has coordinates (EUC_2D, CEIL_2D, EUC_3D), not a matrix."""
    with open(path, "rb") as f:
        for _ in range(50):
            line = f.readline().decode("ascii", "replace")
            m = re.match(r"\s*EDGE_WEIGHT_TYPE\s*:?\s*(\S+)", line)
            if m:
                return m.group(1) in ("EUC_2D", "CEIL_2D", "EUC_3D")
    return False


def instance_n(path):
    """Read DIMENSION from a TSPLIB header without loading the whole file."""
    with open(path, "rb") as f:
        for _ in range(50):
            line = f.readline().decode("ascii", "replace")
            m = re.match(r"\s*DIMENSION\s*:?\s*(\d+)", line)
            if m:
                return int(m.group(1))
    raise ValueError(f"no DIMENSION line in {path}")


def list_instances(instance_dir):
    """All instance files in the directory as (name, path, n), smallest first."""
    found = []
    for p in Path(instance_dir).iterdir():
        if p.is_file() and not p.name.startswith("."):
            try:
                found.append((instance_name(p), p, instance_n(p)))
            except ValueError:
                pass
    return sorted(found, key=lambda t: (t[2], t[0]))


def load_bounds(path):
    """bounds file lines: name HK_bound [optimal_or_best_known]."""
    bounds = {}
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 2:
                bounds[parts[0]] = (float(parts[1]),
                                    float(parts[2]) if len(parts) >= 3 else None)
    return bounds


def parse_tour_length(stdout):
    m = re.search(r"Tour length is\s+(\d+)", stdout)
    return int(m.group(1)) if m else None


def parse_extras(stdout):
    """Solver lines of the form key=value become extra result columns."""
    return dict(re.findall(r"^\s*([A-Za-z_]\w*)=(\S+)\s*$", stdout, re.MULTILINE))


# ---------------------------------------------------------------- chart names

CODE_KEY = ("Codes: bare name = that order only; B = boxes (exact loops + merge); "
            "T = 2-opt; O = Or-opt; R1 = one cleanup round; Rn = rounds until nothing improves")


def config_label(algorithm, level="full"):
    """Code name of a run for charts: <start> [B][T][O] [R1|Rn] [merge t=..].

    B = boxes (exact Held-Karp loops + merge), T = 2-opt, O = Or-opt; a bare
    name means just following that order.  Every name with a cleanup (T or O)
    says how long it ran: R1 = one round of the work queue (the bench's
    twoopt and twoopt_or), Rn = rounds repeated until a whole round makes no
    move (every boxtsp cleanup).  level picks the cleanup stage of a box run
    whose tour was recorded at several stages: "none", "2opt", "full".
    """
    stage = {"none": "", "2opt": "T Rn", "full": "TO Rn"}[level]
    fixed = {
        "greedy": "greedy",                       # DIMACS greedy, no cleanup
        "nn": "nn",
        "spacefill": "hilbert (spacefill)",       # Hilbert order, coarser grid
        "twoopt": "greedy T R1",                  # bench twoopt: greedy start, one round
        "twoopt_or": "greedy TO R1",              # bench twoopt_or
        "greedy+2opt+oropt": "greedy TO Rn",      # cleanup repeated until clean
        "greedy+2opt": "greedy T Rn",
        "hilbert-only": "hilbert",
        "boxes-only": "hilbert B M1 t=10",
        "boxes+2opt": "hilbert BT Rn M1 t=10",
        "boxes+2opt+oropt": "hilbert BTO Rn M1 t=10",
    }
    if algorithm in fixed:
        return fixed[algorithm]
    m = re.fullmatch(r"boxes\+2opt\+oropt-t(\d+)", algorithm)       # box experiment
    if m:
        return f"hilbert BTO Rn M1 t={m.group(1)}"
    m = re.fullmatch(r"box:(\w+):(M\d):t(\d+)", algorithm)          # 2D extension, 3D
    if m:
        return f"{m.group(1)} B{stage} {m.group(2)} t={m.group(3)}"
    m = re.fullmatch(r"(box|kd):(\w+):(M\d):t(\d+)(:k\d+)?(:T R1|:TO R1)?", algorithm)   # 3D
    if m:
        box = "kd-B" if m.group(1) == "kd" else "B"
        st = m.group(6)[1:] if m.group(6) else stage
        return (f"{m.group(2)} {box}{st} {m.group(3)} t={m.group(4)}"
                + (f" K={m.group(5)[2:]}" if m.group(5) else ""))
    m = re.fullmatch(r"greedy TO Rn:k(\d+)", algorithm)
    if m:
        return f"greedy TO Rn K={m.group(1)}"
    if algorithm.startswith("curve:"):
        return algorithm[6:]
    return algorithm
