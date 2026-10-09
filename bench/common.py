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


def pin_to_core0():
    if not _kernel32.SetProcessAffinityMask(_kernel32.GetCurrentProcess(), 1):
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


def run_measured(cmd, timeout=None):
    """Run cmd; return user/kernel/wall seconds, peak MB, exit code, output."""
    env = dict(os.environ, OMP_NUM_THREADS="1")
    start = time.perf_counter()
    proc = subprocess.Popen([str(c) for c in cmd], stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, text=True, env=env)
    timed_out = False
    try:
        out, err = proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        proc.kill()
        out, err = proc.communicate()
        timed_out = True
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
