"""DIMACS machine benchmark: run the "time greedy INSTANCE ITERATIONS" list
from codes/README (pinned to CPU 0, twice each by default) and record the user
times, so readers can rescale your solver times to their machine.

Also writes results/machine.json describing this machine and the compiler
settings used for codes/ and solvers/.

    python bench/machine_bench.py
    python bench/machine_bench.py --max-n 100000 --runs 3
"""

import argparse
import csv
import ctypes
import datetime
import json
import os
import platform
import re
import shutil
import subprocess
import winreg

from common import ROOT, instance_n, load_config, pin_to_core0, resolve, run_measured

README_BENCH = re.compile(r"^time greedy (\S+) (\d+)\s*$")


def benchmark_list():
    with open(ROOT / "codes" / "README") as f:
        return [(m.group(1), int(m.group(2))) for m in map(README_BENCH.match, f) if m]


def cpu_name():
    try:
        key = winreg.OpenKey(winreg.HKEY_LOCAL_MACHINE,
                             r"HARDWARE\DESCRIPTION\System\CentralProcessor\0")
        return winreg.QueryValueEx(key, "ProcessorNameString")[0].strip()
    except OSError:
        return platform.processor()


def ram_gb():
    class MEMORYSTATUSEX(ctypes.Structure):
        _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong),
                    ("ullTotalPhys", ctypes.c_ulonglong), ("ullAvailPhys", ctypes.c_ulonglong),
                    ("ullTotalPageFile", ctypes.c_ulonglong),
                    ("ullAvailPageFile", ctypes.c_ulonglong),
                    ("ullTotalVirtual", ctypes.c_ulonglong),
                    ("ullAvailVirtual", ctypes.c_ulonglong),
                    ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]
    m = MEMORYSTATUSEX()
    m.dwLength = ctypes.sizeof(m)
    ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m))
    return round(m.ullTotalPhys / 2**30, 1)


def gcc_version():
    gcc = shutil.which("gcc") or r"C:\msys64\ucrt64\bin\gcc.exe"
    try:
        return subprocess.run([gcc, "--version"], capture_output=True,
                              text=True).stdout.splitlines()[0]
    except (OSError, IndexError):
        return "unknown"


def makefile_cflags(path):
    m = re.search(r"^CFLAGS\s*=\s*(.*)$", path.read_text(), re.MULTILINE)
    return m.group(1).strip() if m else "unknown"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    ap.add_argument("--runs", type=int, default=2)
    ap.add_argument("--max-n", type=int, default=None,
                    help="skip benchmark instances with more cities than this")
    args = ap.parse_args()

    cfg = load_config(args.config)
    inst_dir = resolve(cfg["instance_dir"])
    greedy = ROOT / "codes" / "greedy.exe"
    pin_to_core0()

    out_dir = ROOT / "results"
    out_dir.mkdir(exist_ok=True)
    rows = []
    for name, iters in benchmark_list():
        path = inst_dir / name
        if not path.exists():
            print(f"{name}: missing (run gen_instances.py), skipped")
            continue
        n = instance_n(path)
        if args.max_n and n > args.max_n:
            continue
        for run in range(1, args.runs + 1):
            r = run_measured([greedy, path, iters])
            if r["exit_code"] != 0:
                print(f"{name}: greedy failed: {r['stderr'].strip()}")
                break
            rows.append({"instance": name, "n": n, "iterations": iters, "run": run,
                         "user_s": f"{r['user_s']:.2f}",
                         "user_s_per_iteration": f"{r['user_s'] / iters:.5f}",
                         "peak_mb": f"{r['peak_mb']:.1f}"})
            print(f"greedy {name} x{iters}  run {run}: user {r['user_s']:.2f}s")

    with open(out_dir / "machine_bench.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["instance", "n", "iterations", "run", "user_s",
                                          "user_s_per_iteration", "peak_mb"])
        w.writeheader()
        w.writerows(rows)

    info = {
        "date": datetime.date.today().isoformat(),
        "cpu": cpu_name(),
        "logical_cpus": os.cpu_count(),
        "ram_gb": ram_gb(),
        "os": platform.platform(),
        "compiler": gcc_version(),
        "cflags_codes": makefile_cflags(ROOT / "codes" / "Makefile"),
        "cflags_solvers": makefile_cflags(ROOT / "solvers" / "Makefile"),
        "pinning": "runner sets affinity mask to CPU 0; children inherit it",
        "time_measure": "GetProcessTimes user time of the whole process",
        "memory_measure": "GetProcessMemoryInfo PeakWorkingSetSize",
    }
    with open(out_dir / "machine.json", "w") as f:
        json.dump(info, f, indent=2)
    print(f"Wrote results/machine_bench.csv and results/machine.json ({info['cpu']})")


if __name__ == "__main__":
    main()
