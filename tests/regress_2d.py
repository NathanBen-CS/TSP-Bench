"""2D regression test: every 2D tour must stay exactly as it was before the
3D restructure.

The solvers are deterministic, so a tour can be compared byte for byte.  The
reference (tests/regress_2d_expected.json) was recorded from a fresh build of
the original solvers/ folder, before any file was moved.  For every run it
keeps the SHA-1 of the tour file and every number the solver prints that does
not depend on the clock (lengths, boxes, glue, moves, rounds).

Runs: on E1k.0, C1k.0, E10k.0, C10k.0
  * boxtsp, 7 curves x merges M1-M5, t = 10, full cleanup   (35 per map)
  * the reference curve with t = 4 and 12, and with K = 8
  * greedy TO Rn, greedy T Rn, greedy alone, plain curve orders
  * nn, spacefill, twoopt (T R1), twoopt_or (TO R1)
and on E100k.0, C100k.0: hilbert M1, sierpinski M5 t=4, greedy TO Rn, twoopt_or.
Then boxtest: its full output must match line for line.

    py tests/regress_2d.py                      # check bin/*2d.exe (exit 0 = identical)
    py tests/regress_2d.py --quick              # 1k maps + boxtest only
    py tests/regress_2d.py --record --exe "DIR/{name}.exe"     # write the reference

Maps come from instance_dir in bench/box_config.json.  A missing map is made with the
DIMACS generators (same seeds), so the test also runs on a fresh checkout.
"""

import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
EXPECTED = ROOT / "tests" / "regress_2d_expected.json"
CURVES = ["hilbert", "moore", "sierpinski", "knopp", "peano", "morton", "strip"]
MERGES = ["M1", "M2", "M3", "M4", "M5"]
SIZES = {"1k": 1000, "10k": 10000, "100k": 100000}


def runs():
    """(key, program, map, options)"""
    for m in ["E1k.0", "C1k.0", "E10k.0", "C10k.0"]:
        for c in CURVES:
            for g in MERGES:
                yield f"box:{c}:{g}:t10", "boxtsp", m, ["--start", "boxes", "--curve", c,
                                                       "--merge", g, "--t", "10", "--cleanup", "oropt"]
            yield f"curve:{c}", "boxtsp", m, ["--start", "curve", "--curve", c, "--cleanup", "none"]
        for t in ("4", "12"):
            for g in ("M1", "M2", "M5"):
                yield f"box:hilbert:{g}:t{t}", "boxtsp", m, ["--merge", g, "--t", t]
        yield "box:hilbert:M4:t10:k8", "boxtsp", m, ["--merge", "M4", "--k", "8"]
        yield "box:hilbert:M1:t10:2opt", "boxtsp", m, ["--cleanup", "2opt"]
        yield "hilbert-only", "boxtsp", m, ["--start", "hilbert", "--cleanup", "none"]
        yield "greedy TO Rn", "boxtsp", m, ["--start", "greedy", "--cleanup", "oropt"]
        yield "greedy T Rn", "boxtsp", m, ["--start", "greedy", "--cleanup", "2opt"]
        yield "greedy", "boxtsp", m, ["--start", "greedy", "--cleanup", "none"]
        yield "greedy TO Rn k8", "boxtsp", m, ["--start", "greedy", "--k", "8"]
        for p in ("nn", "spacefill", "twoopt", "twoopt_or"):
            yield p, p, m, []
    for m in ["E100k.0", "C100k.0"]:
        yield "box:hilbert:M1:t10", "boxtsp", m, []
        yield "box:sierpinski:M5:t4", "boxtsp", m, ["--curve", "sierpinski", "--merge", "M5", "--t", "4"]
        yield "greedy TO Rn", "boxtsp", m, ["--start", "greedy", "--cleanup", "oropt"]
        yield "twoopt_or", "twoopt_or", m, []


def find_map(name, spare):
    folder = Path(json.load(open(ROOT / "bench" / "box_config.json"))["instance_dir"])
    if (folder / name).exists():
        return folder / name
    out = spare / name
    if not out.exists():                          # DIMACS README seeds: .0 is N itself
        n = SIZES[name[1:].split(".")[0]]
        gen = ROOT / "codes" / ("portgen.exe" if name[0] == "E" else "portcgen.exe")
        with open(out, "wb") as f:
            subprocess.run([str(gen), str(n), str(n)], stdout=f, check=True)
    return out


def run_one(exe, path, opts):
    fd, tour = tempfile.mkstemp(prefix="reg_", suffix=".tour")
    os.close(fd)
    try:
        r = subprocess.run([str(exe), str(path), tour] + opts, capture_output=True, text=True)
        if r.returncode:
            return {"error": (r.stderr or r.stdout).strip()[-200:]}
        out = {k: v for k, v in re.findall(r"^(\w+)=(\S+)\s*$", r.stdout, re.M)
               if not k.endswith("_time")}
        out["sha1"] = hashlib.sha1(open(tour, "rb").read()).hexdigest()
        return out
    finally:
        os.remove(tour)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exe", default=str(ROOT / "bin" / "{name}2d.exe"))
    ap.add_argument("--record", action="store_true")
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--jobs", type=int, default=6)
    args = ap.parse_args()
    exe = lambda p: Path(args.exe.replace("{name}", p))

    spare = Path(tempfile.gettempdir()) / "regress_2d_maps"
    spare.mkdir(exist_ok=True)
    todo = [r for r in runs() if not args.quick or "1k." in r[2]]
    paths = {m: find_map(m, spare) for m in sorted({r[2] for r in todo})}
    with ThreadPoolExecutor(args.jobs) as pool:
        got = dict(zip([f"{m} | {k}" for k, _, m, _ in todo],
                       pool.map(lambda r: run_one(exe(r[1]), paths[r[2]], r[3]), todo)))
    bt = subprocess.run([str(exe("boxtest"))], capture_output=True, text=True)
    got["boxtest"] = {"exit": bt.returncode, "lines": bt.stdout.splitlines()}

    if args.record:
        bad = [k for k, v in got.items() if "error" in v] + ([] if bt.returncode == 0 else ["boxtest"])
        if bad:
            sys.exit(f"Not recorded, these failed: {bad}")
        EXPECTED.write_text(json.dumps(got, indent=1))
        print(f"Recorded {len(got) - 1} tours and boxtest "
              f"({sum(l.startswith('PASS') for l in got['boxtest']['lines'])} PASS) in {EXPECTED.name}")
        return

    want = json.load(open(EXPECTED))
    bad = []
    for k, v in got.items():
        if v != want.get(k):
            bad.append(k)
            w = want.get(k, {})
            diff = {x: (w.get(x), v.get(x)) for x in set(v) | set(w) if v.get(x) != w.get(x) and x != "lines"}
            print(f"DIFFERENT  {k}: {diff}")
    missing = [] if args.quick else [k for k in want if k not in got]
    passes = sum(l.startswith("PASS") for l in got["boxtest"]["lines"])
    print(f"{len(got) - 1} tours compared, boxtest {passes} PASS; "
          f"{len(bad)} different, {len(missing)} missing")
    print("2D REGRESSION:", "FAILED" if bad or missing else "PASSED (every tour identical)")
    sys.exit(1 if bad or missing else 0)


if __name__ == "__main__":
    main()
