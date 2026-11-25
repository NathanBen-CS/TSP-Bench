"""Write the run configs for the 3D study.  Same idea as box2_plan.py: each stage
is an ordinary bench/run.py config, so runs are resumable and every tour is
checked by bin/length3d.  All stages write to results/3d/runs.csv.

    py bench/box3_plan.py                  # writes bench/box3_*.json
    py bench/run.py --config bench/box3_baseline.json      (and so on)

  baseline   greedy, T R1, TO R1, T Rn, TO Rn, plain curve orders, nn: every map
  main       every curve x merge M1-M5, t = 10, cleanup Rn (records none, T Rn, TO Rn)
  r1         the same settings with one cleanup round (T R1 and TO R1)
  tsweep     t = 4 ... 12 at 10k and on star10k
  ksweep     K = 8, 16, 24, 32 at 10k and 100k, baseline and box settings
  big*       1M cities and gaia1: hilbert3 with every merge, the M1 and M5 settings of the
             other curves; bigrest (lengths only) goes to runs_untimed.csv
  extra, r1extra, kd, kdscale, bigF    Stage F: moore3, sierpinski3, kd-median boxes

Names: box:<curve>:<merge>:t<t>[:k<K>][:T R1|:TO R1], kd:... for kd-median boxes,
curve:<curve>, and the baselines' 2D names (greedy here = boxtsp3d's own greedy
start, because the DIMACS greedy program is 2D only).
"""

import json

from common import ROOT

CURVES = ["hilbert3", "morton3", "peano3", "snake3"]          # required curves
EXTRA = ["moore3", "sierpinski3"]                             # Stage F curves
MERGES = ["M1", "M2", "M3", "M4", "M5"]
SIZES = {"1k": 3, "10k": 3, "100k": 3, "1M": 1}
STARS = ["star100", "star1k", "star10k", "kj37859", "hyg109399", "hyg119614", "star250k", "gaia1"]
EXE = "bin/boxtsp3d.exe"


def maps(sizes, prefixes=("E3d", "C3d", "CT3d")):
    return [f"{p}{s}.{k}" for s in sizes for p in prefixes for k in range(SIZES[s])]


ALL_MAPS = maps(["1k", "10k", "100k", "1M"]) + STARS
UP_TO_250K = maps(["1k", "10k", "100k"]) + [s for s in STARS if s != "gaia1"]


def entry(name, args=None, exe=EXE, max_n=10000000):
    e = {"name": name, "exe": exe, "randomized": False, "seeds": [],
         "geometric_only": True, "max_n": max_n, "timeout_s": 3600}
    if args:
        e["args"] = [str(a) for a in args]
    return e


def box(curve, merge, t=10, k=16, r1=None, kd=False):
    name = f"{'kd' if kd else 'box'}:{curve}:{merge}:t{t}" + (f":k{k}" if k != 16 else "")
    args = ["--start", "boxes", "--curve", curve, "--merge", merge, "--t", t]
    if kd:
        args += ["--split", "kd"]
    if k != 16:
        args += ["--k", k]
    if r1:
        name += f":{r1} R1"
        args += ["--cleanup", "2opt" if r1 == "T" else "oropt", "--rounds", "1"]
    else:
        args += ["--cleanup", "oropt"]
    return entry(name, args)


def greedy_rn(k=16):
    return entry("greedy+2opt+oropt" if k == 16 else f"greedy TO Rn:k{k}",
                 ["--start", "greedy", "--cleanup", "oropt"] + (["--k", k] if k != 16 else []))


BASE = [entry("greedy", ["--start", "greedy", "--cleanup", "none"]),
        entry("twoopt", exe="bin/twoopt3d.exe"),
        entry("twoopt_or", exe="bin/twoopt_or3d.exe"),
        entry("greedy+2opt", ["--start", "greedy", "--cleanup", "2opt"]),
        greedy_rn()]


def write(stage, instances, algs, results="results/3d/runs.csv"):
    cfg = {"instance_dir": "data3d/instances", "bounds_file": "refs/bounds3d",
           "length_exe": "bin/length3d.exe", "results_csv": results,
           "instances": instances, "algorithms": algs}
    path = ROOT / "bench" / f"box3_{stage}.json"
    path.write_text(json.dumps(cfg, indent=1))
    print(f"Wrote {path.relative_to(ROOT)}: {len(algs)} settings x {len(instances)} maps = "
          f"{len(algs) * len(instances)} runs")


def main():
    plain = [entry(f"curve:{c}", ["--start", "curve", "--curve", c, "--cleanup", "none"])
             for c in CURVES + EXTRA]
    write("baseline", ALL_MAPS, BASE)
    write("orders", ALL_MAPS, plain + [entry("nn", exe="bin/nn3d.exe", max_n=40000)])
    write("main", UP_TO_250K, [box(c, m) for c in CURVES for m in MERGES])
    write("extra", UP_TO_250K, [box(c, m) for c in EXTRA for m in MERGES])
    write("r1", UP_TO_250K, [box(c, m, r1="TO") for c in CURVES for m in MERGES])
    write("r1t", UP_TO_250K, [box(c, m, r1="T") for c in CURVES for m in MERGES])
    write("r1extra", UP_TO_250K, [box(c, m, r1=r) for c in EXTRA for m in MERGES for r in ("T", "TO")])
    ts = [box(c, m, t) for c in ("hilbert3", "peano3", "morton3", "sierpinski3") for m in MERGES for t in range(4, 13)]
    write("tsweep", maps(["10k"]) + ["star10k"], ts)
    ks = [greedy_rn(k) for k in (8, 16, 24, 32)]
    ks += [box(c, m, 10, k) for c in ("hilbert3", "sierpinski3") for m in MERGES for k in (8, 16, 24, 32)]
    write("ksweep", maps(["10k", "100k"]) + ["star10k", "kj37859", "hyg109399"], ks)
    big = [f"{p}1M.0" for p in ("E3d", "C3d", "CT3d")] + ["gaia1"]
    timed = [box("hilbert3", m) for m in MERGES] + [box(c, m) for c in ("morton3", "peano3") for m in ("M1", "M5")]
    write("big", big, timed)
    write("bigr1", big, [box("hilbert3", m, r1=r) for m in MERGES for r in ("TO", "T")])
    write("bigF", big, [box("moore3", m) for m in ("M1", "M5")] + [box("sierpinski3", "M5")]
          + [box("hilbert3", "M5", 12, kd=True)])
    # the rest of the matrix at 1M, for lengths only: several runs at once on other
    # cores, so their times are not comparable and go to a separate file
    write("bigrest", big, [box(c, m) for c in ("morton3", "peano3") for m in ("M2", "M3", "M4")],
          "results/3d/runs_untimed.csv")
    kd = [box("hilbert3", m, t, kd=True) for m in MERGES for t in (4, 6, 8, 10, 12)]
    write("kd", maps(["10k"]) + ["star10k", "kj37859"], kd)
    kd = [box("hilbert3", m, t, kd=True) for m in MERGES for t in (10, 12)]
    write("kdscale", maps(["1k", "100k"]) + ["hyg109399", "star250k"], kd)


if __name__ == "__main__":
    main()
