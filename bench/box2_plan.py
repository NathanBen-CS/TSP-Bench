"""Write the run configs for the 2D box extension (curves x merges x t).

Each stage is an ordinary bench/run.py config, so runs are resumable and
every tour is checked by codes/length.exe.  All stages write to the same
results file, results/box2/runs.csv; a setting that appears in two stages
has the same name and is only run once.

Every box run uses the full cleanup (2-opt, then 2-opt + Or-opt).  The run
also records the tour before cleanup (start_length) and after 2-opt
(after_2opt_length), so one run gives all three cleanup levels.

    py bench/box2_plan.py A                       # screen at N = 10,000
    py bench/box2_plan.py B --curves hilbert,moore,peano --merges M3,M4
    py bench/box2_plan.py C --settings hilbert:M4:10,moore:M4:8
then
    py bench/run.py --config bench/box2_stageA.json      (B, C likewise)

Names: box:<curve>:<merge>:t<t>  (e.g. box:hilbert:M1:t10, the reference),
curve:<curve> (the curve's own order, no boxes, no cleanup), and the
baselines greedy+2opt+oropt plus every existing bench solver.
"""

import argparse
import json

from common import ROOT

CURVES = ["hilbert", "moore", "sierpinski", "knopp", "peano", "morton", "strip"]
MERGES = ["M1", "M2", "M3", "M4", "M5"]
TS = [4, 6, 8, 10, 12]
REFERENCE = ("hilbert", "M1", 10)
SEEDS = {"10k": [0, 1, 2], "100k": [0, 1, 2]}

EXISTING = [
    {"name": "greedy",    "exe": "codes/greedy.exe"},
    {"name": "nn",        "exe": "solvers/nn.exe", "max_n": 320000},
    {"name": "spacefill", "exe": "solvers/spacefill.exe"},
    {"name": "twoopt",    "exe": "solvers/twoopt.exe"},
    {"name": "twoopt_or", "exe": "solvers/twoopt_or.exe"},
]


def entry(name, exe, args=None, max_n=10000000):
    e = {"name": name, "exe": exe, "randomized": False, "seeds": [],
         "geometric_only": True, "max_n": max_n, "timeout_s": 3600}
    if args:
        e["args"] = args
    return e


def box(curve, merge, t):
    return entry(f"box:{curve}:{merge}:t{t}", "solvers/boxtsp.exe",
                 ["--start", "boxes", "--curve", curve, "--merge", merge,
                  "--t", str(t), "--cleanup", "oropt"])


def baselines(with_curves):
    algs = [entry(e["name"], e["exe"], max_n=e.get("max_n", 10000000)) for e in EXISTING]
    algs.append(entry("greedy+2opt+oropt", "solvers/boxtsp.exe",
                      ["--start", "greedy", "--cleanup", "oropt"]))
    if with_curves:
        algs += [entry(f"curve:{c}", "solvers/boxtsp.exe",
                       ["--start", "curve", "--curve", c, "--cleanup", "none"]) for c in CURVES]
    return algs


def write(stage, size, algs):
    cfg = {
        "instance_dir": "C:/Users/Natha/tsp-instances-box",
        "bounds_file": "refs/bounds",
        "length_exe": "codes/length.exe",
        "results_csv": "results/box2/runs.csv",
        "instances": [f"{p}{size}.{s}" for p in "EC" for s in SEEDS[size]],
        "algorithms": algs,
    }
    path = ROOT / "bench" / f"box2_stage{stage}.json"
    path.write_text(json.dumps(cfg, indent=2))
    print(f"Wrote {path.relative_to(ROOT)}: {len(algs)} settings x "
          f"{len(cfg['instances'])} instances = {len(algs) * len(cfg['instances'])} runs")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stage", choices=["A", "B", "C"])
    ap.add_argument("--curves", help="stage B: the 3 best curves, comma-separated")
    ap.add_argument("--merges", help="stage B: the 2 best merges, comma-separated")
    ap.add_argument("--settings", help="stage C: curve:merge:t finalists, comma-separated")
    args = ap.parse_args()

    if args.stage == "A":            # screen: every curve x merge, t = 10
        algs = [box(c, m, 10) for c in CURVES for m in MERGES] + baselines(True)
        write("A", "10k", algs)
    elif args.stage == "B":          # t sweep for the best curves x merges (+ reference)
        curves = args.curves.split(",")
        merges = args.merges.split(",")
        combos = [(c, m) for c in curves for m in merges]
        if (REFERENCE[0], REFERENCE[1]) not in combos:
            combos.append((REFERENCE[0], REFERENCE[1]))
        algs = [box(c, m, t) for c, m in combos for t in TS]
        write("B", "10k", algs)
    else:                            # scale: finalists + reference at N = 100,000
        settings = [s.split(":") for s in args.settings.split(",")]
        if [REFERENCE[0], REFERENCE[1], str(REFERENCE[2])] not in settings:
            settings.append([REFERENCE[0], REFERENCE[1], str(REFERENCE[2])])
        algs = [box(c, m, int(t)) for c, m, t in settings] + baselines(False)
        write("C", "100k", algs)


if __name__ == "__main__":
    main()
