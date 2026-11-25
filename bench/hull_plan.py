"""Write the run configs for the hull-ring study (docs/hull_tsp.md).

Every config is an ordinary bench/run.py config writing to
results/2d/hull/runs.csv, so runs are resumable and every tour is checked by
codes/length.exe.  Runs that already exist in results/2d/box2/runs.csv (same
program, options and map) are copied in by bench/import_runs.py first.

  main     E and C maps, 1k and 10k (.0 .1 .2): every setting below
  big      100k (.0 .1 .2) and 1M (.0): the same settings
  r1       10k and 100k: hull, box and greedy settings with one cleanup round
           (TO R1 and T R1)
  shapes   annulus, circles, blob, thin rectangle at 10k (.0 .1 .2)

Settings (names as in results/2d/box2 where they overlap):
  hull:<merge>          rings + merge H0 or M1-M5 + TO Rn
  curve:hull            hull order alone, no merge optimisation, no cleanup
  box:<curve>:<merge>:t10   the 7 curves x M1-M5, TO Rn
  greedy+2opt+oropt     greedy TO Rn (its start_length is greedy alone)
  curve:hilbert         plain Hilbert order, no cleanup
  <name>:TO R1, <name>:T R1, greedy:TO R1, greedy:T R1   one cleanup round

    py bench/hull_plan.py
    py bench/import_runs.py --from results/2d/box2/runs.csv --to-config bench/hull_main.json --alg ...
    py bench/run.py --config bench/hull_main.json        (then big, r1, shapes)
"""

import json

from common import ROOT

CURVES = ["hilbert", "moore", "sierpinski", "knopp", "peano", "morton", "strip"]
MERGES = ["M1", "M2", "M3", "M4", "M5"]
HULL_MERGES = ["H0"] + MERGES
EXE = "bin/boxtsp2d.exe"
DIRS = {"dimacs": "C:/Users/Natha/tsp-instances-box", "shapes": "C:/Users/Natha/tsp-instances-shapes"}


def entry(name, args):
    return {"name": name, "exe": EXE, "randomized": False, "seeds": [], "geometric_only": True,
            "max_n": 10000000, "timeout_s": 3600, "args": [str(a) for a in args]}


def hull(merge, cleanup="oropt", r1=False):
    suffix = "" if not r1 else (":TO R1" if cleanup == "oropt" else ":T R1")
    return entry(f"hull:{merge}{suffix}", ["--start", "boxes", "--split", "rings", "--merge", merge,
                                           "--cleanup", cleanup] + (["--rounds", "1"] if r1 else []))


def box(curve, merge, cleanup="oropt", r1=False):
    suffix = "" if not r1 else (":TO R1" if cleanup == "oropt" else ":T R1")
    return entry(f"box:{curve}:{merge}:t10{suffix}", ["--start", "boxes", "--curve", curve, "--merge", merge,
                                                      "--t", 10, "--cleanup", cleanup]
                 + (["--rounds", "1"] if r1 else []))


def greedy(cleanup="oropt", r1=False):
    if not r1:
        return entry("greedy+2opt+oropt", ["--start", "greedy", "--cleanup", "oropt"])
    return entry("greedy:TO R1" if cleanup == "oropt" else "greedy:T R1",
                 ["--start", "greedy", "--cleanup", cleanup, "--rounds", "1"])


MAIN = ([hull(m) for m in HULL_MERGES]
        + [entry("curve:hull", ["--start", "hull", "--cleanup", "none"])]
        + [box(c, m) for c in CURVES for m in MERGES]
        + [greedy(), entry("curve:hilbert", ["--start", "curve", "--curve", "hilbert", "--cleanup", "none"])])
R1 = [f(*a, cleanup=cl, r1=True) for cl in ("oropt", "2opt")
      for f, a in [(hull, (m,)) for m in HULL_MERGES] + [(box, (c, m)) for c in CURVES for m in MERGES]]
R1 += [greedy("oropt", True), greedy("2opt", True)]


def maps(sizes, seeds=(0, 1, 2)):
    return [f"{p}{s}.{k}" for s in sizes for p in "EC" for k in (seeds if s != "1M" else (0,))]


def write(stage, folder, instances, algs):
    cfg = {"instance_dir": DIRS[folder], "bounds_file": "refs/bounds", "length_exe": "codes/length.exe",
           "results_csv": "results/2d/hull/runs.csv", "instances": instances, "algorithms": algs}
    path = ROOT / "bench" / f"hull_{stage}.json"
    path.write_text(json.dumps(cfg, indent=2))
    print(f"Wrote {path.relative_to(ROOT)}: {len(algs)} settings x {len(instances)} maps = "
          f"{len(algs) * len(instances)} runs")


def main():
    write("main", "dimacs", maps(["1k", "10k"]), MAIN)
    write("big", "dimacs", maps(["100k", "1M"]), MAIN)
    write("r1", "dimacs", maps(["10k", "100k"]), R1)
    write("shapes", "shapes", [f"{f}10k.{k}" for f in ("ann", "circ", "blob", "rect") for k in range(3)], MAIN)


if __name__ == "__main__":
    main()
