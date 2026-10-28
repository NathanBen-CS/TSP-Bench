"""Build the instance set for the box-TSP experiment (Box-TSP plan, Task 1).

Uniform and clustered instances at 1k, 10k, 100k and 1M cities, at least 3
seeds of each, made with the DIMACS generators in codes/ (portgen = uniform
in the 1,000,000 x 1,000,000 square; portcgen = DIMACS clustered: N/100
centres, normal spread).  Instances that already exist in the main DIMACS
testbed are hard-linked (no extra disk space); missing seeds are generated
with the README's seed pattern (N * 10 + k), e.g.

    portgen 100000 1000002 > E100k.2        portgen 1000000 10000001 > E1M.1

Everything goes into the folder named by instance_dir in bench/box_config.json,
so the DIMACS benchmark itself is untouched.  Then the clumpiness check runs
(bench/check_clumpiness.py).

    py bench/box_instances.py
    py bench/box_instances.py --max-n 100000      # skip the 1M instances
"""

import argparse
import os
import shutil
import subprocess
import sys

from common import ROOT, load_config, resolve

SIZES = [("1k", 1000, 10), ("10k", 10000, 3), ("100k", 100000, 3), ("1M", 1000000, 3)]


def dimacs_seed(n, k):
    """README pattern: seed .0 is N itself, seed .k is N*10 + k."""
    return n if k == 0 else n * 10 + k


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=str(ROOT / "bench" / "box_config.json"))
    ap.add_argument("--main-config", default=str(ROOT / "bench" / "config.json"))
    ap.add_argument("--max-n", type=int, default=None)
    args = ap.parse_args()

    out_dir = resolve(load_config(args.config)["instance_dir"])
    main_dir = resolve(load_config(args.main_config)["instance_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"Box experiment instances: {out_dir}")

    for prefix, gen in (("E", "portgen"), ("C", "portcgen")):
        for label, n, seeds in SIZES:
            if args.max_n and n > args.max_n:
                continue
            for k in range(seeds):
                name = f"{prefix}{label}.{k}"
                dst = out_dir / name
                if dst.exists():
                    continue
                src = main_dir / name
                if src.exists():                       # DIMACS testbed file
                    try:
                        os.link(src, dst)              # hard link: no copy
                        how = "linked"
                    except OSError:
                        shutil.copyfile(src, dst)
                        how = "copied"
                else:                                  # missing seed: generate
                    seed = dimacs_seed(n, k)
                    tmp = dst.with_name(name + ".part")
                    with open(tmp, "wb") as f:
                        subprocess.run([str(ROOT / "codes" / f"{gen}.exe"), str(n), str(seed)],
                                       stdout=f, check=True)
                    tmp.replace(dst)
                    how = f"generated: {gen} {n} {seed}"
                print(f"  {name:10} {how}", flush=True)

    print()
    check = subprocess.run([sys.executable, str(ROOT / "bench" / "check_clumpiness.py"),
                            "--dir", str(out_dir)])
    sys.exit(check.returncode)


if __name__ == "__main__":
    main()
