"""Build the 3D instance set in data3d/instances.

  E3d<size>.<k>   uniform:   codes/portgen3d.c
  C3d<size>.<k>   clustered: codes/portcgen3d.c (standard deviation 1e6 / N^(1/3))
  CT3d<size>.<k>  tight clustered: as C3d with half the standard deviation.  Not part of
                  the original plan: C3d is much less clustered than the 2D C maps
                  (check_clumpiness3d.py), and factor 0.5 gives the 2D C maps' nearest-
                  neighbour ratio (0.44 against 0.42-0.44)
  star instances  hard links to the files bench/get_star.py downloaded

Sizes 1k, 10k, 100k with seeds .0 .1 .2 and 1M with .0.  Seeds follow the
DIMACS pattern of the 2D set: .0 is N itself, .k is N * 10 + k.

    py bench/gen_instances3d.py            # then bench/get_star.py's files are linked in
"""

import csv
import os
import subprocess

from common import ROOT

SIZES = [("1k", 1000, 3), ("10k", 10000, 3), ("100k", 100000, 3), ("1M", 1000000, 1)]
OUT = ROOT / "data3d" / "instances"


def generate(gen, n, seed, out, extra=()):
    if out.exists():
        return
    tmp = out.with_name(out.name + ".part")
    with open(tmp, "wb") as f:
        subprocess.run([str(ROOT / "bin" / f"{gen}.exe"), str(n), str(seed), *map(str, extra)],
                       stdout=f, check=True)
    tmp.replace(out)
    print(f"  {out.name:12} {gen} {n} {seed} {' '.join(map(str, extra))}", flush=True)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for prefix, gen, extra in (("E3d", "portgen3d", ()), ("C3d", "portcgen3d", ()),
                               ("CT3d", "portcgen3d", (0.5,))):
        for label, n, seeds in SIZES:
            for k in range(seeds):
                generate(gen, n, n if k == 0 else n * 10 + k, OUT / f"{prefix}{label}.{k}", extra)
    for r in csv.DictReader(open(ROOT / "refs" / "star_bounds.csv", newline="")):
        src, dst = ROOT / "data3d" / "star" / r["file"], OUT / f"{r['instance']}.tsp"
        if src.exists() and not dst.exists():
            os.link(src, dst)
            print(f"  {dst.name:12} linked from data3d/star/{r['file']}")
    print("3D instances are in", OUT)


if __name__ == "__main__":
    main()
