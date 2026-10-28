"""Check that the uniform (E*) and clustered (C*) instances really differ.

For each size, compares E<size>.0 with C<size>.0 on two measures:

  empty cells   share of the cells of a 100 x 100 grid over the instance's
                bounding box that contain no city (clumpy data leaves more
                of the map empty)
  mean NN       mean distance from a city to its nearest other city, over a
                sample of up to 2,000 cities (clumpy data packs cities
                closer together)

The test FAILS (exit code 1) unless, at every size, the clustered instance
has more empty cells AND a mean nearest-neighbour distance below 0.7 x the
uniform one.  Both numbers are printed either way.

    py bench/check_clumpiness.py --dir C:/Users/Natha/tsp-instances-box
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np

GRID = 100
SAMPLE = 2000
NN_RATIO_MAX = 0.7


def read_coords(path):
    """Coordinates of a TSPLIB NODE_COORD_SECTION file as an (n, 2) array."""
    with open(path) as f:
        for line in f:
            if line.startswith("NODE_COORD_SECTION"):
                break
        data = np.loadtxt(f, usecols=(1, 2), comments="EOF")
    return data


def empty_share(xy):
    lo, hi = xy.min(axis=0), xy.max(axis=0)
    span = np.where(hi > lo, hi - lo, 1.0)
    cell = np.minimum(((xy - lo) / span * GRID).astype(int), GRID - 1)
    occupied = np.unique(cell[:, 0] * GRID + cell[:, 1]).size
    return 1 - occupied / (GRID * GRID)


def mean_nn(xy, rng):
    idx = rng.choice(len(xy), size=min(SAMPLE, len(xy)), replace=False)
    total = 0.0
    for i in idx:                                    # exact NN of each sample city
        d2 = ((xy - xy[i]) ** 2).sum(axis=1)
        d2[i] = np.inf
        total += np.sqrt(d2.min())
    return total / len(idx)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", required=True, help="folder with E*.0 / C*.0 instances")
    args = ap.parse_args()

    folder = Path(args.dir)
    sizes = sorted({m.group(1) for p in folder.iterdir()
                    if (m := re.fullmatch(r"[EC](\d+[kM])\.0", p.name))},
                   key=lambda s: int(s[:-1]) * (1000 if s[-1] == "k" else 10**6))
    rng = np.random.default_rng(12345)
    ok = True
    print(f"{'size':6} {'type':10} {'empty cells':>12} {'mean NN dist':>13}")
    for size in sizes:
        e_path, c_path = folder / f"E{size}.0", folder / f"C{size}.0"
        if not (e_path.exists() and c_path.exists()):
            continue
        res = {}
        for label, path in (("uniform", e_path), ("clustered", c_path)):
            xy = read_coords(path)
            res[label] = (empty_share(xy), mean_nn(xy, rng))
            print(f"{size:6} {label:10} {res[label][0]:12.3f} {res[label][1]:13.0f}")
        ratio = res["clustered"][1] / res["uniform"][1]
        size_ok = res["clustered"][0] > res["uniform"][0] and ratio < NN_RATIO_MAX
        ok &= size_ok
        print(f"{'':6} {'':10} {'':12} {'':13}  NN ratio {ratio:.2f} -> "
              f"{'PASS' if size_ok else 'FAIL'}")
    if not sizes:
        print("No E*.0 / C*.0 pairs found.")
        ok = False
    print("\nClumpiness check:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
