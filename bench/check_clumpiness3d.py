"""How clustered are the 3D instances?  The 3D version of check_clumpiness.py.

For each map: the mean distance to the nearest other city (exact, over a sample
of up to 2,000 cities) divided by what a uniform cloud with the same number of
cities in the same volume would have, 0.554 * (V/N)^(1/3) (2D: 0.5 * (A/N)^(1/2)).
For E and C maps the volume is the 10^6 cube (square), for the stars it is
their bounding box.  1 = as even as uniform, smaller = more clustered.  The C/E ratio
is also printed; in 2D it is 0.40-0.44.  Also prints the share of C cities outside the cube.

    py bench/check_clumpiness3d.py
"""

import sys
from pathlib import Path

import numpy as np

from common import ROOT, load_config

SAMPLE = 2000


def read_coords(path, dims):
    with open(path) as f:
        for line in f:
            if line.startswith("NODE_COORD_SECTION"):
                break
        return np.loadtxt(f, usecols=range(1, dims + 1), comments="EOF")


def mean_nn(x, rng):
    idx = rng.choice(len(x), size=min(SAMPLE, len(x)), replace=False)
    total = 0.0
    for i in idx:
        d2 = ((x - x[i]) ** 2).sum(axis=1)
        d2[i] = np.inf
        total += np.sqrt(d2.min())
    return total / len(idx)


def main():
    rng = np.random.default_rng(12345)
    d3 = ROOT / "data3d" / "instances"
    d2 = Path(load_config(ROOT / "bench" / "box_config.json")["instance_dir"])
    expect = {2: lambda v, n: 0.5 * (v / n) ** 0.5, 3: lambda v, n: 0.55396 * (v / n) ** (1 / 3)}
    nn = {}
    print(f"{'map':10} {'mean NN dist':>12} {'vs even spread':>15} {'outside cube':>13}")
    names = [(d, p + s + ".0") for s in ("1k", "10k", "100k", "1M") for d, p in ((2, "E"), (2, "C"), (3, "E3d"), (3, "C3d"))]
    names += [(3, p.name) for p in sorted(d3.glob("*.tsp"))]
    for dims, name in names:
        path = (d3 if dims == 3 else d2) / name
        if not path.exists():
            continue
        x = read_coords(path, dims)
        star = name.endswith(".tsp")
        vol = np.prod(x.max(axis=0) - x.min(axis=0)) if star else 1e6 ** dims
        nn[name] = mean_nn(x, rng)
        out = "" if star else f"{((x < 0) | (x >= 1e6)).any(axis=1).mean():13.3f}"
        print(f"{name:14} {nn[name]:10.1f} {nn[name] / expect[dims](vol, len(x)):13.2f} {out}", flush=True)
    for s in ("1k", "10k", "100k", "1M"):
        print(f"{s:5} C/E nearest-neighbour ratio:  2D {nn['C'+s+'.0']/nn['E'+s+'.0']:.2f}   3D {nn['C3d'+s+'.0']/nn['E3d'+s+'.0']:.2f}")


if __name__ == "__main__":
    main()
