"""Box diagnostics in 2D and 3D, with the same program (solvers/common/boxdiag.c).

    py bench/box3_diag.py          # writes results/3d/diagnostics.csv (and prints it)

For every map: the greedy TO Rn tour (K = 16), then for each curve and t the boxes exactly as the
solver cuts them.  Columns: boxes, average cities per box, depth-capped boxes,
share of each city's 16 nearest neighbours outside its own box (outside_k), share
of cities whose single nearest neighbour is outside (outside_nearest), and how many edges of the greedy TO Rn tour
cross box borders (share of edges, and share of the tour's length).
Not timed.  2D maps come from the 2D box experiment's folder.
"""

import csv
import re
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from common import ROOT, load_config

D2 = Path(load_config(ROOT / "bench" / "box_config.json")["instance_dir"])
D3 = ROOT / "data3d" / "instances"
MAPS = [(2, f"{p}{s}.0") for s in ("1k", "10k", "100k") for p in ("E", "C")]
MAPS += [(3, f"{p}{s}.0") for s in ("1k", "10k", "100k") for p in ("E3d", "C3d", "CT3d")]
MAPS += [(3, f"{p}1M.0") for p in ("E3d", "C3d", "CT3d")]
MAPS += [(3, s) for s in ("star100", "star1k", "star10k", "kj37859", "hyg109399", "hyg119614", "star250k", "gaia1")]
CURVES = {2: ["hilbert", "sierpinski", "peano"], 3: ["hilbert3", "peano3", "sierpinski3"]}


def run(job):
    dims, name = job
    path = (D2 if dims == 2 else D3) / name
    if not path.exists():
        path = path.with_name(name + ".tsp")
    tour = Path(tempfile.gettempdir()) / f"diag_{dims}_{name}.tour"
    subprocess.run([str(ROOT / "bin" / f"boxtsp{dims}d.exe"), str(path), str(tour), "--start", "greedy"],
                   check=True, capture_output=True)
    rows = []
    for c in CURVES[dims]:
        for t in (4, 10, 12):
            o = subprocess.run([str(ROOT / "bin" / f"boxdiag{dims}d.exe"), str(path), str(tour),
                                "--curve", c, "--t", str(t)], capture_output=True, text=True, check=True).stdout
            d = dict(re.findall(r"^(\w+)=(\S+)", o, re.M))
            rows.append({"dims": dims, "map": name, "curve": c, "t": t,
                         **{k[5:]: v for k, v in d.items() if k.startswith("diag_")}})
    tour.unlink()
    print(name, flush=True)
    return rows


def main():
    with ThreadPoolExecutor(6) as pool:
        rows = [r for rs in pool.map(run, MAPS) for r in rs]
    out = ROOT / "results" / "3d" / "diagnostics.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("Wrote", out)


if __name__ == "__main__":
    main()
