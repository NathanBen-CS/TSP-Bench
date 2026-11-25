"""Our own Held-Karp 1-tree bounds (solvers/common/hkbound.c) for the 3D maps,
and the 2D validation against the published DIMACS bounds.

    py bench/hk_bounds3d.py            # writes refs/hk_own.csv (resumable, parallel)

Every row is OUR bound: a proven lower bound on the optimal tour, at most the
true Held-Karp bound.  The 2D rows have the published Held-Karp bound (refs/bounds) next to
them: ours must never be above it, and shows how close it gets.  The star rows
have the published optimal or best tour next to them.
"""

import csv
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from common import ROOT, load_bounds, load_config

OUT = ROOT / "refs" / "hk_own.csv"
COLS = ["instance", "dims", "n", "hk_own", "rounds", "missing_last", "period", "seconds",
        "greedy_TO_R1", "published_hk_or_tour", "own_vs_published_pct"]
JOBS = [(2, m) for m in ("E1k.0", "E1k.1", "C1k.0", "C1k.1", "E10k.0", "C10k.0", "E100k.0", "C100k.0")]
JOBS += [(3, f"{p}{s}.{k}") for s in ("1k", "10k", "100k") for p in ("E3d", "C3d", "CT3d") for k in range(3)]
JOBS += [(3, m) for m in ("star1k", "star10k", "kj37859")]


def main():
    d2 = Path(load_config(ROOT / "bench" / "box_config.json")["instance_dir"])
    d3 = ROOT / "data3d" / "instances"
    bounds = load_bounds(ROOT / "refs" / "bounds")
    stars = {r["instance"]: float(r["best_tour"])
             for r in csv.DictReader(open(ROOT / "refs" / "star_bounds.csv", newline=""))}
    done = {r["instance"] for r in csv.DictReader(open(OUT, newline=""))} if OUT.exists() else set()

    def run(job):
        dims, name = job
        if name in done:
            return
        path = (d2 if dims == 2 else d3) / name
        if not path.exists():
            path = path.with_name(name + ".tsp")
        period = 500 if "100k" in name else 1000 if name == "kj37859" else 2000
        rounds = 5 if "100k" in name or name == "kj37859" else 8
        r = subprocess.run([str(ROOT / "bin" / f"hkbound{dims}d.exe"), str(path), "--period", str(period),
                            "--rounds", str(rounds)], capture_output=True, text=True)
        o = dict(re.findall(r"^(\w+)=(\S+)", r.stdout, re.M))
        pub = bounds.get(name, (None,))[0] if dims == 2 else stars.get(name)
        miss = [v for k, v in o.items() if k.startswith("missing_edges")]
        row = {"instance": name, "dims": dims, "n": o["n"], "hk_own": o["hk_bound"], "rounds": len(miss),
               "missing_last": miss[-1], "period": period, "seconds": o["time"], "greedy_TO_R1": o["tour"],
               "published_hk_or_tour": pub or "",
               "own_vs_published_pct": f"{(float(o['hk_bound']) / pub - 1) * 100:.4f}" if pub else ""}
        new = not OUT.exists()
        with open(OUT, "a", newline="") as f:
            w = csv.DictWriter(f, fieldnames=COLS)
            if new:
                w.writeheader()
            w.writerow(row)
        print(name, row["hk_own"], row["own_vs_published_pct"], flush=True)

    with ThreadPoolExecutor(8) as pool:
        list(pool.map(run, JOBS))


if __name__ == "__main__":
    main()
