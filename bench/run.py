"""Run every algorithm on every instance and append one row per run to the
results CSV (one row per seed for randomized algorithms, never averaged).

For each run:
  1. launch the solver pinned to CPU 0:  solver <instance> <tour_file> [seed]
  2. record user CPU time and peak memory of the whole program
  3. check the tour with codes/length.exe and record its length
  4. append the row, delete the tour file

Runs already in the CSV are skipped, so an interrupted run can be resumed.

    python bench/run.py
    python bench/run.py --alg nn --inst E1k.0,pr1002
"""

import argparse
import csv
import datetime
import os
import tempfile

from common import (instance_is_geometric, list_instances, load_config, parse_extras,
                    parse_tour_length, pin_to_core0, resolve, run_measured)

COLUMNS = ["algorithm", "instance", "n", "seed", "tour_length", "user_s", "kernel_s",
           "wall_s", "peak_mb", "status", "extras", "timestamp"]


def load_done(csv_path):
    if not csv_path.exists():
        return set()
    with open(csv_path, newline="") as f:
        return {(r["algorithm"], r["instance"], r["seed"]) for r in csv.DictReader(f)}


def append_row(csv_path, row):
    new = not csv_path.exists()
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with open(csv_path, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS)
        if new:
            w.writeheader()
        w.writerow(row)


def run_one(alg, inst_path, seed, length_exe):
    fd, tour = tempfile.mkstemp(prefix="tour_", suffix=".txt")
    os.close(fd)
    try:
        cmd = [resolve(alg["exe"]), inst_path, tour]
        if seed != "":
            cmd.append(seed)
        r = run_measured(cmd, timeout=alg.get("timeout_s"))
        row = {"user_s": f"{r['user_s']:.3f}", "kernel_s": f"{r['kernel_s']:.3f}",
               "wall_s": f"{r['wall_s']:.3f}", "peak_mb": f"{r['peak_mb']:.1f}",
               "tour_length": "",
               "extras": ";".join(f"{k}={v}" for k, v in parse_extras(r["stdout"]).items())}
        if r["timed_out"]:
            row["status"] = "DNF"
        elif r["exit_code"] != 0:
            row["status"] = "FAILED"
            print("    solver error:", (r["stderr"] or r["stdout"]).strip()[-300:])
        else:
            check = run_measured([length_exe, inst_path, tour])
            length = parse_tour_length(check["stdout"])
            if check["exit_code"] == 0 and length is not None:
                row["status"], row["tour_length"] = "OK", length
            else:
                row["status"] = "INVALID"
                print("    length.exe:", (check["stderr"] or check["stdout"]).strip()[-300:])
        return row
    finally:
        os.remove(tour)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    ap.add_argument("--alg", help="comma-separated algorithm names (default: all)")
    ap.add_argument("--inst", help="comma-separated instance names (default: config)")
    args = ap.parse_args()

    cfg = load_config(args.config)
    pin_to_core0()

    instances = list_instances(resolve(cfg["instance_dir"]))
    wanted = args.inst.split(",") if args.inst else cfg.get("instances", "all")
    if wanted != "all":
        instances = [t for t in instances if t[0] in wanted]
    algs = cfg["algorithms"]
    if args.alg:
        algs = [a for a in algs if a["name"] in args.alg.split(",")]

    csv_path = resolve(cfg["results_csv"])
    length_exe = resolve(cfg["length_exe"])
    done = load_done(csv_path)
    print(f"{len(algs)} algorithm(s) x {len(instances)} instance(s), pinned to CPU 0")

    for alg in algs:
        seeds = [str(s) for s in alg.get("seeds", [])] if alg.get("randomized") else [""]
        for name, path, n in instances:
            for seed in seeds:
                if (alg["name"], name, seed) in done:
                    continue
                label = f"{alg['name']:10} {name:12} n={n:<9}" + (f" seed={seed}" if seed else "")
                if n > alg.get("max_n", float("inf")) or (
                        alg.get("geometric_only") and not instance_is_geometric(path)):
                    row = {"status": "SKIPPED"}
                else:
                    print(label, "...", end=" ", flush=True)
                    row = run_one(alg, path, seed, length_exe)
                    print(row["status"], row["tour_length"],
                          f"user={row['user_s']}s peak={row['peak_mb']}MB")
                row.update(algorithm=alg["name"], instance=name, n=n, seed=seed,
                           timestamp=datetime.datetime.now().isoformat(timespec="seconds"))
                append_row(csv_path, row)


if __name__ == "__main__":
    main()
