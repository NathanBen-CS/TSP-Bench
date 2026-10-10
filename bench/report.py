"""Turn results/runs.csv into the report.

    gap over HK  = (tour length - Held-Karp bound) / Held-Karp bound x 100
    gap over opt = (tour length - optimal) / optimal x 100

Both references come from refs/bounds (DIMACS "bounds" file: name, HK bound,
and optimal / best-known length where available).

Writes results/report.md (one table per algorithm + machine benchmark) and
results/report_full.csv (every run with gaps and extras as columns).

    python bench/report.py
"""

import argparse
import csv
import json
from collections import defaultdict

from common import ROOT, load_bounds, load_config, resolve


def gap(length, ref):
    return f"{(length - ref) / ref * 100:.2f}" if ref else ""


def md_table(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    bounds = load_bounds(resolve(cfg["bounds_file"]))
    with open(resolve(cfg["results_csv"]), newline="") as f:
        runs = list(csv.DictReader(f))

    by_alg = defaultdict(list)
    for r in runs:
        r["extras"] = dict(kv.split("=", 1) for kv in r["extras"].split(";") if "=" in kv)
        hk, opt = bounds.get(r["instance"], (None, None))
        ok = r["status"] == "OK"
        r["gap_hk"] = gap(int(r["tour_length"]), hk) if ok else ""
        r["gap_opt"] = gap(int(r["tour_length"]), opt) if ok else ""
        by_alg[r["algorithm"]].append(r)

    lines = ["# TSP benchmark report", ""]
    machine = ROOT / "results" / "machine.json"
    if machine.exists():
        info = json.loads(machine.read_text())
        lines += ["## Machine", "", md_table(["Field", "Value"], list(info.items())), ""]

    full_rows = []
    for alg, rows in by_alg.items():
        rows.sort(key=lambda r: (int(r["n"]), r["instance"], r["seed"]))
        extra_keys = sorted({k for r in rows for k in r["extras"]})
        header = ["Instance", "N", "Seed", "Tour length", "User time (s)",
                  "Peak memory (MB)", "% over HK", "% over opt", "Status"] + extra_keys
        table = [[r["instance"], r["n"], r["seed"], r["tour_length"], r["user_s"],
                  r["peak_mb"], r["gap_hk"], r["gap_opt"], r["status"]]
                 + [r["extras"].get(k, "") for k in extra_keys] for r in rows]
        lines += [f"## {alg}", "", md_table(header, table), ""]
        for r in rows:
            full_rows.append({**{k: v for k, v in r.items() if k != "extras"}, **r["extras"]})

    bench = ROOT / "results" / "machine_bench.csv"
    if bench.exists():
        with open(bench, newline="") as f:
            b = list(csv.DictReader(f))
        lines += ["## Machine benchmark (DIMACS greedy, pinned to one CPU)", "",
                  md_table(["Instance", "N", "Iterations", "Run", "User time (s)",
                            "User s / iteration"],
                           [[r["instance"], r["n"], r["iterations"], r["run"], r["user_s"],
                             r["user_s_per_iteration"]] for r in b]), ""]

    lines += ["## Notes", "",
              "- User time and peak memory cover the whole program, including reading "
              "the instance and building data structures.",
              "- Peak memory is the Windows peak working set (equivalent of Linux max RSS).",
              "- Windows counts process CPU time in ~15.6 ms ticks, so times below about "
              "0.05 s are imprecise (0.000 means less than one tick).",
              "- Randomized algorithms: one row per seed, no averaging.",
              "- `codes/` contains the DIMACS programs with Windows portability fixes only "
              "(ctype.h include, NUL instead of /dev/null, unused sys/times.h removed, "
              "binary-mode file reading with CR stripping). Generated instances match "
              "the DIMACS samples byte for byte apart from line endings.",
              ""]

    out_dir = resolve(cfg["results_csv"]).parent      # report sits next to runs.csv
    out = out_dir / "report.md"
    out.write_text("\n".join(lines))
    keys = list(dict.fromkeys(k for r in full_rows for k in r))
    full = out_dir / "report_full.csv"
    with open(full, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(full_rows)
    print(f"Wrote {out} and {full} ({len(runs)} runs)")


if __name__ == "__main__":
    main()
