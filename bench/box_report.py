"""Summarise the box-TSP experiment (Box-TSP plan, Task 4).

Reads the raw runs (results/2d/box/runs.csv, written by bench/run.py with
bench/box_config.json) and writes

    results/2d/box/summary.md    tables: mean and spread (sample std dev) over seeds
    results/2d/box/summary.csv   the same numbers, one row per algorithm/type/size

References for the gap columns:
  uniform    0.7124 * sqrt(N * A), A = 10^12 (the generator's square).  This
             estimate of the optimal tour is only valid for uniform data.
  clustered  the Held-Karp lower bound from refs/bounds where it exists.
             New seeds (C100k.2, C1M.*) have no bound, so the clustered gap
             at those sizes covers fewer seeds; the "n" column says how many.
"vs base" compares each algorithm with greedy+2opt+oropt on the same
instance: (length / baseline length - 1) * 100, then mean and spread over
seeds.  It needs no reference value, so it covers every seed.

    py bench/box_report.py
"""

import argparse
import csv
import math
import statistics
from collections import defaultdict

from common import ROOT, load_bounds, load_config, resolve

BASE = "greedy+2opt+oropt"
FULL = "boxes+2opt+oropt"
AREA = 1e12
SIZE_ORDER = {"1k": 1, "10k": 2, "100k": 3, "1M": 4}
STAGES = ["read_time", "sort_time", "split_time", "dp_time", "glue_time",
          "neighbor_time", "twoopt_time", "oropt_time"]


def mean_sd(values):
    values = [v for v in values if v is not None]
    if not values:
        return None, None, 0
    sd = statistics.stdev(values) if len(values) > 1 else 0.0
    return statistics.mean(values), sd, len(values)


def fmt(m, sd, digits=2, signed=True):
    if m is None:
        return "–"
    s = f"{m:+.{digits}f}" if signed else f"{m:.{digits}f}"
    return f"{s} ± {sd:.{digits}f}"


def md_table(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=str(ROOT / "bench" / "box_config.json"))
    args = ap.parse_args()

    cfg = load_config(args.config)
    bounds = load_bounds(resolve(cfg["bounds_file"]))
    runs_path = resolve(cfg["results_csv"])
    with open(runs_path, newline="") as f:
        rows = list(csv.DictReader(f))
    algs = [a["name"] for a in cfg["algorithms"]]

    # ---- index the OK runs by (algorithm, instance) and note problems
    ok, problems = {}, defaultdict(int)
    for r in rows:
        if r["status"] == "OK":
            r["extras"] = dict(kv.split("=", 1) for kv in r["extras"].split(";") if "=" in kv)
            ok[(r["algorithm"], r["instance"])] = r
        elif r["status"] != "SKIPPED":
            problems[(r["algorithm"], r["status"])] += 1

    def group_of(inst):                      # "E100k.2" -> ("uniform", "100k")
        return ("uniform" if inst[0] == "E" else "clustered", inst[1:].split(".")[0])

    def reference(inst, n):
        kind, _ = group_of(inst)
        if kind == "uniform":
            return 0.7124 * math.sqrt(n * AREA)
        hk = bounds.get(inst, (None, None))[0]
        return hk

    instances = sorted({r["instance"] for r in rows},
                       key=lambda i: (i[0], SIZE_ORDER.get(group_of(i)[1], 9), i))
    n_of = {r["instance"]: int(r["n"]) for r in rows}
    groups = defaultdict(list)                # (type, size) -> instances
    for inst in instances:
        groups[group_of(inst)].append(inst)

    # ---- per (alg, type, size): gap, vs base, time, box extras
    stats = {}
    for (kind, size), insts in groups.items():
        for alg in algs:
            gaps, vsb, times, starts, after2, glue = [], [], [], [], [], []
            stage = defaultdict(list)
            for inst in insts:
                r = ok.get((alg, inst))
                if not r:
                    continue
                L, n = int(r["tour_length"]), n_of[inst]
                ref = reference(inst, n)
                gaps.append((L / ref - 1) * 100 if ref else None)
                base = ok.get((BASE, inst))
                vsb.append((L / int(base["tour_length"]) - 1) * 100 if base else None)
                times.append(float(r["user_s"]))
                ex = r["extras"]
                if "start_length" in ex and ref:
                    starts.append((float(ex["start_length"]) / ref - 1) * 100)
                if "after_2opt_length" in ex and ref:
                    after2.append((float(ex["after_2opt_length"]) / ref - 1) * 100)
                if "glue_share" in ex:
                    glue.append(float(ex["glue_share"]))
                for s in STAGES:
                    if s in ex:
                        stage[s].append(float(ex[s]))
            if not times:
                continue
            stats[(alg, kind, size)] = {
                "gap": mean_sd(gaps), "vs_base": mean_sd(vsb), "user": mean_sd(times),
                "start": mean_sd(starts), "after2": mean_sd(after2), "glue": mean_sd(glue),
                "stages": {s: mean_sd(v)[0] for s, v in stage.items()},
                "seeds": len(times)}

    sizes = sorted({s for (_, s) in groups}, key=lambda s: SIZE_ORDER.get(s, 9))
    kinds = ["uniform", "clustered"]
    ref_note = {"uniform": "% over 0.7124·√(N·A)", "clustered": "% over Held-Karp bound"}

    lines = ["# Box-TSP experiment: results", "",
             f"Raw runs: `{runs_path.relative_to(ROOT)}` · mean ± sample standard deviation "
             "over seeds · user CPU time of the whole program (file reading included).", ""]
    if problems:
        lines += ["**Runs that did not finish OK:** " + ", ".join(
            f"{a} {s} ×{c}" for (a, s), c in sorted(problems.items())), ""]

    # 1. tour quality
    for kind in kinds:
        lines += [f"## Tour quality, {kind} ({ref_note[kind]})", ""]
        header = ["algorithm"] + [f"{s} (n)" for s in sizes]
        body = []
        for alg in algs:
            cells = []
            for s in sizes:
                st = stats.get((alg, kind, s))
                cells.append("–" if not st else f"{fmt(*st['gap'][:2])} ({st['gap'][2]})")
            body.append([alg] + cells)
        lines += [md_table(header, body), ""]

    # 2. paired comparison with the greedy baseline
    for kind in kinds:
        lines += [f"## Versus {BASE}, {kind} (% longer than the baseline on the same instance; "
                  "negative = better)", ""]
        header = ["algorithm"] + sizes
        body = [[alg] + [fmt(*stats[(alg, kind, s)]["vs_base"][:2])
                         if (alg, kind, s) in stats else "–" for s in sizes] for alg in algs]
        lines += [md_table(header, body), ""]

    # 3. time
    for kind in kinds:
        lines += [f"## User time in seconds, {kind}", ""]
        header = ["algorithm"] + sizes
        body = [[alg] + [fmt(*stats[(alg, kind, s)]["user"][:2], signed=False)
                         if (alg, kind, s) in stats else "–" for s in sizes] for alg in algs]
        lines += [md_table(header, body), ""]

    # 4. box size t: before and after cleanup
    t_algs = [("6", f"{FULL}-t6"), ("8", f"{FULL}-t8"), ("10", FULL), ("12", f"{FULL}-t12")]
    for kind in kinds:
        lines += [f"## Box size t, {kind}: gap of the box tour before cleanup → after 2-opt "
                  "→ after Or-opt", ""]
        header = ["t"] + sizes
        body = []
        for t, alg in t_algs:
            cells = []
            for s in sizes:
                st = stats.get((alg, kind, s))
                if not st or st["start"][0] is None:
                    cells.append("–")
                else:
                    cells.append(f"{st['start'][0]:+.2f} → {st['after2'][0]:+.2f} → "
                                 f"{st['gap'][0]:+.2f}")
            body.append([t] + cells)
        lines += [md_table(header, body), ""]

    # 5. stage times and glue share of the full method
    for kind in kinds:
        lines += [f"## Full method ({FULL}, t = 10): time per stage (s) and glue share, {kind}",
                  ""]
        header = ["size"] + [s.replace("_time", "") for s in STAGES] + ["G / box tour"]
        body = []
        for s in sizes:
            st = stats.get((FULL, kind, s))
            if not st:
                continue
            body.append([s] + [f"{st['stages'].get(x, 0):.2f}" for x in STAGES]
                        + [fmt(*st["glue"][:2], digits=3, signed=False)])
        lines += [md_table(header, body), ""]

    out_dir = runs_path.parent
    (out_dir / "summary.md").write_text("\n".join(lines), encoding="utf-8")
    with open(out_dir / "summary.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["algorithm", "type", "size", "seeds", "gap_mean", "gap_sd", "gap_n",
                    "vs_base_mean", "vs_base_sd", "user_s_mean", "user_s_sd",
                    "start_gap_mean", "after_2opt_gap_mean", "glue_share_mean"]
                   + [f"{s}_mean" for s in STAGES])
        for (alg, kind, s), st in sorted(stats.items(),
                                         key=lambda k: (algs.index(k[0][0]), k[0][1],
                                                        SIZE_ORDER.get(k[0][2], 9))):
            w.writerow([alg, kind, s, st["seeds"], *st["gap"], *st["vs_base"][:2],
                        *st["user"][:2], st["start"][0], st["after2"][0], st["glue"][0]]
                       + [st["stages"].get(x) for x in STAGES])
    print(f"Wrote {out_dir / 'summary.md'} and {out_dir / 'summary.csv'}")


if __name__ == "__main__":
    main()
