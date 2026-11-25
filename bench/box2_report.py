"""Summarise the 2D box extension (curves x merges x t), results/2d/box2/runs.csv.

    py bench/box2_report.py            # writes results/2d/box2/summary.md and .csv

For every setting and cleanup level:
  * paired difference from the reference setting (Hilbert, M1, t = 10) on
    the same instance: (length / reference length - 1) x 100, then mean and
    sample standard deviation over seeds.  Negative = shorter than reference.
  * % above the Held-Karp bound (refs/bounds), where the instance has one.
Cleanup levels come from one run each: start_length (none),
after_2opt_length (2-opt) and the final tour (2-opt + Or-opt).
"""

import csv
import statistics
from collections import defaultdict

from common import ROOT, load_bounds

REF = "box:hilbert:M1:t10"
LEVELS = [("none", "start_length"), ("2opt", "after_2opt_length"), ("full", None)]
CURVES = ["hilbert", "moore", "sierpinski", "knopp", "peano", "morton", "strip"]
MERGES = ["M1", "M2", "M3", "M4", "M5"]
TS = [4, 6, 8, 10, 12]
STAGES = ["sort_time", "split_time", "dp_time", "glue_time", "neighbor_time",
          "twoopt_time", "oropt_time"]


def ms(values):
    values = [v for v in values if v is not None]
    if not values:
        return None, None, 0
    return statistics.mean(values), (statistics.stdev(values) if len(values) > 1 else 0.0), len(values)


def f(m, sd=None, digits=2):
    if m is None:
        return "–"
    return f"{m:+.{digits}f}" + (f" ± {sd:.{digits}f}" if sd is not None else "")


def table(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    return "\n".join(out + ["| " + " | ".join(str(c) for c in r) + " |" for r in rows])


def main():
    bounds = load_bounds(ROOT / "refs" / "bounds")
    path = ROOT / "results" / "2d" / "box2" / "runs.csv"
    rows = list(csv.DictReader(open(path, newline="")))
    problems = defaultdict(int)
    runs = {}                                       # (name, instance) -> lengths/extras
    for r in rows:
        if r["status"] != "OK":
            if r["status"] != "SKIPPED":
                problems[(r["algorithm"], r["status"])] += 1
            continue
        ex = dict(kv.split("=", 1) for kv in r["extras"].split(";") if "=" in kv)
        L = {"full": float(r["tour_length"])}
        for lvl, key in LEVELS:
            if key and key in ex:
                L[lvl] = float(ex[key])
        runs[(r["algorithm"], r["instance"])] = {"L": L, "ex": ex, "user": float(r["user_s"]),
                                                  "n": int(r["n"])}

    instances = sorted({i for _, i in runs})
    kind = lambda i: "uniform" if i[0] == "E" else "clustered"
    size = lambda i: i[1:].split(".")[0]
    groups = defaultdict(list)                      # (kind, size) -> instances
    for i in instances:
        groups[(kind(i), size(i))].append(i)

    def paired(name, insts, lvl):
        vals = []
        for i in insts:
            a, b = runs.get((name, i)), runs.get((REF, i))
            if a and b and lvl in a["L"] and lvl in b["L"]:
                vals.append((a["L"][lvl] / b["L"][lvl] - 1) * 100)
        return ms(vals)

    def gap(name, insts, lvl):
        vals = []
        for i in insts:
            a = runs.get((name, i))
            if a and lvl in a["L"] and i in bounds:
                vals.append((a["L"][lvl] / bounds[i][0] - 1) * 100)
        return ms(vals)

    def extra(name, insts, key, rel=None):
        vals = []
        for i in insts:
            a = runs.get((name, i))
            if a and key in a["ex"]:
                v = float(a["ex"][key])
                vals.append(v / a["L"][rel] if rel else v)
        return ms(vals)

    names = sorted({n for n, _ in runs})
    lines = ["# Box-TSP 2D extension: results", "",
             f"Raw runs: `results/2d/box2/runs.csv` ({len(rows)} rows). Reference setting: "
             f"`{REF}`. \"vs ref\" = % longer than the reference on the same instance "
             "(negative = shorter), mean ± sample std dev over seeds. \"vs HK\" = % above "
             "the Held-Karp bound. Cleanup levels: none (box tour as merged), 2-opt, "
             "full (2-opt then 2-opt + Or-opt).", ""]
    if problems:
        lines += ["**Runs not OK:** " + ", ".join(f"{a} {s} ×{c}" for (a, s), c in problems.items()), ""]

    csv_rows = []
    for (k, sz), insts in sorted(groups.items(), key=lambda g: (g[0][1], g[0][0])):
        lines += [f"## {k}, {sz} cities ({len(insts)} seeds)", ""]

        # --- screen: curve x merge at t = 10, per cleanup level
        if any((f"box:{c}:{m}:t10", i) in runs for c in CURVES for m in MERGES for i in insts):
            for lvl, _ in LEVELS:
                lines += [f"**Curve × merge at t = 10, cleanup = {lvl}: vs ref (%)**", ""]
                body = []
                for c in CURVES:
                    cells = []
                    for m in MERGES:
                        mm, sd, n = paired(f"box:{c}:{m}:t10", insts, lvl)
                        cells.append(f(mm, sd))
                    body.append([c] + cells)
                lines += [table(["curve"] + MERGES, body), ""]

        # --- every setting: gaps, join share, time
        body = []
        for name in names:
            if not any((name, i) in runs for i in insts):
                continue
            row = [name]
            for lvl, _ in LEVELS:
                mm, sd, n = paired(name, insts, lvl)
                row.append(f(mm, sd))
            g = gap(name, insts, "full")
            row.append(f(g[0]) + (f" ({g[2]})" if g[2] else ""))
            j = extra(name, insts, "join_cost", "none")
            row.append("–" if j[0] is None else f"{j[0]:.3f}")
            t = ms([runs[(name, i)]["user"] for i in insts if (name, i) in runs])
            row.append(f"{t[0]:.2f}")
            body.append(row)
            csv_rows.append({"setting": name, "type": k, "size": sz,
                             **{f"vs_ref_{lvl}": paired(name, insts, lvl)[0] for lvl, _ in LEVELS},
                             **{f"vs_ref_{lvl}_sd": paired(name, insts, lvl)[1] for lvl, _ in LEVELS},
                             **{f"vs_hk_{lvl}": gap(name, insts, lvl)[0] for lvl, _ in LEVELS},
                             "join_share": j[0], "user_s": t[0],
                             **{s: extra(name, insts, s)[0] for s in STAGES}})
        lines += ["**Every setting**", "",
                  table(["setting", "vs ref: none", "vs ref: 2-opt", "vs ref: full",
                         "vs HK: full (n)", "join cost / box tour", "user s"], body), ""]

        # --- t sweep (stage B)
        combos = sorted({tuple(n.split(":")[1:3]) for n in names if n.startswith("box:")
                         and sum((f"box:{n.split(':')[1]}:{n.split(':')[2]}:t{t}", i) in runs
                                 for t in TS for i in insts) > len(insts)})
        if combos:
            lines += ["**Box size t: vs ref at cleanup none → full**", ""]
            body = []
            for c, m in combos:
                cells = []
                for t in TS:
                    a = paired(f"box:{c}:{m}:t{t}", insts, "none")[0]
                    b = paired(f"box:{c}:{m}:t{t}", insts, "full")[0]
                    cells.append("–" if a is None else f"{a:+.2f} → {b:+.2f}")
                body.append([f"{c} {m}"] + cells)
            lines += [table(["curve, merge"] + [f"t = {t}" for t in TS], body), ""]

    # --- ranking for stage selection (10k, full cleanup, both data types)
    tenk = groups.get(("uniform", "10k"), []) + groups.get(("clustered", "10k"), [])
    if tenk:
        lines += ["## Ranking at 10k (full cleanup, uniform and clustered pooled)", ""]
        cur = sorted(((ms([paired(f"box:{c}:{m}:t10", tenk, "full")[0] for m in MERGES])[0], c)
                      for c in CURVES if (f"box:{c}:M1:t10", tenk[0]) in runs), key=lambda x: x[0])
        mer = sorted(((ms([paired(f"box:{c}:{m}:t10", tenk, "full")[0] for c in CURVES])[0], m)
                      for m in MERGES if (f"box:hilbert:{m}:t10", tenk[0]) in runs), key=lambda x: x[0])
        lines += ["Curves (mean over merges): " + ", ".join(f"{c} {v:+.2f}" for v, c in cur), "",
                  "Merges (mean over curves): " + ", ".join(f"{m} {v:+.2f}" for v, m in mer), ""]
        print("curves:", [c for _, c in cur])
        print("merges:", [m for _, m in mer])

    out = ROOT / "results" / "2d" / "box2" / "summary.md"
    out.write_text("\n".join(lines), encoding="utf-8")
    with open(ROOT / "results" / "2d" / "box2" / "summary.csv", "w", newline="") as fh:
        if csv_rows:
            w = csv.DictWriter(fh, fieldnames=list(csv_rows[0]))
            w.writeheader()
            w.writerows(csv_rows)
    print(f"Wrote {out} and summary.csv")


if __name__ == "__main__":
    main()
