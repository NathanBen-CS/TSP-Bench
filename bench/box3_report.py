"""Summarise the 3D study, results/3d/runs.csv.

    py bench/box3_report.py        # writes results/3d/summary.md and summary.csv

Main measure: % longer than greedy TO Rn (greedy + 2-opt + Or-opt, rounds until
clean, K = 16) on the same map; negative = the setting beats the baseline.
Synthetic maps: mean +/- sample standard deviation over the 3 seeds.  Stars:
one map each, so no spread.  Also:
  vs best  stars: % above the published optimal or best tour (refs/star_bounds.csv)
  vs HK    synthetic up to 100k: % above OUR OWN Held-Karp 1-tree bound (refs/hk_own.csv)
Cleanup levels of one Rn run: none (start_length), T Rn (after_2opt_length), TO Rn.
"""

import csv
import re
import statistics
from collections import defaultdict

from common import ROOT

BASE = "greedy+2opt+oropt"
MERGES = ["M1", "M2", "M3", "M4", "M5"]
CURVES = ["hilbert3", "moore3", "peano3", "sierpinski3", "morton3", "snake3"]
STARS = ["star100", "star1k", "star10k", "kj37859", "hyg109399", "hyg119614", "star250k", "gaia1"]
SIZES = ["1k", "10k", "100k", "1M"]
KINDS = [("E3d", "uniform"), ("C3d", "clustered"), ("CT3d", "tight clustered")]
OUT = ROOT / "results" / "3d"


# ---------------------------------------------------------------- data

class Data:
    def __init__(self):
        self.runs = {}                   # (setting, map) -> dict
        self.problems = defaultdict(int)
        rows = [(r, True) for r in csv.DictReader(open(OUT / "runs.csv", newline=""))]
        if (OUT / "runs_untimed.csv").exists():      # lengths only (several runs at once)
            rows += [(r, False) for r in csv.DictReader(open(OUT / "runs_untimed.csv", newline=""))]
        for r, timed in rows:
            if (r["algorithm"], r["instance"]) in self.runs:
                continue
            if r["status"] != "OK":
                if r["status"] != "SKIPPED":
                    self.problems[(r["algorithm"], r["instance"], r["status"])] += 1
                continue
            ex = dict(kv.split("=", 1) for kv in r["extras"].split(";") if "=" in kv)
            L = {"full": float(r["tour_length"])}
            if "start_length" in ex:
                L["none"] = float(ex["start_length"])
            if "after_2opt_length" in ex:
                L["2opt"] = float(ex["after_2opt_length"])
            self.runs[(r["algorithm"], r["instance"])] = {
                "L": L, "ex": ex, "user": float(r["user_s"]) if timed else None, "wall": float(r["wall_s"]),
                "mb": float(r["peak_mb"]), "n": int(r["n"])}
        self.best = {r["instance"]: float(r["best_tour"]) for r in
                     csv.DictReader(open(ROOT / "refs" / "star_bounds.csv", newline=""))}
        self.hk = {r["instance"]: float(r["hk_own"]) for r in
                   csv.DictReader(open(ROOT / "refs" / "hk_own.csv", newline="")) if r["dims"] == "3"}
        self.maps = sorted({m for _, m in self.runs})

    def group(self, prefix, size):
        return [m for m in self.maps if re.fullmatch(rf"{prefix}{size}\.\d", m)]

    def length(self, name, m, lvl="full"):
        r = self.runs.get((name, m))
        return r["L"].get(lvl) if r else None

    def vs_base(self, name, maps, lvl="full"):
        """list of % longer than greedy TO Rn, one per map that has both"""
        out = []
        for m in maps:
            a, b = self.length(name, m, lvl), self.length(BASE, m)
            if a and b:
                out.append((a / b - 1) * 100)
        return out

    def vs_hk(self, name, maps, lvl="full"):
        return [(self.length(name, m, lvl) / self.hk[m] - 1) * 100 for m in maps
                if m in self.hk and self.length(name, m, lvl)]

    def vs_best(self, name, m, lvl="full"):
        a = self.length(name, m, lvl)
        return (a / self.best[m] - 1) * 100 if a and m in self.best else None

    def extra(self, name, maps, key):
        return [float(self.runs[(name, m)]["ex"][key]) for m in maps
                if (name, m) in self.runs and key in self.runs[(name, m)]["ex"]]

    def user(self, name, maps):
        return [u for m in maps if (name, m) in self.runs and (u := self.runs[(name, m)]["user"]) is not None]


def mean(v):
    return statistics.mean(v) if v else None


def sd(v):
    return statistics.stdev(v) if len(v) > 1 else None


def f(v, digits=2, plus=True):
    """mean (+/- sd) of a list, or a dash"""
    if v is None or v == []:
        return "–"
    if not isinstance(v, list):
        v = [v]
    s = f"{mean(v):{'+' if plus else ''}.{digits}f}"
    return s + (f" ± {sd(v):.{digits}f}" if len(v) > 1 else "")


def table(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    return "\n".join(out + ["| " + " | ".join(str(c) for c in r) + " |" for r in rows])


# ---------------------------------------------------------------- report

def baseline_tables(d, add):
    names = [("greedy", "greedy"), ("twoopt", "greedy T R1"), ("twoopt_or", "greedy TO R1"),
             ("greedy+2opt", "greedy T Rn"), (BASE, "greedy TO Rn")]
    add("## 1. Baseline (greedy, T, TO, R1, Rn)", "")
    add("**Stars: % above the published optimal or best tour** (best = Waterloo, "
        "refs/star_bounds.csv)", "")
    rows = [[m, d.runs[(BASE, m)]["n"]] + [f(d.vs_best(a, m)) for a, _ in names]
            + [f(d.user(BASE, [m]), plus=False)] for m in STARS if (BASE, m) in d.runs]
    add(table(["map", "cities"] + [x for _, x in names] + ["TO Rn user s"], rows), "")
    add("**Synthetic: % above OUR OWN Held-Karp bound** (mean ± sd over seeds; no bound at 1M)", "")
    rows = []
    for p, kind in KINDS:
        for s in SIZES[:3]:
            g = d.group(p, s)
            if g and any(m in d.hk for m in g):
                rows.append([kind, s] + [f(d.vs_hk(a, g)) for a, _ in names])
    add(table(["maps", "size"] + [x for _, x in names], rows), "")
    add("**Synthetic: absolute lengths (mean over seeds, millions) and user seconds**", "")
    rows = []
    for p, kind in KINDS:
        for s in SIZES:
            g = [m for m in d.group(p, s) if (BASE, m) in d.runs]
            if g:
                rows.append([kind, s] + [f"{mean([d.length(a, m) for m in g]) / 1e6:.3f} ({mean(d.user(a, g)):.2f} s)"
                                         for a, _ in names])
    add(table(["maps", "size"] + [x for _, x in names], rows), "")


def main():
    d = Data()
    lines = []
    add = lambda *x: lines.extend(x)
    add("# Box-TSP in 3D: results", "",
        f"Generated by `bench/box3_report.py` from `results/3d/runs.csv` ({len(d.runs)} OK runs). "
        "\"vs greedy TO Rn\" = % longer than the baseline on the same map (negative = beats it), "
        "mean ± sample sd over seeds.", "")
    if d.problems:
        add("**Runs not OK:** " + ", ".join(f"{a} on {m}: {s}" for (a, m, s) in d.problems), "")
    baseline_tables(d, add)

    groups = [(f"{kind}, {s}", d.group(p, s)) for p, kind in KINDS for s in SIZES]
    groups += [(m, [m]) for m in STARS]
    groups = [(n, g) for n, g in groups if any((BASE, m) in d.runs for m in g)]

    # ---- curve x merge, all four cleanup variants
    add("## 2. Curve × merge at t = 10: % vs greedy TO Rn", "")
    variants = [("TO Rn", "", "full"), ("T Rn", "", "2opt"), ("TO R1", ":TO R1", "full"),
                ("T R1", ":T R1", "full"), ("no cleanup", "", "none")]
    for gname, g in groups:
        if not any((f"box:{c}:{m}:t10", x) in d.runs for c in CURVES for m in MERGES for x in g):
            continue
        add(f"### {gname} ({len(g)} map{'s' if len(g) > 1 else ''})", "")
        for label, suffix, lvl in variants:
            rows = []
            for c in CURVES:
                cells = [d.vs_base(f"box:{c}:{m}:t10{suffix}", g, lvl) for m in MERGES]
                if any(cells):
                    rows.append([c] + [f(x) for x in cells])
            if rows:
                add(f"**{label}**", "", table(["curve"] + MERGES, rows), "")

    # ---- best box setting against the baseline, per group
    add("## 3. Does any box setting beat greedy TO Rn?", "")
    rows = []
    every = sorted({a for a, _ in d.runs if a.startswith(("box:", "kd:"))})
    for gname, g in groups:
        full = [(mean(v), v, a) for a in every
                if len(v := d.vs_base(a, g)) == len([m for m in g if (BASE, m) in d.runs]) and v]
        if not full:
            continue
        full.sort(key=lambda x: x[0])
        beat = sum(x[0] < 0 for x in full)
        wins = sum(v < 0 for x in full for v in x[1])
        runs = sum(len(x[1]) for x in full)
        rows.append([gname, len(full), f"{full[0][2]}: {f(full[0][1])}", beat, f"{wins} of {runs}",
                     f(d.vs_base("twoopt_or", g))])
    add(table(["maps", "box settings run", "best (mean)", "settings with mean < 0", "single runs < 0",
               "greedy TO R1"], rows), "")

    # ---- t sweep
    add("## 4. Box size t: % vs greedy TO Rn, TO Rn (before cleanup in brackets)", "")
    ts = list(range(4, 13))
    for gname, g in groups:
        rows = []
        for c in CURVES:
            for m in MERGES:
                if sum((f"box:{c}:{m}:t{t}", x) in d.runs for t in ts for x in g) > 2 * len(g):
                    rows.append([f"{c} {m}"] + [f"{f(mean(d.vs_base(f'box:{c}:{m}:t{t}', g)))} "
                                                f"({f(mean(d.vs_base(f'box:{c}:{m}:t{t}', g, 'none')), 0)})"
                                                for t in ts])
        if rows:
            add(f"**{gname}**", "", table(["curve, merge"] + [f"t={t}" for t in ts], rows), "")

    # ---- K sweep
    add("## 5. Neighbour-list size K (TO Rn)", "",
        "% vs greedy TO Rn at K = 16 on the same map; user seconds in brackets.", "")
    ks = [8, 16, 24, 32]
    kname = lambda base, k: base if k == 16 else base + f":k{k}"
    for gname, g in groups:
        rows = []
        for label, base in [("greedy TO Rn", None)] + [(f"{c} {m}", f"box:{c}:{m}:t10") for c in CURVES for m in MERGES]:
            names = [(BASE if k == 16 else f"greedy TO Rn:k{k}") if base is None else kname(base, k) for k in ks]
            if sum(any((n, x) in d.runs for x in g) for n in names) < 3:
                continue
            rows.append([label] + [f"{f(d.vs_base(n, g))} ({f(mean(d.user(n, g)), plus=False)} s)" for n in names])
        if rows:
            add(f"**{gname}**", "", table(["setting"] + [f"K={k}" for k in ks], rows), "")

    # ---- kd boxes
    add("## 6. kd-median boxes (hold exactly t cities) against curve cells, hilbert3: % vs greedy TO Rn", "")
    for gname, g in groups:
        rows = []
        for m in MERGES:
            for t in (4, 6, 8, 10, 12):
                a, b = d.vs_base(f"kd:hilbert3:{m}:t{t}", g), d.vs_base(f"box:hilbert3:{m}:t{t}", g)
                if a:
                    rows.append([m, t, f(d.vs_base(f"kd:hilbert3:{m}:t{t}", g, "none")), f(a),
                                 f(d.vs_base(f"box:hilbert3:{m}:t{t}", g, "none")), f(b)])
        if rows:
            add(f"**{gname}**", "", table(["merge", "t", "kd, no cleanup", "kd, TO Rn", "cells, no cleanup",
                                           "cells, TO Rn"], rows), "")

    # ---- plain orders
    add("## 7. Plain curve orders and nn: % vs greedy TO Rn", "")
    names = [f"curve:{c}" for c in CURVES] + ["nn"]
    rows = [[gname] + [f(d.vs_base(a, g), 1) for a in names] for gname, g in groups]
    add(table(["maps"] + names, rows), "")

    # ---- time and memory
    add("## 8. User CPU seconds (mean) at t = 10, TO Rn, and peak memory", "")
    rows = []
    for gname, g in groups:
        g = [m for m in g if (BASE, m) in d.runs]
        cells = [f"{mean(d.user(BASE, g)):.2f}", f"{mean(d.user('twoopt_or', g)):.2f}"]
        for m in MERGES:
            u = d.user(f"box:hilbert3:{m}:t10", g)
            cells.append(f"{mean(u):.2f}" if u else "–")
        mb = lambda a: max([d.runs[(a, m)]["mb"] for m in g if (a, m) in d.runs], default=0)
        cells += [f"{mb(BASE):.0f}", f"{mb('box:hilbert3:M1:t10'):.0f}", f"{mb('box:hilbert3:M5:t10'):.0f}"]
        rows.append([gname] + cells)
    add(table(["maps", "greedy TO Rn", "greedy TO R1"] + [f"hilbert3 {m}" for m in MERGES]
              + ["MB greedy", "MB M1", "MB M5"], rows), "")

    (OUT / "summary.md").write_text("\n".join(lines), encoding="utf-8")

    # ---- summary.csv: one row per setting and group
    out = []
    for gname, g in groups:
        for a in sorted({a for a, m in d.runs if m in g}):
            row = {"maps": gname, "setting": a, "n_maps": len([m for m in g if (a, m) in d.runs])}
            for lvl in ("none", "2opt", "full"):
                v = d.vs_base(a, g, lvl)
                row[f"vs_greedyTORn_{lvl}"] = f"{mean(v):.4f}" if v else ""
                row[f"sd_{lvl}"] = f"{sd(v):.4f}" if sd(v) is not None else ""
            h = d.vs_hk(a, g)
            row["vs_own_hk"] = f"{mean(h):.4f}" if h else ""
            b = [x for m in g if (x := d.vs_best(a, m)) is not None]
            row["vs_best"] = f"{mean(b):.4f}" if b else ""
            u = d.user(a, g)
            row["user_s"] = f"{mean(u):.3f}" if u else ""
            for k in ("sort_time", "split_time", "dp_time", "glue_time", "neighbor_time", "greedy_time",
                      "twoopt_time", "oropt_time", "boxes", "capped_boxes", "glue_share", "join_share"):
                v = d.extra(a, g, k)
                row[k] = f"{mean(v):.4f}" if v else ""
            out.append(row)
    with open(OUT / "summary.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)
    print(f"Wrote {OUT / 'summary.md'} and summary.csv ({len(d.runs)} runs, {len(d.problems)} problems)")


if __name__ == "__main__":
    main()
