"""Summarise the hull-ring study, results/2d/hull/runs.csv.

    py bench/hull_report.py      # writes results/2d/hull/summary.md, summary.csv, hull_clean_vs_raw.csv

Measures:
  vs HK    % above the published DIMACS Held-Karp bound (refs/bounds); only maps that
           have one: 1k and 10k .0-.2, 100k .0-.1, E1M.0 (none for C1M.0 or the shapes)
  vs base  % longer than greedy TO Rn (greedy + 2-opt + Or-opt, rounds until clean) on
           the same map; every map
Mean +/- sample standard deviation over the maps of a group.  Cleanup levels of one
Rn run: raw (start_length, no cleanup), T Rn (after_2opt_length), TO Rn (the tour).
"""

import csv
import math
import re
import statistics
from collections import defaultdict

from common import ROOT, load_bounds

OUT = ROOT / "results" / "2d" / "hull"
BASE = "greedy+2opt+oropt"
CURVES = ["hilbert", "moore", "sierpinski", "knopp", "peano", "morton", "strip"]
MERGES = ["M1", "M2", "M3", "M4", "M5"]
HULL = ["H0"] + MERGES
SIZES = ["1k", "10k", "100k", "1M"]
KINDS = [("E", "uniform"), ("C", "clustered")]
SHAPES = [("ann", "annulus"), ("circ", "concentric circles"), ("blob", "Gaussian blob"),
          ("rect", "thin rectangle 10:1")]


def mean(v):
    return statistics.mean(v) if v else None


def sd(v):
    return statistics.stdev(v) if len(v) > 1 else None


def f(v, dg=2):
    """mean +/- sd of a list, or '-'"""
    if not v:
        return "–"
    s = f"{mean(v):+.{dg}f}"
    return s + (f" ± {sd(v):.{dg}f}" if len(v) > 1 else "")


class Data:
    def __init__(self):
        self.runs = {}
        self.problems = defaultdict(int)
        for r in csv.DictReader(open(OUT / "runs.csv", newline="")):
            key = (r["algorithm"], r["instance"])
            if key in self.runs:
                continue
            if r["status"] != "OK":
                self.problems[(r["algorithm"], r["instance"], r["status"])] += 1
                continue
            ex = dict(kv.split("=", 1) for kv in r["extras"].split(";") if "=" in kv)
            L = {"full": float(r["tour_length"])}
            if "start_length" in ex:
                L["none"] = float(ex["start_length"])
            if "after_2opt_length" in ex:
                L["2opt"] = float(ex["after_2opt_length"])
            self.runs[key] = {"L": L, "ex": ex, "user": float(r["user_s"]), "n": int(r["n"]),
                              "stamp": r["timestamp"]}
        self.hk = {k: v[0] for k, v in load_bounds(ROOT / "refs" / "bounds").items()}
        self.maps = sorted({m for _, m in self.runs})

    def group(self, prefix, size):
        return [m for m in self.maps if re.fullmatch(rf"{prefix}{size}\.\d", m)]

    def length(self, name, m, lvl="full"):
        r = self.runs.get((name, m))
        return r["L"].get(lvl) if r else None

    def vs_hk(self, name, maps, lvl="full"):
        return [(self.length(name, m, lvl) / self.hk[m] - 1) * 100 for m in maps
                if m in self.hk and self.length(name, m, lvl)]

    def vs_base(self, name, maps, lvl="full", base=BASE, base_lvl="full"):
        out = []
        for m in maps:
            a, b = self.length(name, m, lvl), self.length(base, m, base_lvl)
            if a and b:
                out.append((a / b - 1) * 100)
        return out

    def ex(self, name, maps, key):
        return [float(self.runs[(name, m)]["ex"][key]) for m in maps
                if (name, m) in self.runs and key in self.runs[(name, m)]["ex"]]

    def user(self, name, maps):
        return [self.runs[(name, m)]["user"] for m in maps if (name, m) in self.runs]

    def full(self, name, maps):
        return all((name, m) in self.runs for m in maps)

    def best(self, names, maps, score):
        """the name with the lowest mean score over maps (only names run on every map)"""
        best = None
        for n in names:
            if not self.full(n, maps):
                continue
            v = score(n)
            if v and (best is None or mean(v) < mean(best[1])):
                best = (n, v)
        return best


def table(head, rows):
    out = ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


def label(name):
    if name.startswith("hull:"):
        return "hull " + name[5:]
    if name.startswith("box:"):
        c, m = name.split(":")[1:3]
        return f"{c} {m}"
    return {BASE: "greedy TO Rn", "curve:hull": "hull order", "curve:hilbert": "Hilbert order"}.get(name, name)


def chart_lines(d, prefix, size):
    """the lines of the clean-vs-raw figure for one map group: (line, setting)"""
    ms = [m for m in d.group(prefix, size) if m in d.hk]
    lines = []
    for c in CURVES:
        names = [f"box:{c}:{m}:t10" for m in MERGES]
        if not all(d.full(n, ms) for n in names):     # best of M1-M5 only when all 5 ran
            continue                                   #   (strip has M1 only at 1M)
        b = d.best(names, ms, lambda n: d.vs_hk(n, ms))
        if b:
            lines.append((c, b[0]))
    b = d.best([f"hull:{m}" for m in HULL], ms, lambda n: d.vs_hk(n, ms))
    if b:
        lines.append(("hull", b[0]))
    lines.append(("greedy", BASE))
    lines.append(("Hilbert order", "curve:hilbert"))
    lines.append(("hull order", "curve:hull"))
    return ms, lines


def main():
    d = Data()
    md = ["# Hull rings: summary", "",
          "Made by `bench/hull_report.py` from `results/2d/hull/runs.csv`. "
          "Mean ± sample sd over the maps of a group. **vs HK** = % above the published "
          "Held-Karp bound (maps with one only); **vs greedy TO Rn** = % longer than greedy + 2-opt + "
          "Or-opt (Rn) on the same map. Raw = the start tour before any cleanup.", ""]
    csv_rows = []
    groups = [(f"{kind} {s}", p, s) for p, kind in KINDS for s in SIZES]

    # ---- 1. rings
    rows = []
    for g, p, s in groups + [(name, fam, "10k") for fam, name in SHAPES]:
        ms = d.group(p, s)
        if not d.full("hull:H0", ms) or not ms:
            continue
        n = d.runs[("hull:H0", ms[0])]["n"]
        rings = d.ex("hull:H0", ms, "rings")
        ratio = [float(e) / float(nn) for e, nn in zip(d.ex("hull:H0", ms, "ring_edge_mean"),
                                                       d.ex("hull:H0", ms, "nn_mean"))]
        rows.append([g, len(ms), f(rings, 0)[1:], f(d.ex("hull:H0", ms, "ring_mean"), 1)[1:],
                     f(d.ex("hull:H0", ms, "ring_largest"), 0)[1:],
                     f([r / n ** (2 / 3) for r in rings], 3)[1:], f(ratio, 1)[1:],
                     f(d.ex("hull:H0", ms, "ring_time"), 2)[1:]])
    md += ["## 1. Rings", "",
           "Layers, cities per ring, layers / N^(2/3), and the mean ring edge divided by the mean "
           "nearest-neighbour distance (both from the H0 runs; ring construction time in seconds, "
           "from the program's own clock).", "",
           table(["maps", "n maps", "rings", "cities per ring", "largest ring", "rings / N^(2/3)",
                  "ring edge / NN distance", "ring time s"], rows), ""]

    # ---- 2. vs HK, cleaned and raw, every hull merge and the best of each curve
    for p, kind in KINDS:
        for lvl, lvl_name in (("full", "TO Rn"), ("none", "raw, no cleanup")):
            rows = []
            names = [f"hull:{m}" for m in HULL] + ["curve:hull"] + \
                    [f"box:{c}:{m}:t10" for c in CURVES for m in MERGES] + [BASE, "curve:hilbert"]
            for n in names:
                if n == "curve:hull" and lvl == "full" or n == "curve:hilbert" and lvl == "full":
                    continue
                r = [label(n)]
                for s in SIZES:
                    ms = [m for m in d.group(p, s) if m in d.hk]
                    v = d.vs_hk(n, ms, lvl if n not in ("curve:hull", "curve:hilbert") else "full")
                    r.append(f(v) if v and len(v) == len(ms) else "–")
                rows.append(r)
            heads = ["setting"] + [f"{s} ({len([m for m in d.group(p, s) if m in d.hk])})" for s in SIZES]
            md += [f"## 2{'ab'[KINDS.index((p, kind))]}. {kind}, % above Held-Karp, {lvl_name}", "",
                   "Maps with a published bound in brackets. Raw for greedy TO Rn = greedy alone.", "",
                   table(heads, rows), ""]

    # ---- 3. vs greedy TO Rn, all maps
    for p, kind in KINDS:
        rows = []
        for n in [f"hull:{m}" for m in HULL] + [f"box:{c}:{m}:t10" for c in CURVES for m in MERGES]:
            r = [label(n)]
            for s in SIZES:
                ms = d.group(p, s)
                v = d.vs_base(n, ms)
                r.append(f(v) if v and len(v) == len(ms) else "–")
            rows.append(r)
        md += [f"## 3{'ab'[KINDS.index((p, kind))]}. {kind}, % longer than greedy TO Rn (every map)", "",
               table(["setting"] + [f"{s} ({len(d.group(p, s))})" for s in SIZES], rows), ""]

    # ---- 4. the chart's data
    rows = []
    for p, kind in KINDS:
        for s in SIZES:
            ms, lines = chart_lines(d, p, s)
            if not ms:
                continue
            for line, n in lines:
                clean = d.vs_hk(n, ms) if n not in ("curve:hull", "curve:hilbert") else []
                raw = d.vs_hk(n, ms, "none" if n not in ("curve:hull", "curve:hilbert") else "full")
                if (clean and len(clean) != len(ms)) or not raw or len(raw) != len(ms):
                    continue
                rows.append({"kind": kind, "size": s, "n": d.runs[(n, ms[0])]["n"], "line": line,
                             "setting": n, "maps": len(ms),
                             "clean_mean": round(mean(clean), 4) if clean else "",
                             "clean_sd": round(sd(clean), 4) if clean and sd(clean) is not None else "",
                             "raw_mean": round(mean(raw), 4),
                             "raw_sd": round(sd(raw), 4) if sd(raw) is not None else ""})
    with open(OUT / "hull_clean_vs_raw.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    md += ["## 4. Lines of the clean-vs-raw chart", "",
           "Each curve and hull at its best merge per size, chosen on the cleaned (TO Rn) mean; the raw "
           "column is the start tour of the same runs. Data: `hull_clean_vs_raw.csv`.", "",
           table(["maps", "line", "setting", "maps n", "TO Rn vs HK", "raw vs HK"],
                 [[f"{r['kind']} {r['size']}", r["line"], label(r["setting"]), r["maps"],
                   "–" if r["clean_mean"] == "" else f"{r['clean_mean']:+.2f}" + (f" ± {r['clean_sd']:.2f}" if r["clean_sd"] != "" else ""),
                   f"{r['raw_mean']:+.1f}" + (f" ± {r['raw_sd']:.1f}" if r["raw_sd"] != "" else "")] for r in rows]), ""]

    # ---- 5. shapes
    rows = []
    for fam, name in SHAPES:
        ms = d.group(fam, "10k")
        if not ms:
            continue
        hb = d.best([f"hull:{m}" for m in HULL], ms, lambda n: d.vs_base(n, ms))
        bb = d.best([f"box:{c}:{m}:t10" for c in CURVES for m in MERGES], ms, lambda n: d.vs_base(n, ms))
        hraw = d.vs_base(hb[0], ms, "none") if hb else []
        braw = d.vs_base(bb[0], ms, "none") if bb else []
        rows.append([name, len(ms), label(hb[0]) if hb else "–", f(hb[1]) if hb else "–", f(hraw, 1),
                     f(d.vs_base("curve:hull", ms), 1), label(bb[0]) if bb else "–", f(bb[1]) if bb else "–",
                     f(braw, 1), f(d.vs_base(BASE, ms, "none"), 1)])
    for p, kind in KINDS:
        ms = d.group(p, "10k")
        hb = d.best([f"hull:{m}" for m in HULL], ms, lambda n: d.vs_base(n, ms))
        bb = d.best([f"box:{c}:{m}:t10" for c in CURVES for m in MERGES], ms, lambda n: d.vs_base(n, ms))
        if hb and bb:
            rows.append([f"{kind} (E/C) for comparison", len(ms), label(hb[0]), f(hb[1]),
                         f(d.vs_base(hb[0], ms, "none"), 1), f(d.vs_base("curve:hull", ms), 1), label(bb[0]),
                         f(bb[1]), f(d.vs_base(bb[0], ms, "none"), 1), f(d.vs_base(BASE, ms, "none"), 1)])
    md += ["## 5. Map families at 10k, % longer than greedy TO Rn", "",
           "Best hull merge and best box setting (of 7 curves × M1–M5) per family, chosen on the TO Rn mean. "
           "Raw = the same runs' start tour; greedy alone = greedy TO Rn's start.", "",
           table(["maps", "n maps", "best hull", "hull TO Rn", "hull raw", "hull order", "best box",
                  "box TO Rn", "box raw", "greedy alone"], rows), ""]

    # ---- 6. time
    rows = []
    for g, p, s in groups:
        ms = d.group(p, s)
        if not d.full("hull:H0", ms):
            continue
        r = [g, len(ms), f"{mean(d.user(BASE, ms)):.2f}" if d.full(BASE, ms) else "–",
             f"{mean(d.ex('hull:H0', ms, 'ring_time')):.2f}"]
        for n in [f"hull:{m}" for m in HULL] + ["curve:hull"]:
            r.append(f"{mean(d.user(n, ms)):.2f}" if d.full(n, ms) else "–")
        rows.append(r)
    md += ["## 6. Time (user CPU seconds, whole program, mean over maps)", "",
           "Ring time is the program's own clock for the peeling (part of every hull run).", "",
           table(["maps", "n maps", "greedy TO Rn", "ring construction"] + [label(f"hull:{m}") for m in HULL]
                 + ["hull order"], rows), ""]

    # ---- 7. less cleanup
    rows = []
    for g, p, s in groups:
        ms = d.group(p, s)
        if not d.full("greedy:TO R1", ms):
            continue
        r = [g]
        for suf in (":TO R1", ":T R1"):
            gr = "greedy" + suf
            r.append(f(d.vs_base(gr, ms)))
            hb = d.best([f"hull:{m}{suf}" for m in HULL], ms, lambda n: d.vs_base(n, ms))
            bb = d.best([f"box:{c}:{m}:t10{suf}" for c in CURVES for m in MERGES], ms, lambda n: d.vs_base(n, ms))
            r.append(f"{f(hb[1])} ({label(hb[0].split(':T')[0])})" if hb else "–")
            r.append(f"{f(bb[1])} ({label(bb[0].split(':T')[0])})" if bb else "–")
        hb = d.best([f"hull:{m}" for m in HULL], ms, lambda n: d.vs_base(n, ms, "2opt"))
        bb = d.best([f"box:{c}:{m}:t10" for c in CURVES for m in MERGES], ms, lambda n: d.vs_base(n, ms, "2opt"))
        r += [f(d.vs_base(BASE, ms, "2opt")), f"{f(hb[1])} ({label(hb[0])})" if hb else "–",
              f"{f(bb[1])} ({label(bb[0])})" if bb else "–"]
        rows.append(r)
    md += ["## 7. Less cleanup, % longer than greedy TO Rn", "",
           "R1 = one round of the work queue. T Rn = the same Rn runs after their 2-opt stage. "
           "Best hull merge and best box setting for each cleanup, chosen on the mean.", "",
           table(["maps", "greedy TO R1", "best hull TO R1", "best box TO R1", "greedy T R1", "best hull T R1",
                  "best box T R1", "greedy T Rn", "best hull T Rn", "best box T Rn"], rows), ""]

    # ---- csv: every setting x group
    for g, p, s in groups + [(name, fam, "10k") for fam, name in SHAPES]:
        ms = d.group(p, s)
        for n in sorted({a for a, m in d.runs if m in ms}):
            if not d.full(n, ms):
                continue
            hk = d.vs_hk(n, ms)
            vb = d.vs_base(n, ms)
            raw = d.vs_base(n, ms, "none")
            csv_rows.append({"maps": g, "setting": n, "n_maps": len(ms),
                             "vs_hk_mean": round(mean(hk), 4) if hk else "", "vs_hk_maps": len(hk),
                             "vs_greedy_TO_Rn_mean": round(mean(vb), 4) if vb else "",
                             "vs_greedy_TO_Rn_sd": round(sd(vb), 4) if sd(vb) is not None else "",
                             "raw_vs_greedy_TO_Rn_mean": round(mean(raw), 4) if raw else "",
                             "user_s_mean": round(mean(d.user(n, ms)), 3)})
    with open(OUT / "summary.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, list(csv_rows[0]))
        w.writeheader()
        w.writerows(csv_rows)
    (OUT / "summary.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print(f"Wrote {OUT / 'summary.md'}, summary.csv, hull_clean_vs_raw.csv "
          f"({len(d.runs)} runs, {sum(d.problems.values())} problems)")
    for k, v in d.problems.items():
        print("  problem:", k, v)


if __name__ == "__main__":
    main()
