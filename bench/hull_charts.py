"""Charts for the hull-ring study (results/2d/hull/runs.csv).

    py bench/hull_charts.py      # after bench/hull_report.py; writes results/2d/hull_*.png

  hull_clean_vs_raw.png   % above the Held-Karp bound against the number of cities.
                          Left half: after 2-opt + Or-opt (TO Rn); right half: the same
                          runs' start tours with no cleanup.  Built from
                          results/2d/hull/hull_clean_vs_raw.csv (hull_report.py).
  hull_shapes.png         map families at 10k: best hull and best box setting against
                          greedy TO Rn, cleaned and raw
"""

import csv
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from common import CODE_KEY, ROOT
from hull_report import BASE, CURVES, HULL, MERGES, OUT, SHAPES, Data, mean, sd

SURFACE, TEXT, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e0"
HULL_COL = "#eb6834"                      # accent 1: hull rings
GREEDY_COL = TEXT                         # accent 2: greedy, the baseline
CURVE_COL = ["#8f8d88", "#b3b1ab"]        # the 7 curves: quiet greys, alternating
SIZES = {"1k": 1000, "10k": 10000, "100k": 100000, "1M": 1000000}
CHARTS = ROOT / "results" / "2d"


def style():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.edgecolor": MUTED,
                         "axes.labelcolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED,
                         "figure.facecolor": SURFACE})


def clean(ax, grid_axis="y"):
    ax.set_facecolor(SURFACE)
    ax.grid(True, axis=grid_axis, color=GRID, linewidth=1, which="major")
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def spread(items, ylim, log, frac=0.042):
    """label positions at least frac of the axis height apart, kept inside the axis;
    sorts items by y"""
    key = (lambda v: math.log10(v)) if log else (lambda v: v)
    back = (lambda v: 10 ** v) if log else (lambda v: v)
    lo, hi = key(ylim[0]), key(ylim[1])
    gap = frac * (hi - lo)
    items.sort(key=lambda e: e[0])
    pos = [key(e[0]) for e in items]
    for i in range(1, len(pos)):
        pos[i] = max(pos[i], pos[i - 1] + gap)
    over = pos[-1] - (hi - gap / 2)
    if over > 0:                                   # stack runs off the top: push it down
        pos = [p - over for p in pos]
        for i in range(len(pos) - 2, -1, -1):      # and keep the gaps from the top down
            pos[i] = min(pos[i], pos[i + 1] - gap)
    return [back(p) for p in pos]


def line_style(line):
    if line == "hull":
        return dict(color=HULL_COL, lw=2.6, zorder=5)
    if line == "greedy":
        return dict(color=GREEDY_COL, lw=2.6, zorder=6)
    if line == "hull order":
        return dict(color=HULL_COL, lw=1.8, ls=":", zorder=4)
    if line == "Hilbert order":
        return dict(color=MUTED, lw=1.8, ls=":", zorder=4)
    return dict(color=CURVE_COL[CURVES.index(line) % 2], lw=1.3, zorder=3)


ORDER = CURVES + ["Hilbert order", "hull order", "hull", "greedy"]     # drawing order, every panel


def clean_vs_raw():
    rows = list(csv.DictReader(open(OUT / "hull_clean_vs_raw.csv", newline="")))
    fig = plt.figure(figsize=(19, 7.4))
    gs = fig.add_gridspec(1, 5, width_ratios=[1, 1, 0.16, 1, 1], left=0.045, right=0.985, top=0.80,
                          bottom=0.15, wspace=0.42)
    axes = [fig.add_subplot(gs[0, i]) for i in (0, 1, 3, 4)]
    axes[3].sharey(axes[2])          # raw: one log scale; cleaned: each panel its own range
    panels = [(axes[0], "uniform", "clean"), (axes[1], "clustered", "clean"),
              (axes[2], "uniform", "raw"), (axes[3], "clustered", "raw")]
    for ax, kind, lvl in panels:
        clean(ax)
        raw = lvl == "raw"
        if raw:
            ax.set_yscale("log")
        ax.set_xscale("log")
        ends = []
        for line in ORDER:
            pts = sorted([r for r in rows if r["kind"] == kind and r["line"] == line],
                         key=lambda r: int(r["n"]))
            if raw:
                pts = [r for r in pts if r["raw_mean"]]
                ys = [float(r["raw_mean"]) for r in pts]
                es = [float(r["raw_sd"]) if r["raw_sd"] else 0 for r in pts]
            else:
                pts = [r for r in pts if r["clean_mean"]]
                ys = [float(r["clean_mean"]) for r in pts]
                es = [float(r["clean_sd"]) if r["clean_sd"] else 0 for r in pts]
            if not pts:
                continue
            xs = [int(r["n"]) for r in pts]
            st = line_style(line)
            ax.errorbar(xs, ys, yerr=es, marker="o", ms=4.5 if st["lw"] > 2 else 3.2, capsize=2.5,
                        elinewidth=1, **st)
            name = {"greedy": "greedy TO Rn" if not raw else "greedy alone",
                    "hull": "hull BTO Rn" if not raw else "hull B"}.get(line, line)
            ends.append([ys[-1], f"{name} {ys[-1]:+.1f}" if not raw else f"{name} {ys[-1]:.0f}%",
                         st["color"] if line not in CURVES else MUTED,
                         "bold" if line in ("hull", "greedy") else "normal", xs[-1]])
        cap = OUT / "reversal_cap_test.csv"            # one untimed run without the reversal cap
        if kind == "uniform" and not raw and cap.exists():
            for r in csv.DictReader(open(cap, newline="")):
                if r["setting"] == "hull:M1" and int(r["max_reversal"]) > 1000000:
                    v = float(r["vs_hk_pct"])
                    ax.plot([1000000], [v], marker="o", ms=8, mfc=SURFACE, mec=HULL_COL, mew=2, zorder=7)
                    ax.annotate(f"hull M1 without the\nreversal cap: {v:+.1f}", (1000000, v), xytext=(-12, -30),
                                textcoords="offset points", ha="right", fontsize=8.3, color=HULL_COL)
        top_x = max(e[4] for e in ends)
        for (y, txt, col, w, x), q in zip(ends, spread(ends, ax.get_ylim(), raw)):
            ax.text(top_x * 1.45, q, txt, color=col, fontsize=8.6, fontweight=w, va="center")
        ax.set_xticks(list(SIZES.values()), list(SIZES))
        ax.minorticks_off()
        ax.set_xlim(650, top_x * 14)
        ax.set_xlabel("cities (log scale)")
        ax.set_title(f"{kind}" + ("" if kind == "uniform" else " (no 1M bound)"), loc="left", color=TEXT,
                     fontsize=11.5, fontweight="bold")
    axes[0].set_ylabel("% above the Held-Karp bound (linear)")
    axes[2].set_ylabel("% above the Held-Karp bound (log scale)")
    axes[2].yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}%"))
    for x, txt in ((0.045, "With cleanup: start + 2-opt + Or-opt, rounds until clean (TO Rn)"),
                   (0.555, "Without cleanup: the start tours of the same runs")):
        fig.text(x, 0.865, txt, color=TEXT, fontsize=12.5, fontweight="bold")
    fig.suptitle("2D hull rings against the 7 curves and greedy: lower is better", x=0.045, ha="left",
                 color=TEXT, fontsize=15, fontweight="bold", y=0.975)
    fig.text(0.045, 0.925, "Each curve (grey) and hull (orange) at its best merge per size, chosen on the "
             "cleaned mean (t = 10, K = 16); the raw panels show the same runs before cleanup. Mean ± 1 sd over "
             "the maps with a published bound: 3 at 1k and 10k, 2 at 100k, 1 at 1M (uniform only).",
             color=MUTED, fontsize=9.5)
    fig.text(0.045, 0.045, CODE_KEY + ". For hull, B = rings + merge.", color=MUTED, fontsize=9)
    fig.text(0.045, 0.018, "Dotted: an order alone, no merge optimisation and no cleanup. The two halves use "
             "different y scales; a curve's best merge may change from one size to the next. Every run uses the 2-opt "
             "reversal cap (100,000 cities) except the hollow point.", color=MUTED, fontsize=9)
    fig.savefig(CHARTS / "hull_clean_vs_raw.png", dpi=150, facecolor=SURFACE)
    plt.close(fig)


def shapes(d):
    fams = [(p, name) for p, name in SHAPES] + [("E", "uniform"), ("C", "clustered")]
    labels, hull_c, hull_r, box_c, box_r, notes = [], [], [], [], [], []
    for p, name in fams:
        ms = d.group(p, "10k")
        hb = d.best([f"hull:{m}" for m in HULL], ms, lambda n: d.vs_base(n, ms))
        bb = d.best([f"box:{c}:{m}:t10" for c in CURVES for m in MERGES], ms, lambda n: d.vs_base(n, ms))
        if not hb or not bb:
            continue
        labels.append(name)
        hull_c.append((mean(hb[1]), sd(hb[1]) or 0))
        box_c.append((mean(bb[1]), sd(bb[1]) or 0))
        hull_r.append(mean(d.vs_base(hb[0], ms, "none")))
        box_r.append(mean(d.vs_base(bb[0], ms, "none")))
        notes.append((hb[0].split(":")[1], bb[0].split(":")[1] + " " + bb[0].split(":")[2]))
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(15, 5.6), sharey=True)
    ys = list(range(len(labels)))
    for ax, hv, bv, raw in ((a1, hull_c, box_c, False), (a2, hull_r, box_r, True)):
        clean(ax, "x")
        h = [v[0] if not raw else v for v in hv]
        b = [v[0] if not raw else v for v in bv]
        ax.barh([y - 0.2 for y in ys], b, height=0.36, color="#a9a7a1", label="best box setting")
        ax.barh([y + 0.2 for y in ys], h, height=0.36, color=HULL_COL, label="best hull merge")
        if not raw:
            ax.errorbar(b, [y - 0.2 for y in ys], xerr=[v[1] for v in bv], fmt="none", ecolor=MUTED, capsize=2)
            ax.errorbar(h, [y + 0.2 for y in ys], xerr=[v[1] for v in hv], fmt="none", ecolor=MUTED, capsize=2)
        eb = [0 if raw else v[1] for v in bv]
        eh = [0 if raw else v[1] for v in hv]
        for y, v, e in zip(ys, b, eb):                 # value after the error bar
            ax.text(v + e, y - 0.2, f"  {v:+.1f}", va="center", fontsize=8.5, color=TEXT)
        for y, v, e in zip(ys, h, eh):
            ax.text(v + e, y + 0.2, f"  {v:+.1f}", va="center", fontsize=8.5, color=TEXT, fontweight="bold")
        ax.axvline(0, color=TEXT, lw=1.4)
        ax.set_xlabel("% longer than greedy TO Rn on the same map" + (" (raw start tours)" if raw else ""))
        ax.set_title("after 2-opt + Or-opt (TO Rn)" if not raw else "before cleanup (raw)", loc="left",
                     color=TEXT, fontsize=11.5, fontweight="bold")
        if raw:
            ax.set_xscale("symlog", linthresh=10)
            ax.set_xlim(0, max(h + b) * 4)
            ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:+g}"))
        else:
            ax.set_xlim(min(0, min(h + b)) - 0.5, max(x + e for x, e in zip(h + b, eh + eb)) * 1.2)
    a1.set_yticks(ys, [f"{l}\n(hull {n[0]}, box {n[1]})" for l, n in zip(labels, notes)])
    a1.invert_yaxis()
    a1.legend(loc="upper right", frameon=False, labelcolor=TEXT)
    fig.suptitle("Where hull rings do well: map families at 10k cities", x=0.04, ha="left", color=TEXT,
                 fontsize=14, fontweight="bold", y=0.975)
    fig.text(0.04, 0.905, "Best hull merge and best box setting (7 curves × M1–M5, t = 10) per family, mean "
             "± 1 sd over 3 maps; 0 = greedy TO Rn. Raw: the same runs' start tours (symmetric log scale).",
             color=MUTED, fontsize=9.5)
    fig.text(0.04, 0.02, CODE_KEY, color=MUTED, fontsize=8.5)
    fig.subplots_adjust(left=0.17, right=0.98, top=0.82, bottom=0.16, wspace=0.08)
    fig.savefig(CHARTS / "hull_shapes.png", dpi=150, facecolor=SURFACE)
    plt.close(fig)


def main():
    style()
    clean_vs_raw()
    shapes(Data())
    print("Wrote results/2d/hull_clean_vs_raw.png and hull_shapes.png")


if __name__ == "__main__":
    main()
