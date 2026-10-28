"""Charts for the 2D box extension (results/box2/runs.csv).

    py bench/box2_charts.py

Writes four PNGs to results/box2/:
  box2_overview.png             every solver and the box finalists, % above the
                                Held-Karp bound, at 10k and 100k cities
  box2_t_sweep.png              final result vs the reference for each t, per
                                curve and merge, uniform and clustered
  box2_clustered_vs_uniform.png each solver's loss on clustered data (10k)
  box2_curve_merge.png          curve x merge grid at t = 10, after full
                                cleanup and before any cleanup
Gaps use the Held-Karp bound (refs/bounds): all 10k seeds have one, 2 of
the 3 100k seeds.  "vs reference" = % longer than box:hilbert:M1:t10 on
the same instance.
"""

import csv
import statistics
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

from common import ROOT, load_bounds

SURFACE, TEXT, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e0"
UNIFORM, CLUSTERED = "#2a78d6", "#eb6834"
REF = "box:hilbert:M1:t10"
OUT = ROOT / "results" / "box2"
CURVES = ["hilbert", "moore", "sierpinski", "knopp", "peano", "morton", "strip"]
MERGES = ["M1", "M2", "M3", "M4", "M5"]
FINALISTS = ["box:knopp:M1:t10", "box:sierpinski:M1:t10", "box:sierpinski:M2:t12",
             "box:sierpinski:M5:t4", REF]
TS = [4, 6, 8, 10, 12]


def label(name):
    if name.startswith("box:"):
        _, c, m, t = name.split(":")
        return f"box {c} {m} t={t[1:]}" + ("  (reference)" if name == REF else "")
    if name.startswith("curve:"):
        return f"{name[6:]} order (no boxes)"
    return name


def style():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.edgecolor": MUTED, "axes.labelcolor": MUTED,
                         "xtick.color": MUTED, "ytick.color": MUTED})


def clean(ax, grid_axis="x"):
    ax.set_facecolor(SURFACE)
    ax.grid(True, axis=grid_axis, color=GRID, linewidth=1)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def load():
    bounds = load_bounds(ROOT / "refs" / "bounds")
    L = {}                                   # (name, instance) -> {level: length}
    for r in csv.DictReader(open(OUT / "runs.csv", newline="")):
        if r["status"] != "OK":
            continue
        ex = dict(kv.split("=", 1) for kv in r["extras"].split(";") if "=" in kv)
        d = {"full": float(r["tour_length"])}
        if "start_length" in ex: d["none"] = float(ex["start_length"])
        if "after_2opt_length" in ex: d["2opt"] = float(ex["after_2opt_length"])
        L[(r["algorithm"], r["instance"])] = d
    return bounds, L


def insts(L, kind, size):
    p = "E" if kind == "uniform" else "C"
    return sorted({i for _, i in L if i.startswith(p + size + ".")})


def mean(v):
    return statistics.mean(v) if v else None


def hk_gap(L, bounds, name, instances, lvl="full"):
    return mean([(L[(name, i)][lvl] / bounds[i][0] - 1) * 100
                 for i in instances if (name, i) in L and i in bounds and lvl in L[(name, i)]])


def vs_ref(L, name, instances, lvl="full"):
    return mean([(L[(name, i)][lvl] / L[(REF, i)][lvl] - 1) * 100
                 for i in instances if (name, i) in L and (REF, i) in L and lvl in L[(name, i)]])


def overview(bounds, L):
    rows = ["spacefill", "curve:sierpinski", "nn", "greedy", REF, "box:sierpinski:M5:t4",
            "box:sierpinski:M2:t12", "box:sierpinski:M1:t10", "box:knopp:M1:t10",
            "twoopt", "twoopt_or", "greedy+2opt+oropt"]
    fig, axes = plt.subplots(1, 2, figsize=(15, 6.6), facecolor=SURFACE, sharey=True)
    for ax, size in zip(axes, ["10k", "100k"]):
        clean(ax)
        u_i, c_i = insts(L, "uniform", size), insts(L, "clustered", size)
        nb = len([i for i in u_i if i in bounds])
        shown = [r for r in rows if hk_gap(L, bounds, r, u_i) is not None]
        for k, r in enumerate(rows):
            if r not in shown:
                continue
            for off, data, col in ((-0.2, u_i, UNIFORM), (0.2, c_i, CLUSTERED)):
                g = hk_gap(L, bounds, r, data)
                if g is None:
                    continue
                ax.barh(k + off, g, height=0.38, color=col)
                ax.text(g + 0.6, k + off, f"{g:.1f}%", va="center", fontsize=8.5, color=TEXT)
        ax.set_yticks(range(len(rows)), [label(r) for r in rows])
        ax.invert_yaxis()
        ax.set_xlabel("% above Held-Karp bound (lower is better)")
        ax.set_title(f"{size} cities  ·  {nb} seeds per type with a bound", loc="left",
                     color=TEXT, fontsize=11.5, fontweight="bold")
        ax.set_xlim(0, 64)
    for lbl in axes[0].get_yticklabels():
        lbl.set_color(TEXT)
    fig.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=UNIFORM, label="uniform"),
                        plt.Rectangle((0, 0), 1, 1, color=CLUSTERED, label="clustered")],
               loc="lower left", bbox_to_anchor=(0.2, 0.0), ncol=2, frameon=False, labelcolor=TEXT)
    fig.suptitle("Every solver and the box finalists: tour quality", x=0.2, ha="left",
                 color=TEXT, fontsize=14, fontweight="bold", y=0.985)
    fig.text(0.2, 0.915, "Mean % above the Held-Karp bound. Box settings are run with full cleanup "
             "(2-opt, then 2-opt + Or-opt). The plain curve orders were only run at 10k.", color=MUTED, fontsize=9.5)
    fig.subplots_adjust(left=0.2, right=0.98, top=0.85, bottom=0.14, wspace=0.08)
    fig.savefig(OUT / "box2_overview.png", dpi=170, facecolor=SURFACE)


def t_sweep(L):
    combos = [("hilbert", "M1"), ("sierpinski", "M1"), ("knopp", "M1"), ("sierpinski", "M2"),
              ("moore", "M2"), ("knopp", "M2"), ("sierpinski", "M5"), ("moore", "M5"), ("knopp", "M5")]
    colours = {"hilbert": "#52514e", "sierpinski": "#2a78d6", "moore": "#1baf7a", "knopp": "#eb6834"}
    styles = {"M1": "-", "M2": "--", "M5": ":"}
    markers = {"M1": "o", "M2": "s", "M5": "^"}
    fig, axes = plt.subplots(1, 2, figsize=(15, 6.2), facecolor=SURFACE)
    for ax, kind in zip(axes, ["uniform", "clustered"]):
        clean(ax, "y")
        ii = insts(L, kind, "10k")
        ends = []
        for c, m in combos:
            ys = [vs_ref(L, f"box:{c}:{m}:t{t}", ii) for t in TS]
            if sum(y is not None for y in ys) < 2:
                continue
            xs = [t for t, y in zip(TS, ys) if y is not None]
            ys = [y for y in ys if y is not None]
            ax.plot(xs, ys, color=colours[c], ls=styles[m], marker=markers[m], lw=2, ms=6)
            ends.append([ys[-1], f"{c} {m}  {ys[-1]:+.1f}", colours[c]])
        ax.axhline(0, color=MUTED, lw=1)
        g = vs_ref(L, "greedy+2opt+oropt", ii)
        ax.axhline(g, color=MUTED, lw=1.2, ls=(0, (4, 3)))
        ax.text(4.1, g, f"greedy start + same cleanup  {g:+.1f}", va="bottom", fontsize=8.5, color=MUTED)
        ends.sort(key=lambda e: e[0])
        lo, hi = ax.get_ylim()
        pos = [e[0] for e in ends]
        for i in range(1, len(pos)):
            pos[i] = max(pos[i], pos[i - 1] + (hi - lo) * 0.045)
        for (y, text, col), p in zip(ends, pos):
            ax.annotate(text, xy=(12, y), xytext=(12.35, p), color=col, fontsize=8.8,
                        va="center", fontweight="bold", annotation_clip=False)
        ax.set_xticks(TS, [f"t = {t}" for t in TS])
        ax.set_xlim(3.6, 12.4)
        ax.set_xlabel("box size t")
        ax.set_ylabel("% vs reference, after full cleanup (lower is better)")
        ax.set_title(f"{kind}, 10k cities (3 seeds)", loc="left", color=TEXT, fontsize=11.5, fontweight="bold")
    handles = [plt.Line2D([], [], color=MUTED, ls=styles[m], marker=markers[m], lw=2, label=m) for m in styles]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.05, 0.0), ncol=3, frameon=False,
               labelcolor=TEXT, title="merge (line style)", title_fontsize=9)
    fig.suptitle("Box size t after cleanup: small, inconsistent effects", x=0.05, ha="left",
                 color=TEXT, fontsize=14, fontweight="bold", y=0.985)
    fig.text(0.05, 0.9, "% longer than the reference (Hilbert, M1, t = 10) on the same instance; 0 = reference. "
             "Colour = curve, line style = merge. Seed-to-seed spread is about ±0.3 (uniform) and ±1–3 (clustered).",
             color=MUTED, fontsize=9.5)
    fig.subplots_adjust(left=0.05, right=0.87, top=0.83, bottom=0.2, wspace=0.42)
    fig.savefig(OUT / "box2_t_sweep.png", dpi=170, facecolor=SURFACE)


def clustered_vs_uniform(bounds, L):
    rows = ["greedy+2opt+oropt", "twoopt_or", "twoopt", "box:knopp:M1:t10", "box:sierpinski:M1:t10",
            "box:sierpinski:M2:t12", "box:sierpinski:M5:t4", REF, "greedy", "nn",
            "curve:sierpinski", "curve:knopp", "curve:moore", "curve:hilbert", "curve:peano",
            "curve:morton"]
    u_i, c_i = insts(L, "uniform", "10k"), insts(L, "clustered", "10k")
    fig, ax = plt.subplots(figsize=(12.5, 8.2), facecolor=SURFACE)
    clean(ax)
    for k, r in enumerate(rows):
        u, c = hk_gap(L, bounds, r, u_i), hk_gap(L, bounds, r, c_i)
        if u is None or c is None:
            continue
        ax.plot([u, c], [k, k], color=GRID, lw=5, solid_capstyle="round", zorder=1)
        ax.scatter([u], [k], s=70, color=UNIFORM, zorder=3, edgecolor=SURFACE, linewidth=1.5)
        ax.scatter([c], [k], s=70, color=CLUSTERED, zorder=3, edgecolor=SURFACE, linewidth=1.5)
        ax.annotate(f"{c - u:+.1f} pts", xy=(max(u, c), k), xytext=(9, 0), textcoords="offset points",
                    va="center", fontsize=9.5, fontweight="bold", color=TEXT)
    ax.set_yticks(range(len(rows)), [label(r) for r in rows])
    for lbl in ax.get_yticklabels():
        lbl.set_color(TEXT)
    ax.invert_yaxis()
    ax.set_xlim(0, 125)
    ax.set_xlabel("% above Held-Karp bound (lower is better)")
    fig.legend(handles=[plt.Line2D([], [], marker="o", ls="", ms=8, color=UNIFORM, label="uniform"),
                        plt.Line2D([], [], marker="o", ls="", ms=8, color=CLUSTERED, label="clustered"),
                        plt.Line2D([], [], ls="", label="label = clustered minus uniform")],
               loc="lower left", bbox_to_anchor=(0.27, 0.0), ncol=3, frameon=False, labelcolor=TEXT)
    fig.suptitle("Clustered vs uniform, 10k cities: how much each solver loses", x=0.27, ha="left",
                 color=TEXT, fontsize=14, fontweight="bold", y=0.985)
    fig.text(0.27, 0.925, "Mean % above the Held-Karp bound, 3 seeds per type. Box settings after full "
             "cleanup;\ncurve orders have no boxes and no cleanup. Strip is off the scale and left out.",
             color=MUTED, fontsize=9.5)
    fig.subplots_adjust(left=0.27, right=0.97, top=0.9, bottom=0.1)
    fig.savefig(OUT / "box2_clustered_vs_uniform.png", dpi=170, facecolor=SURFACE)


def curve_merge(L):
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), facecolor=SURFACE)
    for row, lvl in enumerate(["full", "none"]):
        for col, kind in enumerate(["uniform", "clustered"]):
            ax = axes[row][col]
            ii = insts(L, kind, "10k")
            grid = [[vs_ref(L, f"box:{c}:{m}:t10", ii, lvl) for m in MERGES] for c in CURVES]
            lim = 12 if lvl == "full" else 30
            ax.imshow([[max(-lim, min(lim, v)) for v in r] for r in grid], cmap="RdBu_r",
                      norm=TwoSlopeNorm(vcenter=0, vmin=-lim, vmax=lim), aspect="auto")
            for i, r in enumerate(grid):
                for j, v in enumerate(r):
                    txt = f"{v:+.1f}" if abs(v) < 100 else f"{v:+.0f}"
                    ax.text(j, i, txt, ha="center", va="center", fontsize=9.5,
                            color="white" if abs(v) > lim * 0.6 else TEXT)
            ax.set_xticks(range(len(MERGES)), MERGES)
            ax.set_yticks(range(len(CURVES)), CURVES)
            ax.tick_params(length=0)
            for s in ax.spines.values():
                s.set_visible(False)
            title = "after full cleanup" if lvl == "full" else "box tour, before cleanup"
            ax.set_title(f"{kind}: {title}", loc="left", color=TEXT, fontsize=11.5, fontweight="bold")
    fig.suptitle("Curve × merge at t = 10, 10k cities: % vs the reference (Hilbert, M1)", x=0.06,
                 ha="left", color=TEXT, fontsize=14, fontweight="bold", y=0.985)
    fig.text(0.06, 0.945, "Blue = shorter than the reference, red = longer; mean of 3 seeds. Colours are "
             "capped at ±12 (after cleanup) and ±30 (before); numbers are exact.", color=MUTED, fontsize=9.5)
    fig.subplots_adjust(left=0.08, right=0.98, top=0.9, bottom=0.04, hspace=0.25, wspace=0.18)
    fig.savefig(OUT / "box2_curve_merge.png", dpi=170, facecolor=SURFACE)


def main():
    style()
    bounds, L = load()
    curve_merge(L)
    t_sweep(L)
    clustered_vs_uniform(bounds, L)
    overview(bounds, L)
    print("Wrote", ", ".join(p.name for p in sorted(OUT.glob("box2_*.png"))))


if __name__ == "__main__":
    main()
