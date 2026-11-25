"""Charts for the 2D box extension (results/2d/box2/runs.csv).

    py bench/box2_charts.py

Writes six PNGs to results/2d/box2/:
  box2_overview.png             every solver and the box finalists, % above the
                                Held-Karp bound, at 10k and 100k cities
  box2_t_sweep.png              final result vs the reference for each t, per
                                curve and merge, uniform and clustered
  box2_clustered_vs_uniform.png each solver's loss on clustered data (10k)
  box2_raw_clustered_vs_uniform.png  the same with no cleanup: greedy, plain
                                curve orders and boxes only (B), every merge
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

from common import CODE_KEY, ROOT, config_label, load_bounds

SURFACE, TEXT, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e0"
UNIFORM, CLUSTERED = "#2a78d6", "#eb6834"
REF = "box:hilbert:M1:t10"
OUT = ROOT / "results" / "2d" / "box2"
CURVES = ["hilbert", "moore", "sierpinski", "knopp", "peano", "morton", "strip"]
MERGES = ["M1", "M2", "M3", "M4", "M5"]
FINALISTS = ["box:knopp:M1:t10", "box:sierpinski:M1:t10", "box:sierpinski:M2:t12",
             "box:sierpinski:M5:t4", REF]
TS = [4, 6, 8, 10, 12]


def label(name, level="full"):
    """Code name (common.config_label), marking the reference setting."""
    return config_label(name, level) + ("  (reference)" if name == REF else "")


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
    fig.text(0.2, 0.895, "Mean % above the Held-Karp bound. The plain curve orders were only run at 10k.\n"
             + CODE_KEY + ".", color=MUTED, fontsize=9.5)
    fig.subplots_adjust(left=0.2, right=0.98, top=0.83, bottom=0.14, wspace=0.08)
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
        ax.text(4.1, g, f"greedy TO Rn  {g:+.1f}", va="bottom", fontsize=8.5, color=MUTED)
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
    fig.text(0.05, 0.87, "Every line is <curve> BTO Rn <merge> at each t: % longer than the reference "
             "hilbert BTO Rn M1 t=10 on the same instance (0 = reference).\nColour = curve, line style = "
             "merge. Seed-to-seed spread is about ±0.3 (uniform) and ±1–3 (clustered).\n"
             + CODE_KEY + ".", color=MUTED, fontsize=9.5)
    fig.subplots_adjust(left=0.05, right=0.88, top=0.79, bottom=0.2, wspace=0.4)
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
    fig.suptitle("Clustered vs uniform, 10k cities: how much each solver loses", x=0.03, ha="left",
                 color=TEXT, fontsize=14, fontweight="bold", y=0.985)
    fig.text(0.03, 0.915, "Mean % above the Held-Karp bound, 3 seeds per type. Strip is off the scale "
             "and left out.\n" + CODE_KEY + ".", color=MUTED, fontsize=9.5)
    fig.subplots_adjust(left=0.27, right=0.97, top=0.89, bottom=0.1)
    fig.savefig(OUT / "box2_clustered_vs_uniform.png", dpi=170, facecolor=SURFACE)


def raw_clustered_vs_uniform(bounds, L):
    """No cleanup at all: greedy, each curve's plain order, and that curve's
    boxes (B) with every merge at t = 10, uniform vs clustered, 10k cities."""
    curves = [c for c in CURVES if c != "strip"]           # strip is off the scale
    rows = [("greedy", "full")]
    for c in curves:
        rows += [None, (f"curve:{c}", "full")] + [(f"box:{c}:{m}:t10", "none") for m in MERGES]
    u_i, c_i = insts(L, "uniform", "10k"), insts(L, "clustered", "10k")
    fig, ax = plt.subplots(figsize=(12.5, 13), facecolor=SURFACE)
    clean(ax)
    ypos, labels, y = [], [], 0
    for r in rows:
        if r is None:
            y += 0.6
            continue
        name, lvl = r
        u, c = hk_gap(L, bounds, name, u_i, lvl), hk_gap(L, bounds, name, c_i, lvl)
        ax.plot([u, c], [y, y], color=GRID, lw=5, solid_capstyle="round", zorder=1)
        ax.scatter([u], [y], s=70, color=UNIFORM, zorder=3, edgecolor=SURFACE, linewidth=1.5)
        ax.scatter([c], [y], s=70, color=CLUSTERED, zorder=3, edgecolor=SURFACE, linewidth=1.5)
        ax.annotate(f"{c - u:+.1f} pts", xy=(max(u, c), y), xytext=(9, 0), textcoords="offset points",
                    va="center", fontsize=9.5, fontweight="bold", color=TEXT)
        ypos.append(y); labels.append(config_label(name, lvl)); y += 1
    for kind, col in (("uniform", UNIFORM), ("clustered", CLUSTERED)):
        g = hk_gap(L, bounds, "greedy", insts(L, kind, "10k"))
        ax.axvline(g, color=col, lw=1, ls=(0, (4, 3)), zorder=0)
        left = kind == "uniform"                          # labels either side of the two lines
        ax.text(g - 0.8 if left else g + 0.8, y - 0.2, f"greedy, {kind}\n{g:.1f}%", color=col,
                fontsize=8.5, va="top", ha="right" if left else "left")
    ax.set_yticks(ypos, labels)
    for lbl in ax.get_yticklabels():
        lbl.set_color(TEXT)
    ax.set_ylim(y + 1.6, -0.8)
    ax.set_xlim(0, 125)
    ax.set_xlabel("% above Held-Karp bound (lower is better)")
    fig.legend(handles=[plt.Line2D([], [], marker="o", ls="", ms=8, color=UNIFORM, label="uniform"),
                        plt.Line2D([], [], marker="o", ls="", ms=8, color=CLUSTERED, label="clustered"),
                        plt.Line2D([], [], ls="", label="label = clustered minus uniform")],
               loc="lower left", bbox_to_anchor=(0.2, 0.0), ncol=3, frameon=False, labelcolor=TEXT)
    fig.suptitle("No cleanup: raw greedy, raw curve orders and raw boxes, clustered vs uniform (10k)",
                 x=0.03, ha="left", color=TEXT, fontsize=14, fontweight="bold", y=0.99)
    fig.text(0.03, 0.965, "No 2-opt, no Or-opt. Each curve: its plain order, then its boxes (B) with each "
             "merge, all t = 10. Mean % above the Held-Karp bound, 3 seeds per type.\nStrip is off the "
             "scale (1,000–7,000%) and left out.\n" + CODE_KEY + ".", color=MUTED, fontsize=9.5, va="top")
    fig.subplots_adjust(left=0.2, right=0.97, top=0.915, bottom=0.06)
    fig.savefig(OUT / "box2_raw_clustered_vs_uniform.png", dpi=170, facecolor=SURFACE)


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
            stage = "after cleanup (BTO Rn)" if lvl == "full" else "boxes only (B)"
            ax.set_title(f"{kind}, {stage}",
                         loc="left", color=TEXT, fontsize=11.5, fontweight="bold")
    fig.suptitle("Curve × merge at t = 10, 10k cities: % vs the reference at the same stage", x=0.06,
                 ha="left", color=TEXT, fontsize=14, fontweight="bold", y=0.985)
    fig.text(0.06, 0.955, "Each cell: <curve> <merge> t=10, % longer than hilbert M1 t=10 at the same stage. "
             "Blue = shorter than the reference, red = longer; mean of 3 seeds.\nColours are capped at ±12 "
             "(top) and ±30 (bottom); numbers are exact.\n" + CODE_KEY + ".", color=MUTED, fontsize=9.5,
             va="top")
    fig.subplots_adjust(left=0.08, right=0.98, top=0.87, bottom=0.04, hspace=0.25, wspace=0.18)
    fig.savefig(OUT / "box2_curve_merge.png", dpi=170, facecolor=SURFACE)


def stage_a(bounds, L):
    """Every Stage A result: each box setting at B, BT and BTO, every plain
    curve order and every baseline, uniform and clustered, 10k cities."""
    stage_colour = {"none": "#b8d4f5", "2opt": "#5b9be0", "full": "#0d2f5c"}
    rows = []                                          # (label, kind, name)
    for c in CURVES:
        for m in MERGES:
            rows.append((f"{c}  {m}", "box", f"box:{c}:{m}:t10"))
        rows.append(None)                              # gap between curves
    rows += [(config_label(f"curve:{c}"), "plain", f"curve:{c}") for c in CURVES]
    rows.append(None)
    for b in ["spacefill", "nn", "greedy", "twoopt", "twoopt_or", "greedy+2opt+oropt"]:
        rows.append((config_label(b), "plain", b))
    ypos, labels, y = [], [], 0
    for r in rows:                                     # rows top to bottom
        if r is None:
            y += 0.6
            continue
        ypos.append(y); labels.append(r[0]); y += 1

    fig, axes = plt.subplots(1, 2, figsize=(15, 15.5), facecolor=SURFACE, sharey=True)
    for ax, kind in zip(axes, ["uniform", "clustered"]):
        clean(ax)
        ax.set_xscale("log")
        ii = insts(L, kind, "10k")
        k = 0
        for r in rows:
            if r is None:
                continue
            label, typ, name = r
            yk = ypos[k]; k += 1
            if typ == "box":
                pts = [(lvl, hk_gap(L, bounds, name, ii, lvl)) for lvl in ("none", "2opt", "full")]
                xs = [v for _, v in pts if v is not None]
                ax.plot(xs, [yk] * len(xs), color=GRID, lw=3, zorder=1, solid_capstyle="round")
                for lvl, v in pts:
                    ax.scatter([v], [yk], s=34, color=stage_colour[lvl], zorder=3,
                               edgecolor=SURFACE, linewidth=0.8)
                fin = pts[-1][1]
                ax.text(fin / 1.12, yk, f"{fin:.1f}", ha="right", va="center", fontsize=7.5, color=TEXT)
            else:
                v = hk_gap(L, bounds, name, ii, "full")
                col = "#e34948" if name == "greedy+2opt+oropt" else MUTED
                ax.scatter([v], [yk], s=40, marker="D", color=col, zorder=3, edgecolor=SURFACE, linewidth=0.8)
                ax.text(v * 1.12, yk, f"{v:.1f}", ha="left", va="center", fontsize=7.5, color=TEXT)
        g = hk_gap(L, bounds, "greedy+2opt+oropt", ii, "full")
        ax.axvline(g, color="#e34948", lw=1, ls=(0, (4, 3)), zorder=0)
        ax.text(g * 1.05, -1.2, f"greedy TO Rn {g:.1f}%", color="#e34948", fontsize=8.5, va="center")
        ax.set_xlim(2.5, 15000)
        ax.set_xticks([3, 5, 10, 20, 50, 100, 1000, 10000], ["3", "5", "10", "20", "50", "100", "1,000", "10,000"])
        ax.minorticks_off()
        ax.set_xlabel("% above Held-Karp bound, log scale (lower is better)")
        ax.set_title(f"{kind}, 10k cities (3 seeds)", loc="left", color=TEXT, fontsize=12, fontweight="bold")
    axes[0].set_yticks(ypos, labels)
    axes[0].invert_yaxis()
    axes[0].set_ylim(y, -2)
    for lbl in axes[0].get_yticklabels():
        lbl.set_color(TEXT); lbl.set_fontsize(9)
    handles = [plt.Line2D([], [], marker="o", ls="", ms=7, color=stage_colour[l], label=t)
               for l, t in (("none", "B: boxes only"), ("2opt", "BT Rn: + 2-opt"), ("full", "BTO Rn: + Or-opt (number shown)"))]
    handles += [plt.Line2D([], [], marker="D", ls="", ms=6, color=MUTED, label="plain order or baseline (one stage)")]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.14, -0.004), ncol=4, frameon=False,
               labelcolor=TEXT, fontsize=9.5)
    fig.suptitle("Stage A: every result at 10k cities", x=0.14, ha="left", color=TEXT,
                 fontsize=15, fontweight="bold", y=0.993)
    fig.text(0.14, 0.975, "Box rows (\"<curve>  <merge>\", all t = 10) show the same run at three stages: "
             "B, then BT Rn, then BTO Rn. Mean of 3 seeds.\n" + CODE_KEY + ".", color=MUTED, fontsize=9.5, va="top")
    fig.subplots_adjust(left=0.14, right=0.98, top=0.925, bottom=0.07, wspace=0.06)
    fig.savefig(OUT / "box2_stageA_all.png", dpi=150, facecolor=SURFACE)


def main():
    style()
    bounds, L = load()
    stage_a(bounds, L)
    curve_merge(L)
    t_sweep(L)
    clustered_vs_uniform(bounds, L)
    raw_clustered_vs_uniform(bounds, L)
    overview(bounds, L)
    print("Wrote", ", ".join(p.name for p in sorted(OUT.glob("box2_*.png"))))


if __name__ == "__main__":
    main()
