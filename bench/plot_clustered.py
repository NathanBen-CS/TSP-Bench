"""Chart: how much worse each solver does on clustered data than on uniform data.

    py bench/plot_clustered.py               # writes results/2d/clustered_vs_uniform.png

Both data types are measured against the same reference, the Held-Karp lower
bound from refs/bounds, so the two dots of a row are directly comparable.
Bounds exist for the same seeds of both types at 1k (10), 10k (3) and
100k (2: .0 and .1); seeds without a bound (E100k.2, C100k.2, all 1M
clustered) are left out, so 1M is not shown.
"""

import csv
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from common import CODE_KEY, ROOT, config_label, load_bounds

SOLVERS = [                                   # best first, as in solver_performance.png;
    "twoopt_or", "twoopt", "boxes+2opt+oropt",  # labels from config_label (common.py)
    "greedy", "nn", "boxes-only", "spacefill",
]
SOLVERS = [(k, config_label(k)) for k in SOLVERS]
SIZES = [("1k", 1000), ("10k", 10000), ("100k", 100000)]
UNIFORM, CLUSTERED = "#2a78d6", "#eb6834"
SURFACE, TEXT, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e0"


def main():
    bounds = load_bounds(ROOT / "refs" / "bounds")
    gaps = defaultdict(list)                  # (solver, type, n) -> [% above HK]
    with open(ROOT / "results" / "2d" / "box" / "runs.csv", newline="") as f:
        for r in csv.DictReader(f):
            inst = r["instance"]
            if r["status"] != "OK" or inst[0] not in "EC" or inst not in bounds:
                continue
            hk = bounds[inst][0]
            kind = "uniform" if inst[0] == "E" else "clustered"
            gaps[(r["algorithm"], kind, int(r["n"]))].append((int(r["tour_length"]) / hk - 1) * 100)

    mean = lambda v: sum(v) / len(v) if v else None
    seeds = {n: len(gaps[("greedy", "uniform", n)]) for _, n in SIZES}
    xmax = max(mean(v) for v in gaps.values() if v) * 1.22      # room for the labels

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED})
    fig, axes = plt.subplots(1, len(SIZES), figsize=(15, 5.9), facecolor=SURFACE, sharey=True)

    for ax, (label, n) in zip(axes, SIZES):
        ax.set_facecolor(SURFACE)
        ax.grid(True, axis="x", color=GRID, linewidth=1)
        ax.set_axisbelow(True)
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        ax.tick_params(axis="y", length=0)
        ax.set_xlim(0, xmax)
        ax.set_xlabel("% above Held-Karp bound (lower is better)", color=MUTED)
        ax.set_title(f"{label} cities  ·  {seeds[n]} seeds each", loc="left",
                     color=TEXT, fontsize=11.5, fontweight="bold")

        for row, (key, name) in enumerate(SOLVERS):
            u = mean(gaps[(key, "uniform", n)])
            c = mean(gaps[(key, "clustered", n)])
            if u is None or c is None:
                continue
            y = row
            ax.plot([u, c], [y, y], color=GRID, lw=5, solid_capstyle="round", zorder=1)
            ax.scatter([u], [y], s=70, color=UNIFORM, zorder=3, edgecolor=SURFACE, linewidth=1.5)
            ax.scatter([c], [y], s=70, color=CLUSTERED, zorder=3, edgecolor=SURFACE, linewidth=1.5)
            diff = c - u
            ax.annotate(f"{diff:+.1f} pts", xy=(max(u, c), y), xytext=(9, 0),
                        textcoords="offset points", va="center", fontsize=9.5,
                        color=TEXT, fontweight="bold")

        ax.set_yticks(range(len(SOLVERS)), [name for _, name in SOLVERS])
        ax.invert_yaxis()

    for lbl in axes[0].get_yticklabels():
        lbl.set_color(TEXT)
        lbl.set_fontsize(10.5)

    handles = [plt.Line2D([], [], marker="o", ls="", ms=8, color=UNIFORM, label="uniform"),
               plt.Line2D([], [], marker="o", ls="", ms=8, color=CLUSTERED, label="clustered"),
               plt.Line2D([], [], ls="", label="label = clustered minus uniform, in percentage points")]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.15, 0.0), ncol=3,
               frameon=False, fontsize=10, labelcolor=TEXT)
    fig.suptitle("Clustered vs uniform: how much each solver loses on clumpy data",
                 x=0.15, ha="left", color=TEXT, fontsize=14, fontweight="bold", y=0.985)
    fig.text(0.15, 0.885, "Mean % above the Held-Karp lower bound, the same reference for both "
             "data types. Only seeds with a bound for both types are used, so 1M is not shown.\n"
             + CODE_KEY + ".", color=MUTED, fontsize=9.5)
    fig.subplots_adjust(left=0.15, right=0.975, top=0.78, bottom=0.2, wspace=0.18)

    out = ROOT / "results" / "2d" / "clustered_vs_uniform.png"
    fig.savefig(out, dpi=170, facecolor=SURFACE)
    print(f"Wrote {out}")
    for key, name in SOLVERS:
        print(f"  {name:11}", "  ".join(
            f"{lab}: {mean(gaps[(key, 'uniform', n)]):5.1f} -> {mean(gaps[(key, 'clustered', n)]):5.1f}"
            for lab, n in SIZES))


if __name__ == "__main__":
    main()
