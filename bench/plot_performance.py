"""One chart of solver performance: tour quality and running time against
instance size, on the uniform instances of the box experiment.

    py bench/plot_performance.py                 # writes results/solver_performance.png

Left: mean % above 0.7124 * sqrt(N * A) (the uniform-data estimate of the
optimal tour), lower is better.  Right: mean user CPU seconds, log scale.
Each point is the mean over the seeds at that size (10 at 1k, 3 otherwise).
"""

import csv
import math
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from common import CODE_KEY, ROOT, config_label

SOLVERS = [                          # (name in the CSV, colour, line style); labels
    ("spacefill",         "#e87ba4", "-"),   # come from config_label (common.py)
    ("boxes-only",        "#4a3aa7", "--"),  # same method as the BTO line,
    ("nn",                "#eda100", "-"),   #   before cleanup
    ("greedy",            "#eb6834", "-"),
    ("boxes+2opt+oropt",  "#4a3aa7", "-"),
    ("twoopt",            "#1baf7a", "-"),
    ("twoopt_or",         "#2a78d6", "-"),
]
SOLVERS = [(k, config_label(k), c, s) for k, c, s in SOLVERS]
SURFACE, TEXT, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e0"
AREA = 1e12
TIMER_FLOOR = 0.01            # Windows CPU time comes in ~0.016 s ticks


def main():
    gaps, times = defaultdict(list), defaultdict(list)
    with open(ROOT / "results" / "box" / "runs.csv", newline="") as f:
        for r in csv.DictReader(f):
            if r["status"] != "OK" or not r["instance"].startswith("E"):
                continue
            n = int(r["n"])
            ref = 0.7124 * math.sqrt(n * AREA)
            gaps[(r["algorithm"], n)].append((int(r["tour_length"]) / ref - 1) * 100)
            times[(r["algorithm"], n)].append(max(float(r["user_s"]), TIMER_FLOOR))

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.edgecolor": MUTED, "axes.labelcolor": MUTED,
                         "xtick.color": MUTED, "ytick.color": MUTED})
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16.5, 6.2), facecolor=SURFACE)

    for ax in (ax1, ax2):
        ax.set_facecolor(SURFACE)
        ax.set_xscale("log")
        ax.grid(True, which="major", color=GRID, linewidth=1)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.set_xlabel("cities (log scale)")
        ax.set_xticks([1e3, 1e4, 1e5, 1e6], ["1k", "10k", "100k", "1M"])
        ax.minorticks_off()

    ends1, ends2 = [], []
    for key, label, colour, style in SOLVERS:
        ns = sorted(n for (a, n) in gaps if a == key)
        if not ns:
            continue
        g = [sum(gaps[(key, n)]) / len(gaps[(key, n)]) for n in ns]
        t = [sum(times[(key, n)]) / len(times[(key, n)]) for n in ns]
        lw = 2.6 if key.startswith("boxes") else 2
        for ax, ys, ends in ((ax1, g, ends1), (ax2, t, ends2)):
            ax.plot(ns, ys, color=colour, lw=lw, ls=style, marker="o", ms=4.5,
                    solid_capstyle="round", solid_joinstyle="round", label=label)
            ends.append([ns[-1], ys[-1], label, colour])

    # direct labels at the line ends, nudged apart so they don't overlap
    def label_ends(ax, ends, fmt, log=False):
        tr = (lambda v: math.log10(v)) if log else (lambda v: v)
        lo, hi = ax.get_ylim()
        span = tr(hi) - tr(lo)
        ends.sort(key=lambda e: tr(e[1]))
        pos = [tr(e[1]) for e in ends]
        for i in range(1, len(pos)):
            pos[i] = max(pos[i], pos[i - 1] + span * 0.055)
        for (x, y, label, colour), p in zip(ends, pos):
            yy = 10 ** p if log else p
            ax.annotate(f"{label}  {fmt(y)}", xy=(x, y), xytext=(x * 1.35, yy),
                        color=colour, fontsize=9.5, va="center", fontweight="bold",
                        annotation_clip=False)

    ax1.set_ylim(0, max(e[1] for e in ends1) * 1.12)
    ax1.set_ylabel("% above optimal-tour estimate (lower is better)")
    ax1.set_title("Tour quality", loc="left", color=TEXT, fontsize=12, fontweight="bold")
    label_ends(ax1, ends1, lambda v: f"{v:.1f}%")

    ax2.set_yscale("log")
    ax2.set_ylim(TIMER_FLOOR * 0.8, max(e[1] for e in ends2) * 3)
    ax2.set_ylabel("user CPU seconds (log scale)")
    ax2.set_title("Running time", loc="left", color=TEXT, fontsize=12, fontweight="bold")
    label_ends(ax2, ends2, lambda v: f"{v:.2f} s" if v < 1 else f"{v:.1f} s", log=True)

    ax1.legend(loc="upper left", bbox_to_anchor=(0, -0.13), ncol=4, frameon=False,
               fontsize=9.5, labelcolor=TEXT)
    fig.suptitle("TSP solvers on uniform random instances", x=0.06, ha="left",
                 color=TEXT, fontsize=14, fontweight="bold", y=0.985)
    fig.text(0.06, 0.895, "Mean over seeds (10 at 1k, 3 at 10k–1M). nn is not run above "
             "320k cities. Runs faster than the 0.016 s timer tick sit on the floor.\n"
             + CODE_KEY + ".", color=MUTED, fontsize=9.5)
    fig.subplots_adjust(left=0.055, right=0.82, top=0.82, bottom=0.24, wspace=0.9)

    out = ROOT / "results" / "solver_performance.png"
    fig.savefig(out, dpi=170, facecolor=SURFACE)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
