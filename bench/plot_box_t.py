"""Chart: how the box size t affects the box method, alone and after cleanup.

    py bench/plot_box_t.py                   # writes results/box_t_sweep.png

Left panel: the box tour on its own (no cleanup).  Right panel: the same box
tour after 2-opt + Or-opt.  One line per data type and size, t = 6, 8, 10, 12.

Where the numbers come from (results/box/runs.csv):
  * right panel: the tour length of the full-method runs
    (boxes+2opt+oropt-t6, -t8, boxes+2opt+oropt (t = 10), -t12);
  * left panel: the same runs' start_length, the box tour before cleanup.
    It equals the boxes-only run's length on all 38 instances at t = 10,
    because the box tour does not depend on what happens after it.
Both are measured against the Held-Karp lower bound (refs/bounds), so the
two data types are comparable; seeds without a bound are left out (1M, and
the third 100k seed).  Both panels use the same number of percentage points
per unit of height, so how steep a line is compares directly between them.
"""

import csv
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from common import ROOT, load_bounds

T_RUNS = {"boxes+2opt+oropt-t6": 6, "boxes+2opt+oropt-t8": 8,
          "boxes+2opt+oropt": 10, "boxes+2opt+oropt-t12": 12}
TS = [6, 8, 10, 12]
CHOSEN_T = 10
TYPES = [("uniform", "E", "#2a78d6"), ("clustered", "C", "#eb6834")]
SIZES = [("1k", 1000, 0.45, "o"), ("10k", 10000, 0.72, "s"), ("100k", 100000, 1.0, "^")]
SURFACE, TEXT, MUTED, GRID, BAND = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e0", "#efeee9"


def main():
    bounds = load_bounds(ROOT / "refs" / "bounds")
    alone, cleaned = defaultdict(list), defaultdict(list)   # (type, n, t) -> [% above HK]
    with open(ROOT / "results" / "box" / "runs.csv", newline="") as f:
        for r in csv.DictReader(f):
            t = T_RUNS.get(r["algorithm"])
            inst = r["instance"]
            if t is None or r["status"] != "OK" or inst not in bounds:
                continue
            hk = bounds[inst][0]
            kind = "uniform" if inst[0] == "E" else "clustered"
            ex = dict(kv.split("=", 1) for kv in r["extras"].split(";") if "=" in kv)
            key = (kind, int(r["n"]), t)
            alone[key].append((float(ex["start_length"]) / hk - 1) * 100)
            cleaned[key].append((int(r["tour_length"]) / hk - 1) * 100)

    mean = lambda v: sum(v) / len(v)
    series = {}                                   # panel -> list of (label, colour, alpha, marker, ys)
    for panel, data in (("alone", alone), ("cleaned", cleaned)):
        series[panel] = []
        for kind, _, colour in TYPES:
            for size, n, alpha, marker in SIZES:
                if all((kind, n, t) in data for t in TS):
                    ys = [mean(data[(kind, n, t)]) for t in TS]
                    series[panel].append((f"{kind} {size}", colour, alpha, marker, ys))

    # same vertical span (percentage points) in both panels
    pad = 1.5
    ranges = {p: (min(min(s[4]) for s in series[p]), max(max(s[4]) for s in series[p])) for p in series}
    span = max(hi - lo for lo, hi in ranges.values()) + 2 * pad

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.edgecolor": MUTED, "axes.labelcolor": MUTED,
                         "xtick.color": MUTED, "ytick.color": MUTED})
    fig, axes = plt.subplots(1, 2, figsize=(15, 6.2), facecolor=SURFACE)
    titles = {"alone": "Box tour alone (no cleanup)",
              "cleaned": "Box tour + 2-opt + Or-opt"}

    for ax, panel in zip(axes, ("alone", "cleaned")):
        ax.set_facecolor(SURFACE)
        lo = ranges[panel][0] - pad
        ax.set_ylim(lo, lo + span)
        ax.set_xlim(5.4, 12.6)
        ax.axvspan(CHOSEN_T - 0.45, CHOSEN_T + 0.45, color=BAND, zorder=0)
        ax.text(CHOSEN_T, lo + span, "t = 10, used in\nthe other charts", ha="center", va="top",
                fontsize=8.5, color=MUTED)
        ax.grid(True, axis="y", color=GRID, linewidth=1)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.set_xticks(TS, [f"t = {t}" for t in TS])
        ax.set_xlabel("box size t (most cities allowed in one box)")
        ax.set_ylabel("% above Held-Karp bound (lower is better)")
        ax.set_title(titles[panel], loc="left", color=TEXT, fontsize=12, fontweight="bold")

        ends = []
        for label, colour, alpha, marker, ys in series[panel]:
            ax.plot(TS, ys, color=colour, alpha=alpha, lw=2.2, marker=marker, ms=6,
                    solid_capstyle="round", solid_joinstyle="round")
            change = round(ys[-1] - ys[0], 1) + 0.0     # + 0.0 turns -0.0 into 0.0
            ends.append([ys[-1], f"{label}  {ys[-1]:.1f}%  ({change:+.1f})", colour, alpha])

        # end labels, nudged apart so they don't overlap
        ends.sort(key=lambda e: e[0])
        pos = [e[0] for e in ends]
        gap = span * 0.052
        for i in range(1, len(pos)):
            pos[i] = max(pos[i], pos[i - 1] + gap)
        for (y, text, colour, alpha), p in zip(ends, pos):
            ax.annotate(text, xy=(TS[-1], y), xytext=(TS[-1] + 0.35, p), color=colour,
                        alpha=max(alpha, 0.7), fontsize=9, va="center", fontweight="bold",
                        annotation_clip=False)

    handles = [plt.Line2D([], [], color=c, lw=2.2, label=k) for k, _, c in TYPES]
    handles += [plt.Line2D([], [], color=MUTED, alpha=a, marker=m, ls="", ms=6, label=s)
                for s, _, a, m in SIZES]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.055, 0.0), ncol=5,
               frameon=False, fontsize=10, labelcolor=TEXT)
    fig.suptitle("Box size t: a small gain on the box tour, and none left after cleanup",
                 x=0.055, ha="left", color=TEXT, fontsize=14, fontweight="bold", y=0.985)
    fig.text(0.055, 0.875, "Mean % above the Held-Karp bound (10 seeds at 1k, 3 at 10k, 2 at 100k). "
             "Both panels use the same vertical scale, so slopes compare directly.\n"
             "Labels: value at t = 12, and in brackets the change since t = 6 (negative = better).", color=MUTED, fontsize=9.5)
    fig.subplots_adjust(left=0.055, right=0.87, top=0.82, bottom=0.2, wspace=0.55)

    out = ROOT / "results" / "box_t_sweep.png"
    fig.savefig(out, dpi=170, facecolor=SURFACE)
    print(f"Wrote {out}")
    for panel in ("alone", "cleaned"):
        print(" ", panel)
        for label, _, _, _, ys in series[panel]:
            print(f"    {label:16}", "  ".join(f"t={t}: {y:5.2f}" for t, y in zip(TS, ys)))


if __name__ == "__main__":
    main()
