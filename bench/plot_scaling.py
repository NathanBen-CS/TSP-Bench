"""Charts for the 1-million-city experiment: every solver except nn.

    py bench/plot_scaling.py      # writes the three PNGs below

Every chart uses one measure, "% longer than greedy TO Rn": how much longer a
solver's tour is than the greedy TO Rn tour on the same map (greedy start,
then 2-opt + Or-opt; the best solver in the bench).  0 = as good as greedy
TO, lower is better.  It is used because the 1M maps we generated have no
Held-Karp bound, so a gap to Held-Karp can't be given there.

  results/scale_1M.png        one row per solver, at 1,000,000 cities:
                              quality (left) and seconds (right),
                              blue = uniform maps, orange = clustered maps
  results/scale_growth.png    one small panel per cleanup solver: its
                              quality at 1k, 10k, 100k and 1M cities
  results/tsplib_summary.png  one bar per solver on the TSPLIB maps
                              (% above Held-Karp, which TSPLIB maps all have)

Scaling data: the uniform + clustered set (results/box2/runs.csv; 10 maps of
each type at 1k, 3 at 10k, 100k and 1M).  TSPLIB data: results/runs.csv.
A value is the mean over the maps of that size.

"hilbert B M1 t=10" (the box tour with no cleanup) is not a separate run: it
is the reference box run's start_length, and its time is that run's user
time minus the 2-opt and Or-opt times it reports.

nn is left out on purpose: it is O(N^2) and does not reach 1M.
"""

import csv
import re
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from common import CODE_KEY, ROOT, config_label, load_bounds

# (algorithm, level); level "none" = the box tour before cleanup
CLEANED = [("twoopt_or", "full"), ("twoopt", "full"),
           ("box:hilbert:M5:t10", "full"), ("box:sierpinski:M5:t4", "full"),
           ("box:sierpinski:M1:t10", "full"),
           ("box:knopp:M1:t10", "full"), ("box:sierpinski:M2:t12", "full"),
           ("box:hilbert:M1:t10", "full")]
ORDER_ONLY = [("greedy", "full"), ("box:hilbert:M1:t10", "none"), ("curve:sierpinski", "full"),
              ("curve:moore", "full"), ("spacefill", "full")]
REF = ("greedy+2opt+oropt", "full")             # greedy TO Rn
TYPES = [("uniform", "E", "#2a78d6"), ("clustered", "C", "#eb6834")]
SIZES = [("1k", 1e3), ("10k", 1e4), ("100k", 1e5), ("1M", 1e6)]
SURFACE, TEXT, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e0"


def load(path):
    """runs[(alg, level)][instance] = (length, user seconds), OK runs only."""
    runs = defaultdict(dict)
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            if r["status"] != "OK":
                continue
            ex = dict(kv.split("=", 1) for kv in r["extras"].split(";") if "=" in kv)
            runs[(r["algorithm"], "full")][r["instance"]] = (int(r["tour_length"]), float(r["user_s"]))
            if r["algorithm"].startswith("box:") and "start_length" in ex:
                t = float(r["user_s"]) - float(ex.get("twoopt_time", 0)) - float(ex.get("oropt_time", 0))
                runs[(r["algorithm"], "none")][r["instance"]] = (float(ex["start_length"]), t)
    return runs


def maps(runs, prefix, size):
    """The maps of one type and size, e.g. E1M.0, E1M.1, E1M.2."""
    return sorted(i for i in runs[REF] if re.fullmatch(f"{prefix}{size}\\.\\d+", i))


def mean_vs_ref(runs, key, insts):
    """(mean % longer than greedy TO Rn, mean seconds), or None if a map is missing."""
    if not insts or not all(i in runs[key] for i in insts):
        return None
    pct = [100 * (runs[key][i][0] / runs[REF][i][0] - 1) for i in insts]
    sec = [runs[key][i][1] for i in insts]
    return sum(pct) / len(pct), sum(sec) / len(sec)


def plain_axes(ax, grid_axis):
    ax.set_facecolor(SURFACE)
    ax.grid(True, axis=grid_axis, color=GRID, linewidth=1)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def setup():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10.5,
                         "axes.edgecolor": MUTED, "axes.labelcolor": MUTED,
                         "xtick.color": MUTED, "ytick.color": MUTED})


def chart_1m(runs, out):
    """Two panels, one row per solver: % longer than greedy TO Rn, and seconds."""
    rows = [REF] + CLEANED + ORDER_ONLY
    vals = {(k, kind): mean_vs_ref(runs, k, maps(runs, p, "1M")) for k in rows for kind, p, _ in TYPES}
    rows = [k for k in rows if any(vals[(k, kind)] for kind, _, _ in TYPES)]

    fig, (axq, axt) = plt.subplots(1, 2, figsize=(15, 0.62 * len(rows) + 2.6), facecolor=SURFACE,
                                   gridspec_kw={"width_ratios": [1.15, 1]})
    h = 0.36                                             # bar height; two bars per row
    for ax, col, xlabel, fmt in ((axq, 0, "% longer than greedy TO Rn (lower is better)", "{:+.1f}%"),
                                 (axt, 1, "user CPU seconds (lower is faster)", "{:.1f} s")):
        plain_axes(ax, "x")
        top = max(v[col] for v in vals.values() if v)
        for y, key in enumerate(rows):
            for j, (kind, _, colour) in enumerate(TYPES):
                v = vals[(key, kind)]
                yy = y - h / 2 + j * h
                if v is None:
                    ax.text(0, yy, "  not run", va="center", color=MUTED, fontsize=8.5)
                    continue
                ax.barh(yy, v[col], height=h * 0.88, color=colour)
                ax.text(max(v[col], 0) + top * 0.01, yy, fmt.format(v[col]), va="center",
                        color=TEXT, fontsize=8.5)
        ax.set_xlim(0, top * 1.14)
        ax.set_xlabel(xlabel)
        ax.set_ylim(len(rows) - 0.4, -0.6)
        # a line between the cleanup solvers and the order-only ones
        split = 1 + len([k for k in CLEANED if k in rows])
        ax.axhline(split - 0.5, color=MUTED, lw=1)
    axq.set_yticks(range(len(rows)), [config_label(*k) + ("  (= 0)" if k == REF else "") for k in rows],
                   color=TEXT)
    axt.set_yticks(range(len(rows)), [""] * len(rows))
    axq.annotate("with cleanup", xy=(1, 0), xycoords=("axes fraction", "data"), xytext=(-4, 0),
                 textcoords="offset points", ha="right", va="center", color=MUTED, fontsize=9.5,
                 fontweight="bold")
    axq.annotate("order only (no cleanup)", xy=(1, split), xycoords=("axes fraction", "data"),
                 xytext=(-4, 0), textcoords="offset points", ha="right", va="center", color=MUTED,
                 fontsize=9.5, fontweight="bold")
    axq.set_title("How good is the tour?", loc="left", color=TEXT, fontsize=12.5, fontweight="bold")
    axt.set_title("How long does it take?", loc="left", color=TEXT, fontsize=12.5, fontweight="bold")

    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for _, _, c in TYPES]
    fig.legend(handles, [f"{k} maps" for k, _, _ in TYPES], loc="upper right", ncol=2,
               frameon=False, bbox_to_anchor=(0.985, 0.985), fontsize=10.5, labelcolor=TEXT)
    fig.suptitle("Every solver at 1,000,000 cities", x=0.015, ha="left", color=TEXT,
                 fontsize=15, fontweight="bold", y=0.985)
    fig.text(0.015, 1 - 1.15 / fig.get_figheight(),
             "Mean over 3 uniform and 3 clustered maps of 1M cities. Quality is how much longer the "
             "tour is than greedy TO Rn's tour on the same map\n(greedy TO Rn is the best solver, so it is 0). "
             + CODE_KEY + ".", color=MUTED, fontsize=9.5)
    fig.subplots_adjust(left=0.17, right=0.97, top=1 - 1.75 / fig.get_figheight(),
                        bottom=0.75 / fig.get_figheight(), wspace=0.06)
    fig.savefig(out, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    print(f"Wrote {out.relative_to(ROOT)}")
    for key in rows:
        cells = "   ".join(f"{kind}: " + (f"{v[0]:+6.2f}% {v[1]:7.2f}s" if v else "   not run")
                           for kind, _, _ in TYPES for v in [vals[(key, kind)]])
        print(f"  {config_label(*key):24} {cells}")


def chart_growth(runs, out):
    """One small panel per cleanup solver: quality at each size, two lines."""
    cols = 3
    nrows = -(-(len(CLEANED) + 1) // cols)
    fig, axes = plt.subplots(nrows, cols, figsize=(15, 4.1 * nrows + 1.6), facecolor=SURFACE,
                             sharey=True)
    axes = list(axes.flat)
    xs_all = [x for _, x in SIZES]
    hi = 0
    for ax, key in zip(axes, CLEANED):
        plain_axes(ax, "y")
        ax.axhline(0, color=MUTED, lw=1)
        ends = []
        for kind, p, colour in TYPES:
            pts = [(x, v[0]) for s, x in SIZES for v in [mean_vs_ref(runs, key, maps(runs, p, s))] if v]
            if not pts:
                continue
            ax.plot([x for x, _ in pts], [y for _, y in pts], color=colour, lw=2.2, marker="o", ms=8,
                    markeredgecolor=SURFACE, markeredgewidth=1.5)
            ends.append(pts[-1])
            hi = max(hi, max(y for _, y in pts))
        # label the last value of each line; push the two labels apart if they are close
        ends.sort(key=lambda e: e[1])
        close = len(ends) == 2 and abs(ends[1][1] - ends[0][1]) < 0.5
        for j, (x, y) in enumerate(ends):
            ax.annotate(f"{y:+.1f}%", xy=(x, y), xytext=(7, (j * 2 - 1) * 6 if close else 0),
                        textcoords="offset points", va="center", color=TEXT, fontsize=9)
        ax.set_xscale("log")
        ax.set_xticks(xs_all, [s for s, _ in SIZES])
        ax.minorticks_off()
        ax.set_xlim(6e2, 2.6e6)
        ax.set_title(config_label(*key), loc="left", color=TEXT, fontsize=12, fontweight="bold")
    for ax in axes[len(CLEANED):]:
        ax.axis("off")
    last = axes[len(CLEANED)]
    last.text(0.02, 0.95, "How to read a panel:\n\n"
              "• x: number of cities in the map\n"
              "• y: % longer than greedy TO Rn on the same\n   map; the grey line at 0 is greedy TO Rn\n"
              "• blue: uniform maps, orange: clustered\n"
              "• the label is the value at 1M\n\n"
              "A flat line means the solver keeps the\nsame quality as the maps grow.",
              transform=last.transAxes, va="top", color=MUTED, fontsize=10.5)
    for ax in fig.axes[::cols]:
        ax.set_ylabel("% longer than greedy TO Rn")
    fig.axes[0].set_ylim(-0.6, hi * 1.12)

    handles = [plt.Line2D([], [], color=c, lw=2.2, marker="o", ms=8) for _, _, c in TYPES]
    fig.legend(handles, [f"{k} maps" for k, _, _ in TYPES], loc="upper right", ncol=2,
               frameon=False, bbox_to_anchor=(0.985, 0.99), fontsize=10.5, labelcolor=TEXT)
    fig.suptitle("Does each solver stay as good as the maps grow?", x=0.04, ha="left", color=TEXT,
                 fontsize=15, fontweight="bold", y=0.99)
    fig.text(0.04, 1 - 1.0 / fig.get_figheight(),
             "Solvers with cleanup only (the order-only solvers are 10-50% longer at every size). "
             "Mean over 10 maps at 1k and 3 at 10k, 100k and 1M.\n" + CODE_KEY + ".",
             color=MUTED, fontsize=9.5)
    fig.subplots_adjust(left=0.06, right=0.98, top=1 - 1.5 / fig.get_figheight(), bottom=0.06,
                        hspace=0.38, wspace=0.12)
    fig.savefig(out, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    print(f"Wrote {out.relative_to(ROOT)}")


def chart_tsplib(runs, bounds, out):
    """Mean % above Held-Karp on the TSPLIB maps every solver has finished."""
    series = [REF] + CLEANED + ORDER_ONLY
    tsplib = [i for i in runs[REF] if not re.match(r"[ECM]\d", i) and i in bounds]
    common = [i for i in tsplib if all(i in runs[k] for k in series)]
    if not common:
        print("No TSPLIB map has every solver finished yet; skipped", out.name)
        return
    rows = sorted((sum(100 * (runs[k][i][0] / bounds[i][0] - 1) for i in common) / len(common), k)
                  for k in series)

    fig, ax = plt.subplots(figsize=(11, 0.45 * len(rows) + 2), facecolor=SURFACE)
    plain_axes(ax, "x")
    for y, (g, key) in enumerate(rows):
        ax.barh(y, g, height=0.62, color="#2a78d6" if key in [REF] + CLEANED else "#9ec5f4")
        ax.text(g + 0.4, y, f"{g:.1f}%", va="center", color=TEXT, fontsize=9.5)
    ax.set_yticks(range(len(rows)), [config_label(*k) for _, k in rows], color=TEXT)
    ax.set_ylim(len(rows) - 0.4, -0.6)
    ax.set_xlim(0, rows[-1][0] * 1.1)
    ax.set_xlabel("mean % above the Held-Karp lower bound (lower is better)")
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in ("#2a78d6", "#9ec5f4")]
    fig.legend(handles, ["with cleanup", "order only (no cleanup)"], loc="upper right", ncol=2,
               frameon=False, bbox_to_anchor=(0.97, 0.985), labelcolor=TEXT)
    fig.suptitle("TSPLIB maps: every solver except nn", x=0.02, ha="left", color=TEXT,
                 fontsize=14, fontweight="bold", y=0.985)
    fig.text(0.02, 1 - 0.95 / fig.get_figheight(),
             f"Mean over the {len(common)} geometric TSPLIB maps (1,000 to 85,900 cities).\n"
             f"{CODE_KEY}.", color=MUTED, fontsize=9.5)
    fig.subplots_adjust(left=0.22, right=0.97, top=1 - 1.5 / fig.get_figheight(),
                        bottom=0.75 / fig.get_figheight())
    fig.savefig(out, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    print(f"Wrote {out.relative_to(ROOT)}")
    for g, key in rows:
        print(f"  {config_label(*key):24} {g:6.2f}%")


def main():
    setup()
    box = load(ROOT / "results" / "box2" / "runs.csv")
    chart_1m(box, ROOT / "results" / "scale_1M.png")
    chart_growth(box, ROOT / "results" / "scale_growth.png")
    chart_tsplib(load(ROOT / "results" / "runs.csv"), load_bounds(ROOT / "refs" / "bounds"),
                 ROOT / "results" / "tsplib_summary.png")


if __name__ == "__main__":
    main()
