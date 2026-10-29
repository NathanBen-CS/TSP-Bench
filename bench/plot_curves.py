"""Pictures of every space-filling curve the box solver can follow.

    py bench/plot_curves.py          # writes results/curves/<curve>.png and all_curves.png

For each curve, one PNG with three panels:
  1. the curve itself: a small grid of points visited in curve order
  2. its path through 1,000 uniform cities (E1k.0)
  3. its path through 1,000 clustered cities (C1k.0)
The line is shaded light to dark from start to end so the order can be
followed; a closed curve's last step back to the start is drawn dashed.

Every path comes from the real solver (solvers/boxtsp.exe --start curve
--curve NAME --cleanup none), so the pictures show exactly the order the
benchmark used; nothing is re-implemented here.  Following the naming used
in the charts, a bare curve name means "that order only" (no boxes, no
cleanup).
"""

import subprocess
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap

from common import ROOT, load_bounds
from plot_tour import read_instance, read_tour, tour_length

CURVES = [  # name, grid points per side for the picture of the curve, closed?, description
    ("hilbert",    16, False, "squares, 4 children per cell"),
    ("moore",      16, True,  "four Hilbert curves joined into a loop"),
    ("sierpinski", 32, True,  "triangles, 2 children per cell; closed"),
    ("knopp",      32, False, "one triangle covering the square"),
    ("peano",      27, False, "3 x 3 squares, 9 children per cell"),
    ("morton",     16, False, "Z-order: jumps between cells (negative control)"),
    ("strip",      16, False, "snake rows"),
]
INSTANCES = [("E1k.0", "1,000 uniform cities"), ("C1k.0", "1,000 clustered cities")]
BOX_DIR = Path("C:/Users/Natha/tsp-instances-box")
SURFACE, TEXT, MUTED = "#fcfcfb", "#0b0b0b", "#52514e"
START, END = "#1baf7a", "#e34948"
SHADE = LinearSegmentedColormap.from_list("progress", ["#b8d4f5", "#2a78d6", "#0d2f5c"])
OUT = ROOT / "results" / "curves"


def grid_instance(path, side, nudge):
    """A TSPLIB file with one city at every point of a side x side grid.

    nudge: move each point by a fixed pseudo-random amount, up to 30% of the
    spacing.  Used for the triangle curves: on a perfect grid whole rows and
    diagonals lie exactly on the curves' cut lines and all hit the same tie
    rule, which real cities almost never do.                              """
    with open(path, "w") as f:
        f.write(f"NAME : grid{side}\nTYPE : TSP\nDIMENSION : {side * side}\n"
                "EDGE_WEIGHT_TYPE : EUC_2D\nNODE_COORD_SECTION\n")
        for y in range(side):
            for x in range(side):
                dx = dy = 0
                if nudge:
                    h = (x * 73856093) ^ (y * 19349663)           # fixed, repeatable
                    dx = (h % 601) - 300
                    dy = ((h // 601) % 601) - 300
                f.write(f"{y * side + x + 1} {x * 1000 + dx} {y * 1000 + dy}\n")


def curve_order(instance, curve, workdir):
    """Run the solver to get the curve's visiting order on an instance."""
    tour_file = Path(workdir) / f"{curve}.tour"
    subprocess.run([str(ROOT / "solvers" / "boxtsp.exe"), str(instance), str(tour_file),
                    "--start", "curve", "--curve", curve, "--cleanup", "none"],
                   check=True, capture_output=True)
    return tour_file


def draw_path(ax, coords, tour, closed, lw, dots):
    pts = [coords[c] for c in tour]
    segs = [[pts[i], pts[i + 1]] for i in range(len(pts) - 1)]
    lc = LineCollection(segs, cmap=SHADE, linewidths=lw, capstyle="round", joinstyle="round")
    lc.set_array(list(range(len(segs))))                 # shade by position along the path
    ax.add_collection(lc)
    ax.plot([pts[-1][0], pts[0][0]], [pts[-1][1], pts[0][1]],   # back to the start
            color=MUTED, lw=lw * 0.8, ls=(0, (3, 3)) if closed else (0, (1, 3)), alpha=0.8)
    if dots:
        ax.scatter([p[0] for p in pts], [p[1] for p in pts], s=dots, color=MUTED, zorder=3, linewidths=0)
    ax.scatter([pts[0][0]], [pts[0][1]], s=60, color=START, zorder=4, edgecolor=SURFACE, linewidth=1.2)
    ax.scatter([pts[-1][0]], [pts[-1][1]], s=60, color=END, zorder=4, edgecolor=SURFACE, linewidth=1.2)
    ax.set_aspect("equal")
    ax.autoscale()
    ax.axis("off")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    bounds = load_bounds(ROOT / "refs" / "bounds")
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10})
    overview = []

    with tempfile.TemporaryDirectory() as work:
        for name, side, closed, desc in CURVES:
            fig, axes = plt.subplots(1, 3, figsize=(16, 5.9), facecolor=SURFACE)

            # 1. the curve on a grid
            nudge = name in ("sierpinski", "knopp")       # triangle curves
            grid = Path(work) / f"grid{side}{'n' if nudge else ''}.tsp"
            if not grid.exists():
                grid_instance(grid, side, nudge)
            _, _, gcoords = read_instance(grid)
            gtour = read_tour(curve_order(grid, name, work), len(gcoords))
            draw_path(axes[0], gcoords, gtour, closed, lw=2.2, dots=6)
            axes[0].set_title(f"{name} on a {side} x {side} grid", loc="left", color=TEXT,
                              fontsize=11.5, fontweight="bold")
            overview.append((name, side, closed, gcoords, gtour))

            # 2-3. the curve's path through real cities
            for ax, (inst, what) in zip(axes[1:], INSTANCES):
                path = BOX_DIR / inst
                _, ewt, coords = read_instance(path)
                tour = read_tour(curve_order(path, name, work), len(coords))
                L = tour_length(tour, coords, ewt)
                hk = bounds.get(inst, (None, None))[0]
                draw_path(ax, coords, tour, closed, lw=0.9, dots=3)
                gap = f"  ·  {100 * (L / hk - 1):.1f}% above Held-Karp" if hk else ""
                ax.set_title(f"{name} on {inst} ({what})\nlength {L:,}{gap}", loc="left",
                             color=TEXT, fontsize=10.5, fontweight="bold")

            fig.suptitle(f"{name}: {desc}", x=0.03, ha="left", color=TEXT, fontsize=14,
                         fontweight="bold", y=0.985)
            fig.text(0.03, 0.915, "Line shaded light to dark from start (green) to end (red). "
                     f"Dashed grey: the step back to the start{' (closed curve: the ends are neighbours)' if closed else ''}. "
                     "Order from solvers/boxtsp.exe --start curve, no boxes, no cleanup.",
                     color=MUTED, fontsize=9.5)
            fig.subplots_adjust(left=0.03, right=0.98, top=0.8, bottom=0.03, wspace=0.12)
            out = OUT / f"{name}.png"
            fig.savefig(out, dpi=150, facecolor=SURFACE)
            plt.close(fig)
            print(f"Wrote {out.relative_to(ROOT)}")

    # all curves side by side
    fig, axes = plt.subplots(2, 4, figsize=(17, 9.4), facecolor=SURFACE)
    for ax, (name, side, closed, coords, tour) in zip(axes.flat, overview):
        draw_path(ax, coords, tour, closed, lw=1.8, dots=0)
        ax.set_title(f"{name}{'  (closed)' if closed else ''}", loc="left", color=TEXT,
                     fontsize=12, fontweight="bold")
    axes.flat[-1].axis("off")
    axes.flat[-1].text(0.02, 0.6, "Light to dark: start (green) to end (red).\n"
                       "Dashed grey: the step back to the start.\n\n"
                       "Same grid sizes as the single pictures:\n16 x 16 for the square curves,\n"
                       "27 x 27 for Peano, 32 x 32 for the\ntriangle curves (Sierpinski, Knopp),\n"
                       "whose points are nudged slightly off\nthe exact grid (see the script).",
                       transform=axes.flat[-1].transAxes, color=MUTED, fontsize=10.5, va="top")
    fig.suptitle("Every curve the box solver can follow", x=0.03, ha="left", color=TEXT,
                 fontsize=15, fontweight="bold", y=0.985)
    fig.subplots_adjust(left=0.02, right=0.99, top=0.9, bottom=0.02, wspace=0.08, hspace=0.15)
    fig.savefig(OUT / "all_curves.png", dpi=150, facecolor=SURFACE)
    print(f"Wrote {(OUT / 'all_curves.png').relative_to(ROOT)}")


if __name__ == "__main__":
    main()
