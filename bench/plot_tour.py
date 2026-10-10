"""Draw a TSP tour on top of its cities.

    py bench/plot_tour.py <instance.tsp> <tour_file>            # opens a window
    py bench/plot_tour.py <instance.tsp> <tour_file> -o tour.png # saves an image

<instance.tsp>  a TSPLIB file with coordinates (EUC_2D or CEIL_2D), e.g. E1k.0
                or pr1002.tsp.  Matrix instances (M*.k) have no coordinates,
                so there is nothing to draw.
<tour_file>     the DIMACS solver output: city numbers 1..N, one per line,
                in visiting order (the return to the first city is implied).

The tour length in the title uses the same distance rules as codes/length.exe
(EUC_2D: nearest integer, CEIL_2D: rounded up), so the two must agree.
Standalone: only needs matplotlib.  Best for small/medium tours (<= ~20k cities).
"""

import argparse
import math
import sys
from pathlib import Path

# colours: reference palette (dataviz skill), light mode
SURFACE = "#fcfcfb"
TEXT = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
TOUR = "#2a78d6"     # series slot 1, blue
START = "#eb6834"    # series slot 2, orange


def read_instance(path):
    """Return (name, edge_weight_type, coords) where coords[i] = (x, y), i = 1..N."""
    header, n, coords = {}, None, {}
    with open(path) as f:
        lines = iter(f)
        for line in lines:
            line = line.strip()
            if not line:
                continue
            if line.startswith("NODE_COORD_SECTION"):
                for _ in range(n):
                    i, x, y = next(lines).split()[:3]
                    coords[int(i)] = (float(x), float(y))
                break
            if line.startswith("EDGE_WEIGHT_SECTION") or line.startswith("EOF"):
                break
            key, _, value = line.partition(":")
            if not _:                                  # "KEY VALUE" without colon
                key, _, value = line.partition(" ")
            header[key.strip()] = value.strip()
            if key.strip() == "DIMENSION":
                n = int(value)

    ewt = header.get("EDGE_WEIGHT_TYPE", "")
    if not coords:
        sys.exit(f"{path} has no coordinates (EDGE_WEIGHT_TYPE {ewt}); "
                 "matrix instances can't be drawn.")
    if ewt not in ("EUC_2D", "CEIL_2D"):
        sys.exit(f"Unsupported EDGE_WEIGHT_TYPE {ewt} (expected EUC_2D or CEIL_2D).")
    if len(coords) != n:
        sys.exit(f"DIMENSION says {n} cities but {len(coords)} coordinates were read.")
    return header.get("NAME", Path(path).name), ewt, coords


def read_tour(path, n):
    """Read and validate a DIMACS tour: exactly N numbers, each city 1..N once."""
    try:
        tour = [int(t) for t in Path(path).read_text().split()]
    except ValueError as e:
        sys.exit(f"{path}: not a list of city numbers ({e}).")
    if len(tour) != n:
        sys.exit(f"{path}: has {len(tour)} cities, instance has {n}.")
    seen = set()
    for c in tour:
        if not 1 <= c <= n:
            sys.exit(f"{path}: city {c} is out of range 1..{n}.")
        if c in seen:
            sys.exit(f"{path}: city {c} appears twice.")
        seen.add(c)
    return tour


def tour_length(tour, coords, ewt):
    rnd = math.ceil if ewt == "CEIL_2D" else (lambda d: int(d + 0.5))  # as Instance.c
    total = 0
    for a, b in zip(tour, tour[1:] + tour[:1]):
        (xa, ya), (xb, yb) = coords[a], coords[b]
        total += rnd(math.hypot(xa - xb, ya - yb))
    return total


def plot(name, coords, tour, length, tour_path, out):
    import matplotlib
    if out:
        matplotlib.use("Agg")          # no window needed when saving
    import matplotlib.pyplot as plt

    n = len(tour)
    xs = [coords[c][0] for c in tour + tour[:1]]   # close the loop
    ys = [coords[c][1] for c in tour + tour[:1]]

    line_w = 1.5 if n <= 1000 else 0.8 if n <= 10000 else 0.4
    dot = 5 if n <= 200 else 2.5 if n <= 2000 else 0

    # size the figure to the map's shape so wide/tall instances fill the image
    w = (max(xs) - min(xs)) or 1
    h = (max(ys) - min(ys)) or 1
    fig_w = 10
    fig_h = min(max(fig_w * h / w, 4), 13) + 0.9       # + room for the title
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    ax.plot(xs, ys, color=TOUR, linewidth=line_w, solid_joinstyle="round", zorder=1)
    if dot:
        ax.plot(xs[:-1], ys[:-1], "o", color=TOUR, markersize=dot,
                markeredgecolor=SURFACE, markeredgewidth=0.8, zorder=2)

    # start city
    x0, y0 = coords[tour[0]]
    ax.plot([x0], [y0], "o", color=START, markersize=max(dot * 1.8, 8),
            markeredgecolor=SURFACE, markeredgewidth=1.5, zorder=4)
    if n <= 60:                         # city numbers only when they stay readable
        for c in tour[1:]:
            ax.annotate(str(c), coords[c], xytext=(4, 4), textcoords="offset points",
                        color=TEXT_SECONDARY, fontsize=7)
        ax.annotate(f"{tour[0]} · start", (x0, y0), xytext=(6, -6),
                    textcoords="offset points", va="top", color=TEXT, fontsize=8,
                    fontweight="bold", zorder=5)
    else:
        right_half = x0 > (min(xs) + max(xs)) / 2       # label toward the open side
        top_half = y0 > (min(ys) + max(ys)) / 2
        ax.annotate(f"start (city {tour[0]})", (x0, y0),
                    xytext=(-8 if right_half else 8, -8 if top_half else 8),
                    textcoords="offset points", ha="right" if right_half else "left",
                    va="top" if top_half else "bottom", color=TEXT, fontsize=9, zorder=5)

    ax.set_aspect("equal")
    ax.axis("off")
    top = 1 - 0.85 / fig_h                              # title band: fixed 0.85 in
    fig.text(0.02, 1 - 0.2 / fig_h, name, color=TEXT, fontsize=15,
             fontweight="bold", va="top")
    fig.text(0.02, 1 - 0.5 / fig_h, f"{Path(tour_path).name}  ·  {n:,} cities  ·  "
             f"tour length {length:,}", color=TEXT_SECONDARY, fontsize=10, va="top")
    fig.subplots_adjust(left=0.02, right=0.98, bottom=0.02, top=top)

    if out:
        fig.savefig(out, dpi=200, facecolor=SURFACE)
        print(f"Saved {out}")
    else:
        plt.show()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("instance", help="TSPLIB instance file (with coordinates)")
    ap.add_argument("tour", help="tour file: city numbers 1..N, one per line")
    ap.add_argument("-o", "--out", help="save to this image file (.png, .svg, .pdf) "
                                        "instead of opening a window")
    args = ap.parse_args()

    name, ewt, coords = read_instance(args.instance)
    tour = read_tour(args.tour, len(coords))
    length = tour_length(tour, coords, ewt)
    print(f"{name}: {len(tour)} cities, tour length {length}")
    if len(tour) > 20000:
        print("Large tour: drawing may be slow and the picture dense.")
    plot(name, coords, tour, length, args.tour, args.out)


if __name__ == "__main__":
    main()
