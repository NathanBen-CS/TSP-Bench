"""Rank every algorithm on the real-world TSPLIB maps, by type of problem.

    py bench/real_world.py      # writes results/2d/real_world.md and results/2d/real_world.png

The 34 maps in big/ (TSPLIB) come from real applications.  MAPS below
records what each one is, from the COMMENT line in its own file header
(plus the University of Waterloo TSP site for rl5915, whose header doesn't say).
They are grouped into types of problem: cities and towns, circuit-board
drilling, chip layout, origin not stated, and two that are set apart
(dsj1000 is randomly generated; si1032 is a distance table with no
coordinates, so only 4 solvers can read it).

No new runs: every algorithm already ran on these maps in results/2d/runs.csv
(the bench used byte-identical copies in tsp-instances).  Every solver is
deterministic, so one run per map is the whole result.

Ranking: for each type, every algorithm's mean "% above optimal" over the
maps of that type.  Optimal = the proven optimal tour length in refs/bounds
(all 34 TSPLIB maps are solved), so 0% would be a perfect tour.  Lower
mean = better rank.  An algorithm is only ranked on a type if it finished
every map of that type.
"""

import csv
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from common import CODE_KEY, ROOT, config_label, load_bounds

# map -> (type, what the data is, where that comes from)
MAPS = {
    "usa13509": ("cities", "US towns with at least 500 people", "header"),
    "brd14051": ("cities", "West Germany, 1989 borders", "header"),
    "d15112":   ("cities", "Germany", "header"),
    "d18512":   ("cities", "Germany including the former East Germany", "header"),
    "fnl4461":  ("cities", "the five former East German states", "header"),
    "nrw1379":  ("cities", "1379 towns in North Rhine-Westphalia", "header"),
    "d1291":    ("drilling", "circuit-board drilling (Reinelt)", "header"),
    "d1655":    ("drilling", "circuit-board drilling (Reinelt)", "header"),
    "d2103":    ("drilling", "circuit-board drilling (Reinelt)", "header"),
    "fl1400":   ("drilling", "circuit-board drilling (Reinelt)", "header"),
    "fl1577":   ("drilling", "circuit-board drilling (Reinelt)", "header"),
    "fl3795":   ("drilling", "circuit-board drilling (Reinelt)", "header"),
    "pcb1173":  ("drilling", "circuit-board drilling (Juenger/Reinelt)", "header"),
    "pcb3038":  ("drilling", "circuit-board drilling (Juenger/Reinelt)", "header"),
    "u1060":    ("drilling", "circuit-board drilling (Reinelt)", "header"),
    "u1432":    ("drilling", "circuit-board drilling (Reinelt)", "header"),
    "u1817":    ("drilling", "circuit-board drilling (Reinelt)", "header"),
    "u2152":    ("drilling", "circuit-board drilling (Reinelt)", "header"),
    "u2319":    ("drilling", "circuit-board drilling (Reinelt)", "header"),
    "rl5915":   ("drilling", "circuit-board drilling, 5915 holes (Reinelt)", "Waterloo TSP site"),
    "pla7397":  ("chip", "programmed logic array layout (Johnson)", "header"),
    "pla33810": ("chip", "programmed logic array layout (Johnson)", "header"),
    "pla85900": ("chip", "programmed logic array layout (Johnson)", "header"),
    "rl1304":   ("unstated", "'1304-city TSP' (Reinelt)", "header"),
    "rl1323":   ("unstated", "'1323-city TSP' (Reinelt)", "header"),
    "rl1889":   ("unstated", "'1889-city TSP' (Reinelt)", "header"),
    "rl5934":   ("unstated", "'5934-city TSP' (Reinelt)", "header"),
    "rl11849":  ("unstated", "'11849-city TSP' (Reinelt)", "header"),
    "pr1002":   ("unstated", "'1002-city problem' (Padberg/Rinaldi)", "header"),
    "pr2392":   ("unstated", "'2392-city problem' (Padberg/Rinaldi)", "header"),
    "vm1084":   ("unstated", "'1084-city problem' (Reinelt)", "header"),
    "vm1748":   ("unstated", "'1784-city problem' (Reinelt; the header says 1784)", "header"),
    "dsj1000":  ("synthetic", "randomly generated clustered points (Johnson)", "header"),
    "si1032":   ("table", "distance table only, no coordinates (Hofmeister)", "header"),
}
TYPES = [  # key, title
    ("cities", "Cities and towns"),
    ("drilling", "Circuit-board drilling"),
    ("chip", "Chip layout"),
    ("unstated", "Origin not stated"),
    ("synthetic", "Not real data: dsj1000"),
    ("table", "Distance table only: si1032"),
]
ALGS = ["greedy+2opt+oropt", "twoopt_or", "twoopt", "greedy", "nn", "spacefill",
        "box:hilbert:M1:t10", "box:hilbert:M5:t10", "box:knopp:M1:t10", "box:sierpinski:M1:t10",
        "box:sierpinski:M2:t12", "box:sierpinski:M5:t4", "curve:sierpinski", "curve:moore"]
TOP = 5
SURFACE, TEXT, MUTED, GRID, BAR, BAR1 = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e0", "#86b6ef", "#2a78d6"


def load():
    """runs[alg][map] = (tour length, user seconds) for the real-world maps, OK runs only."""
    runs = defaultdict(dict)
    with open(ROOT / "results" / "2d" / "runs.csv", newline="") as f:
        for r in csv.DictReader(f):
            if r["instance"] in MAPS and r["status"] == "OK" and r["algorithm"] in ALGS:
                runs[r["algorithm"]][r["instance"]] = (int(r["tour_length"]), float(r["user_s"]))
    return runs


def rank(runs, opt, maps):
    """Sorted rows (mean %, best-on count, mean seconds, alg) for algorithms with every map."""
    gap = {a: {m: 100 * (runs[a][m][0] / opt[m] - 1) for m in maps} for a in ALGS
           if all(m in runs[a] for m in maps)}
    best = {m: min(g[m] for g in gap.values()) for m in maps}
    rows = []
    for a, g in gap.items():
        wins = sum(g[m] <= best[m] + 1e-9 for m in maps)       # ties count for everyone tied
        secs = sum(runs[a][m][1] for m in maps) / len(maps)
        rows.append((sum(g.values()) / len(maps), wins, secs, a))
    return sorted(rows), gap


def chart(results, out):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10.5,
                         "axes.edgecolor": MUTED, "axes.labelcolor": MUTED,
                         "xtick.color": MUTED, "ytick.color": MUTED})
    fig, axes = plt.subplots(2, 3, figsize=(17, 8.6), facecolor=SURFACE)
    for ax, (key, title) in zip(axes.flat, TYPES):
        rows, maps = results[key]
        top = rows[:TOP]
        ax.set_facecolor(SURFACE)
        for y, (g, wins, secs, a) in enumerate(top):
            ax.barh(y, g, height=0.62, color=BAR1 if y == 0 else BAR)
            ax.text(g + top[-1][0] * 0.02, y, f"{g:.2f}%", va="center", color=TEXT, fontsize=10)
        ax.set_yticks(range(len(top)), [f"{i + 1}. {config_label(a)}" for i, (*_, a) in enumerate(top)],
                      color=TEXT)
        ax.set_ylim(len(top) - 0.4, -0.6)
        ax.set_xlim(0, top[-1][0] * 1.25)
        ax.grid(True, axis="x", color=GRID)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        n = f"{len(maps)} maps" if len(maps) > 1 else "1 map only (small sample)"
        extra = f"; only {len(rows)} solvers can read it" if len(rows) < len(ALGS) else ""
        ax.set_title(title, loc="left", color=TEXT, fontsize=11.5, fontweight="bold", pad=18)
        ax.text(0, 1.015, n + extra, transform=ax.transAxes, color=MUTED, fontsize=9)
        ax.set_xlabel("mean % above the optimal tour")
    fig.suptitle("Real-world maps (TSPLIB): the top 5 algorithms for each type of problem",
                 x=0.015, ha="left", color=TEXT, fontsize=15, fontweight="bold", y=0.985)
    fig.text(0.015, 0.905, "Rank = mean % above the proven optimal tour over the maps of that type "
             "(lower is better; 0% = perfect). All 14 algorithms were ranked; the full "
             "rankings are in results/2d/real_world.md.\n" + CODE_KEY + ".", color=MUTED, fontsize=9.5)
    fig.subplots_adjust(left=0.14, right=0.985, top=0.84, bottom=0.07, wspace=0.75, hspace=0.45)
    fig.savefig(out, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    print(f"Wrote {out.relative_to(ROOT)}")


def main():
    bounds = load_bounds(ROOT / "refs" / "bounds")
    opt = {m: bounds[m][1] for m in MAPS}
    runs = load()
    results = {}
    lines = ["# Real-world maps: algorithm rankings by type of problem", "",
             "Generated by `bench/real_world.py` from `results/2d/runs.csv`. Rank = mean % above",
             "the proven optimal tour (`refs/bounds`) over the maps of that type; lower is",
             "better. \"Best on\" = how many of the maps it had the shortest tour on (ties",
             "count for each). " + CODE_KEY + ".", ""]
    for key, title in TYPES:
        maps = sorted((m for m, v in MAPS.items() if v[0] == key), key=lambda m: int("".join(filter(str.isdigit, m))))
        rows, gap = rank(runs, opt, maps)
        results[key] = (rows, maps)
        lines += [f"## {title}", "",
                  "Maps: " + ", ".join(f"{m} ({MAPS[m][1]})" for m in maps), "",
                  "| Rank | Algorithm | Mean % above optimal | Best on | Mean seconds |",
                  "|---|---|---|---|---|"]
        lines += [f"| {i + 1} | {config_label(a)} | {g:.2f} | {w} of {len(maps)} | {s:.2f} |"
                  for i, (g, w, s, a) in enumerate(rows)]
        missing = [config_label(a) for a in ALGS if a not in gap]
        if missing:
            lines += ["", "Not ranked (can't read this map): " + ", ".join(missing) + "."]
        if len(maps) > 1:
            lines += ["", "Per map, % above optimal:", "",
                      "| Algorithm | " + " | ".join(maps) + " |", "|---|" + "---|" * len(maps)]
            lines += [f"| {config_label(a)} | " + " | ".join(f"{gap[a][m]:.1f}" for m in maps) + " |"
                      for _, _, _, a in rows]
        lines.append("")
        print(f"{title}:")
        for i, (g, w, s, a) in enumerate(rows[:TOP]):
            print(f"  {i + 1}. {config_label(a):24} {g:6.2f}%  best on {w}/{len(maps)}  {s:.2f}s")

    out = ROOT / "results" / "2d" / "real_world.md"
    out.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {out.relative_to(ROOT)}")
    chart(results, ROOT / "results" / "2d" / "real_world.png")


if __name__ == "__main__":
    main()
