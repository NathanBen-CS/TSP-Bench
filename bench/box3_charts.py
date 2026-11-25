"""Charts for the 3D study (results/3d/runs.csv), in the style of the 2D ones.

    py bench/box3_charts.py        # writes results/3d/box3_*.png

  box3_curve_merge.png        curve x merge grid at t = 10, 100k cities and the stars:
                              after full cleanup (TO Rn) and before any cleanup
  box3_overview.png           baselines and the best box settings: % above our own
                              Held-Karp bound (100k) and above the best tour (stars)
  box3_uniform_vs_clustered.png  how far behind greedy TO Rn each merge ends, by map type
  box3_t_sweep.png            box size t, 10k cities (hilbert3)
  box3_t_sweep_sierpinski3.png  the same for sierpinski3, whose boxes come close to t
  box3_k_sweep.png            neighbour-list size K, quality
  box3_k_time.png             neighbour-list size K, time
  box3_scale.png              1k to 1M cities
  box3_1M.png                 every setting run at 1 million cities and on gaia1
  box3_kd.png                 kd-median boxes against curve cells

Main measure everywhere but the overview: % longer than greedy TO Rn (K = 16) on the same map; below 0
would beat it.  Mean over the 3 seeds (stars: one map each).
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

from box3_report import BASE, CURVES, MERGES, OUT, SIZES, Data, mean, sd
from common import CODE_KEY, config_label

SURFACE, TEXT, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e0"
KINDS = [("E3d", "uniform", "#2a78d6"), ("C3d", "clustered", "#eb6834"),
         ("CT3d", "tight clustered", "#1baf7a"), ("stars", "stars", "#4a3aa7")]
BIG_STARS = ["star10k", "kj37859", "hyg109399", "star250k"]
MERGE_COL = {"M1": "#2a78d6", "M2": "#1baf7a", "M3": "#eda100", "M4": "#eb6834", "M5": "#4a3aa7"}


def style():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.edgecolor": MUTED,
                         "axes.labelcolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED,
                         "figure.facecolor": SURFACE})


def clean(ax, grid_axis="x"):
    ax.set_facecolor(SURFACE)
    ax.grid(True, axis=grid_axis, color=GRID, linewidth=1)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def head(fig, title, sub, x=0.05, y=0.985, gap=0.04):
    fig.suptitle(title, x=x, ha="left", color=TEXT, fontsize=14, fontweight="bold", y=y)
    fig.text(x, y - gap, sub, color=MUTED, fontsize=9.5, va="top")


def maps_of(d, prefix, size):
    if prefix == "stars":
        return [m for m in BIG_STARS if (BASE, m) in d.runs]
    return [m for m in d.group(prefix, size) if (BASE, m) in d.runs]


def vs(d, name, maps, lvl="full"):
    v = d.vs_base(name, maps, lvl)
    return mean(v) if v and len(v) == len(maps) else None


# ---------------------------------------------------------------- charts

def curve_merge(d):
    curves = [c for c in CURVES if vs(d, f"box:{c}:M1:t10", maps_of(d, "E3d", "100k")) is not None]
    fig, axes = plt.subplots(2, 4, figsize=(17, 8.6))
    for row, lvl in enumerate(["full", "none"]):
        for ax, (p, kind, _) in zip(axes[row], KINDS):
            ms = maps_of(d, p, "100k")
            grid = [[vs(d, f"box:{c}:{m}:t10", ms, lvl) for m in MERGES] for c in curves]
            lim = 8 if lvl == "full" else 60
            ax.imshow([[0 if v is None else max(-lim, min(lim, v)) for v in r] for r in grid], cmap="RdBu_r",
                      norm=TwoSlopeNorm(vcenter=0, vmin=-lim, vmax=lim), aspect="auto")
            for i, r in enumerate(grid):
                for j, v in enumerate(r):
                    txt = "–" if v is None else f"{v:+.1f}" if abs(v) < 100 else f"{v:+,.0f}"
                    ax.text(j, i, txt, ha="center", va="center", fontsize=9.5,
                            color="white" if v is not None and abs(v) > lim * 0.65 else TEXT)
            ax.set_xticks(range(len(MERGES)), MERGES)
            ax.set_yticks(range(len(curves)), curves if p == "E3d" else [""] * len(curves))
            for lbl in ax.get_yticklabels():
                lbl.set_color(TEXT)
            ax.tick_params(length=0)
            for s in ax.spines.values():
                s.set_visible(False)
            size = "100k" if p != "stars" else "mean of 4 maps"
            ax.set_title(f"{kind}, {size}: " + ("after cleanup (BTO Rn)" if lvl == "full" else "boxes only (B)"),
                         loc="left", color=TEXT, fontsize=10.5, fontweight="bold")
    head(fig, "3D, curve × merge at t = 10: % longer than greedy TO Rn",
         "Every cell is a box setting against greedy TO Rn on the same maps: above 0 (red) = longer, below 0 (blue) "
         "= would beat it. Mean of 3 seeds; stars: star10k, kj37859, hyg109399, star250k.\n"
         "Colours are capped at ±8 (top) and ±60 (bottom); numbers are exact. " + CODE_KEY + ".", x=0.06)
    fig.subplots_adjust(left=0.075, right=0.99, top=0.85, bottom=0.04, hspace=0.25, wspace=0.06)
    fig.savefig(OUT / "box3_curve_merge.png", dpi=160, facecolor=SURFACE)


def best_settings(d, k=4):
    """the k box settings with the best mean over uniform, clustered, tight (100k) and the stars"""
    score = []
    for c in CURVES:
        for m in MERGES:
            n = f"box:{c}:{m}:t10"
            v = [vs(d, n, maps_of(d, p, "100k")) for p, _, _ in KINDS]
            if all(x is not None for x in v):
                score.append((mean(v), n))
    return [n for _, n in sorted(score)[:k]]


def overview(d):
    best = best_settings(d)
    rows = ["curve:hilbert3", "greedy"]
    rows += [n for n in ["box:hilbert3:M1:t10", "box:sierpinski3:M5:t10"] + best[::-1] if n not in rows]
    rows += ["twoopt", "greedy+2opt", "twoopt_or", BASE]
    fig, axes = plt.subplots(1, 2, figsize=(15, 7), sharey=True)
    panels = [("100k cities: % above OUR OWN Held-Karp bound", [("E3d", "100k"), ("C3d", "100k")], "hk"),
              ("stars: % above the published best tour (mean of 4 maps)",
               [("stars", "")], "best")]
    for ax, (title, sets, kind) in zip(axes, panels):
        clean(ax)
        top = 0
        for i, r in enumerate(rows):
            for j, (p, s) in enumerate(sets):
                ms = maps_of(d, p, s)
                g = mean(d.vs_hk(r, ms)) if kind == "hk" else mean([d.vs_best(r, m) for m in ms if d.vs_best(r, m) is not None])
                if g is None:
                    continue
                off = 0 if len(sets) == 1 else (-0.2 if j == 0 else 0.2)
                ax.barh(i + off, min(g, 10), height=0.36 if len(sets) > 1 else 0.6,
                        color=[c for q, _, c in KINDS if q == p][0])
                ax.text(min(g, 10) + 0.12, i + off, f"{g:.1f}%" + (" (cut)" if g > 10 else ""),
                        va="center", fontsize=8.8, color=TEXT)
        ax.set_yticks(range(len(rows)), [config_label(r) for r in rows])
        ax.invert_yaxis()
        ax.set_xlim(0, 12)
        ax.set_xlabel("% above (lower is better)")
        ax.set_title(title, loc="left", color=TEXT, fontsize=10.5, fontweight="bold")
    for lbl in axes[0].get_yticklabels():
        lbl.set_color(TEXT)
    fig.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=c, label=k) for p, k, c in KINDS if p != "CT3d"],
               loc="lower left", bbox_to_anchor=(0.19, 0.0), ncol=3, frameon=False, labelcolor=TEXT)
    head(fig, "3D: every solver and the best box settings", "Left: our own 1-tree bound, not a published one "
         "(up to about 0.3% too low on clustered maps). Right: Waterloo tours of star10k, kj37859, hyg109399, star250k.\n" + CODE_KEY + ".", x=0.19)
    fig.subplots_adjust(left=0.19, right=0.985, top=0.84, bottom=0.13, wspace=0.06)
    fig.savefig(OUT / "box3_overview.png", dpi=160, facecolor=SURFACE)


def uniform_vs_clustered(d):
    names = [f"box:{c}:{m}:t10" for c in ("hilbert3", "peano3", "moore3", "sierpinski3", "morton3") for m in MERGES]
    names = [n for n in names if vs(d, n, maps_of(d, "E3d", "100k")) is not None]
    fig, ax = plt.subplots(figsize=(12.5, 0.36 * len(names) + 2.6))
    clean(ax)
    y, ticks, labels, prev = 0, [], [], None
    for n in names:
        if prev and n.split(":")[1] != prev:
            y += 0.6
        prev = n.split(":")[1]
        xs = []
        for p, _, col in KINDS:
            v = vs(d, n, maps_of(d, p, "100k"))
            if v is None:
                continue
            xs.append(v)
            ax.scatter([v], [y], s=62, color=col, zorder=3, edgecolor=SURFACE, linewidth=1.5)
        ax.plot([min(xs), max(xs)], [y, y], color=GRID, lw=5, zorder=1, solid_capstyle="round")
        ax.annotate(f"{min(xs):+.1f} to {max(xs):+.1f}", xy=(max(xs), y), xytext=(9, 0), textcoords="offset points",
                    va="center", fontsize=8.8, color=TEXT)
        ticks.append(y)
        labels.append(config_label(n))
        y += 1
    ax.axvline(0, color=TEXT, lw=1.2)
    ax.text(0.04, -1.1, "greedy TO Rn", fontsize=9, color=TEXT)
    ax.set_yticks(ticks, labels)
    for lbl in ax.get_yticklabels():
        lbl.set_color(TEXT)
    ax.set_ylim(y, -1.8)
    ax.set_xlim(-0.3, ax.get_xlim()[1] + 0.5)
    ax.set_xlabel("% longer than greedy TO Rn on the same map (lower is better)")
    ax.legend(handles=[plt.Line2D([], [], marker="o", ls="", ms=8, color=c, label=k + (", 100k" if p != "stars" else " (4 maps)"))
                       for p, k, c in KINDS], loc="upper right", frameon=False, labelcolor=TEXT)
    head(fig, "3D, by map type: how far behind greedy TO Rn each box setting ends", "100k cities, mean of 3 seeds; stars: "
         "star10k, kj37859, hyg109399, star250k. snake3 is left out (about +20 on tight clustered).\n" + CODE_KEY + ".",
         x=0.03, gap=0.025)
    fig.subplots_adjust(left=0.23, right=0.96, top=0.91, bottom=0.07)
    fig.savefig(OUT / "box3_uniform_vs_clustered.png", dpi=150, facecolor=SURFACE)


def t_sweep(d, curve="hilbert3"):
    ts = list(range(4, 13))
    fig, axes = plt.subplots(1, 4, figsize=(17, 5.4), sharey=True)
    for ax, (p, kind, _) in zip(axes, KINDS):
        clean(ax, "y")
        ms = ["star10k"] if p == "stars" else maps_of(d, p, "10k")
        ends = []
        for m in MERGES:
            ys = [vs(d, f"box:{curve}:{m}:t{t}", ms) for t in ts]
            if any(v is None for v in ys):
                continue
            ax.plot(ts, ys, color=MERGE_COL[m], lw=2, marker="o", ms=5)
            ends.append([ys[-1], m])
        ends.sort()
        pos = [e[0] for e in ends]
        for i in range(1, len(pos)):
            pos[i] = max(pos[i], pos[i - 1] + 0.18)
        for (v, m), q in zip(ends, pos):
            ax.text(12.3, q, m, color=MERGE_COL[m], fontsize=9.5, fontweight="bold", va="center")
        ax.axhline(0, color=TEXT, lw=1.2)
        ax.set_xticks(ts)
        ax.set_xlim(3.6, 13.4)
        ax.set_xlabel("box size t")
        ax.set_title("star10k (1 map)" if p == "stars" else f"{kind}, 10k (3 seeds)", loc="left", color=TEXT,
                     fontsize=11, fontweight="bold")
    axes[0].set_ylabel("% longer than greedy TO Rn")
    axes[0].set_ylim(bottom=-0.25)
    axes[0].text(4, 0.06, "greedy TO Rn", fontsize=9, color=TEXT)
    head(fig, f"3D, box size t, {curve} BTO Rn", "Every line is one merge, after full cleanup. The zero line is greedy "
         "TO Rn. " + CODE_KEY + ".", x=0.05, gap=0.055)
    fig.subplots_adjust(left=0.05, right=0.99, top=0.82, bottom=0.11, wspace=0.08)
    fig.savefig(OUT / ("box3_t_sweep.png" if curve == "hilbert3" else f"box3_t_sweep_{curve}.png"),
                dpi=160, facecolor=SURFACE)


def t_sweep_sierpinski3(d):
    t_sweep(d, "sierpinski3")


K_ROWS = [("greedy TO Rn", None, TEXT), ("hilbert3 M1", "box:hilbert3:M1:t10", MERGE_COL["M1"]),
          ("hilbert3 M4", "box:hilbert3:M4:t10", MERGE_COL["M4"]), ("hilbert3 M5", "box:hilbert3:M5:t10", MERGE_COL["M5"]),
          ("sierpinski3 M5", "box:sierpinski3:M5:t10", MERGE_COL["M5"])]   # drawn dashed


def k_name(base, k):
    if base is None:
        return BASE if k == 16 else f"greedy TO Rn:k{k}"
    return base if k == 16 else f"{base}:k{k}"


def k_sweep(d):
    ks = [8, 16, 24, 32]
    fig, axes = plt.subplots(1, 4, figsize=(17, 5.4), sharey=True)
    for ax, (p, kind, _) in zip(axes, KINDS):
        clean(ax, "y")
        ms = ["kj37859", "hyg109399"] if p == "stars" else maps_of(d, p, "100k")
        ends = []
        for label, base, col in K_ROWS:
            ys = [vs(d, k_name(base, k), ms) for k in ks]
            if any(v is None for v in ys):
                continue
            ax.plot(ks, ys, color=col, lw=2.4 if base is None else 2, marker="o", ms=6,
                    ls="--" if "sierpinski3" in label else "-")
            ends.append([ys[-1], f"{label.replace('hilbert3 ', '').replace('sierpinski3', 'sierp3').replace('greedy TO Rn', 'greedy')} {ys[-1] if abs(ys[-1]) >= 0.05 else 0:+.1f}", col])
        ends.sort(key=lambda e: e[0])
        pos = [e[0] for e in ends]
        for i in range(1, len(pos)):
            pos[i] = max(pos[i], pos[i - 1] + 0.22)
        for (v, txt, col), q in zip(ends, pos):
            ax.text(33, q, txt, color=col, fontsize=9, fontweight="bold", va="center")
        ax.axhline(0, color=GRID, lw=1.2, zorder=0.5)
        ax.set_xticks(ks, [f"K={k}" for k in ks])
        ax.set_xlim(6.5, 44)
        ax.set_title("stars (kj37859, hyg109399)" if p == "stars" else f"{kind}, 100k (3 seeds)", loc="left",
                     color=TEXT, fontsize=11, fontweight="bold")
    axes[0].set_ylabel("% longer than greedy TO Rn with K = 16")
    head(fig, "3D, neighbour-list size K: quality", "K sets the lists of the local search and of merges M3–M5, and greedy's "
         "round-1 candidates (min(10, K)). BTO Rn at t = 10 (hilbert3 solid, sierpinski3 dashed); black = greedy TO Rn at the same K.", x=0.05, gap=0.055)
    fig.subplots_adjust(left=0.05, right=0.99, top=0.82, bottom=0.08, wspace=0.08)
    fig.savefig(OUT / "box3_k_sweep.png", dpi=160, facecolor=SURFACE)


def k_time(d):
    ks = [8, 16, 24, 32]
    shades = {8: "#b8d4f5", 16: "#5b9be0", 24: "#2a5f9e", 32: "#0d2f5c"}
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), sharey=True, sharex=True)
    for ax, (p, kind, _) in zip(axes, KINDS[:3]):
        clean(ax)
        ms = maps_of(d, p, "100k")
        for i, (label, base, _) in enumerate(K_ROWS):
            for j, k in enumerate(ks):
                u = d.user(k_name(base, k), ms)
                if u:
                    ax.barh(i + (j - 1.5) * 0.21, mean(u), height=0.18, color=shades[k])
                    ax.text(mean(u) + 0.1, i + (j - 1.5) * 0.21, f"{mean(u):.1f} s", va="center", fontsize=8, color=TEXT)
        ax.set_yticks(range(len(K_ROWS)), [r[0] for r in K_ROWS])
        ax.set_title(f"{kind}, 100k", loc="left", color=TEXT, fontsize=11, fontweight="bold")
        ax.set_xlabel("user CPU seconds")
    axes[0].invert_yaxis()
    for lbl in axes[0].get_yticklabels():
        lbl.set_color(TEXT)
    fig.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=shades[k], label=f"K = {k}") for k in ks], loc="upper right",
               bbox_to_anchor=(0.985, 0.97), ncol=4, frameon=False, labelcolor=TEXT)
    head(fig, "3D, neighbour-list size K: time", "Whole program, TO Rn cleanup, mean of 3 seeds, 100k cities.", x=0.05, gap=0.06)
    fig.subplots_adjust(left=0.09, right=0.98, top=0.8, bottom=0.12, wspace=0.06)
    fig.savefig(OUT / "box3_k_time.png", dpi=160, facecolor=SURFACE)


def scale(d):
    rows = [("twoopt_or", "greedy TO R1"), ("box:hilbert3:M1:t10", None), ("box:hilbert3:M4:t10", None),
            ("box:hilbert3:M5:t10", None), ("box:moore3:M5:t10", None), ("box:peano3:M5:t10", None),
            ("box:sierpinski3:M5:t10", None)]
    rows = [(n, l or config_label(n)) for n, l in rows if vs(d, n, maps_of(d, "E3d", "100k")) is not None]
    fig, axes = plt.subplots(1, len(rows), figsize=(3.0 * len(rows), 5), sharey=True)
    for ax, (n, label) in zip(axes, rows):
        clean(ax, "y")
        ends = []
        for p, kind, col in KINDS[:3]:
            ys = [vs(d, n, maps_of(d, p, s)) if maps_of(d, p, s) else None for s in SIZES]
            xs = [i for i, v in enumerate(ys) if v is not None]
            ys = [v for v in ys if v is not None]
            ax.plot(xs, ys, color=col, lw=2, marker="o", ms=6)
            if ys:
                ends.append([ys[-1], xs[-1], col])
        ends.sort()
        pos = [e[0] for e in ends]
        for i in range(1, len(pos)):
            pos[i] = max(pos[i], pos[i - 1] + 0.2)
        for (v, x, col), q in zip(ends, pos):
            ax.text(x + 0.12, q, f"{v:+.1f}", color=col, fontsize=8.8, fontweight="bold", va="center")
        ax.axhline(0, color=TEXT, lw=1.2)
        ax.set_xticks(range(4), SIZES)
        ax.set_xlim(-0.3, 3.8)
        ax.set_title(label, loc="left", color=TEXT, fontsize=10, fontweight="bold")
    axes[0].set_ylabel("% longer than greedy TO Rn")
    fig.legend(handles=[plt.Line2D([], [], marker="o", lw=2, ms=6, color=c, label=k) for _, k, c in KINDS[:3]],
               loc="lower left", bbox_to_anchor=(0.05, 0.0), ncol=3, frameon=False, labelcolor=TEXT)
    head(fig, "3D, 1k to 1M cities: does any setting catch up with greedy TO Rn?", "The black zero line is greedy TO Rn. "
         "Mean of 3 seeds (1M: one map per type).\n" + CODE_KEY + ".", x=0.05, gap=0.06)
    fig.subplots_adjust(left=0.05, right=0.99, top=0.77, bottom=0.16, wspace=0.08)
    fig.savefig(OUT / "box3_scale.png", dpi=160, facecolor=SURFACE)


def million(d):
    cols = [("E3d1M.0", "uniform 1M"), ("C3d1M.0", "clustered 1M"), ("CT3d1M.0", "tight clustered 1M"), ("gaia1", "gaia1 (2.08M stars)")]
    names = sorted({a for a, m in d.runs if m in ("E3d1M.0", "gaia1") and not a.startswith("curve:")
                    and a not in ("nn", "greedy") and ":k" not in a and not a.endswith("R1")})
    key = lambda n: mean([v for m, _ in cols for v in d.vs_base(n, [m])]) or 0
    names = sorted(names, key=key)
    fig, axes = plt.subplots(1, 4, figsize=(18, 0.33 * len(names) + 2.4), sharey=True)
    colours = [KINDS[0][2], KINDS[1][2], KINDS[2][2], KINDS[3][2]]
    for ax, (m, title), col in zip(axes, cols, colours):
        clean(ax)
        for i, n in enumerate(names):
            v = d.vs_base(n, [m])
            if not v:
                continue
            ax.barh(i, min(v[0], 12), height=0.62, color=col)
            ax.text(min(max(v[0], 0), 12) + 0.1, i, f"{v[0]:+.1f}%" + "".join(f"   {u:.0f} s" for u in d.user(n, [m])), va="center", fontsize=8.5, color=TEXT)
        ax.axvline(0, color=TEXT, lw=1.2)
        ax.set_xlim(-0.2, 13)
        ax.set_title(title, loc="left", color=TEXT, fontsize=11, fontweight="bold")
        ax.set_xlabel("% longer than greedy TO Rn")
    axes[0].set_yticks(range(len(names)), [config_label(n) for n in names])
    axes[0].invert_yaxis()
    for lbl in axes[0].get_yticklabels():
        lbl.set_color(TEXT)
    head(fig, "3D at 1 million cities and more", "One map per column. Each bar: % longer than greedy TO Rn, then the user "
         "CPU seconds of the whole program.\n" + CODE_KEY + ".", x=0.05, gap=0.05)
    fig.subplots_adjust(left=0.15, right=0.99, top=0.85, bottom=0.1, wspace=0.06)
    fig.savefig(OUT / "box3_1M.png", dpi=150, facecolor=SURFACE)


def kd(d):
    ts = [4, 6, 8, 10, 12]
    fig, axes = plt.subplots(2, 4, figsize=(17, 8.2), sharex=True, sharey="row")
    for row, lvl in enumerate(["none", "full"]):
        for ax, (p, kind, _) in zip(axes[row], KINDS):
            clean(ax, "y")
            ms = ["star10k"] if p == "stars" else maps_of(d, p, "10k")
            for m in ("M1", "M4", "M5"):
                a = [vs(d, f"kd:hilbert3:{m}:t{t}", ms, lvl) for t in ts]
                b = [vs(d, f"box:hilbert3:{m}:t{t}", ms, lvl) for t in ts]
                if any(v is None for v in a + b):
                    continue
                ax.plot(ts, a, color=MERGE_COL[m], lw=2, marker="o", ms=5)
                ax.plot(ts, b, color=MERGE_COL[m], lw=1.6, ls=(0, (4, 3)), marker="o", ms=4, mfc=SURFACE)
            ax.axhline(0, color=TEXT, lw=1.2)
            ax.set_xticks(ts, [f"t={t}" for t in ts])
            ax.set_title(("star10k" if p == "stars" else f"{kind}, 10k") +
                         (": boxes only (B)" if lvl == "none" else ": after cleanup (BTO Rn)"),
                         loc="left", color=TEXT, fontsize=10.5, fontweight="bold")
        axes[row][0].set_ylabel("% longer than greedy TO Rn")
    handles = [plt.Line2D([], [], color=MERGE_COL[m], lw=2, marker="o", label=m) for m in ("M1", "M4", "M5")]
    handles += [plt.Line2D([], [], color=MUTED, lw=2, label="kd boxes (hold exactly t)"),
                plt.Line2D([], [], color=MUTED, lw=1.6, ls=(0, (4, 3)), label="curve cells (hold at most t)")]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.05, 0.0), ncol=5, frameon=False, labelcolor=TEXT)
    head(fig, "3D, kd-median boxes against curve-cell boxes (hilbert3 order)", "The black zero line is greedy TO Rn. " + CODE_KEY + ".",
         x=0.05, gap=0.04)
    fig.subplots_adjust(left=0.05, right=0.99, top=0.87, bottom=0.1, wspace=0.08, hspace=0.2)
    fig.savefig(OUT / "box3_kd.png", dpi=150, facecolor=SURFACE)


def main():
    style()
    d = Data()
    for chart in (curve_merge, overview, uniform_vs_clustered, t_sweep, t_sweep_sierpinski3, k_sweep, k_time, scale, million, kd):
        chart(d)
        plt.close("all")
        print("  ", chart.__name__, flush=True)


if __name__ == "__main__":
    main()
