"""Hero figure of the project report: raw tour -> finished tour, 2D 100k maps.

    py report/scripts/fig_hero.py        # writes report/figures/fig_hero.png

Numbers come from the same run files as report/data/matrix_at_a_glance_sd.md
(one hollow dot = the raw tour, one filled dot = after 2-opt + Or-opt to a local
optimum; % above the Held-Karp bound, mean over the maps of the group, log axis).
Boxes are the mean of the 25 settings (5 curves x 5 merges), hull rings of the 6
merges, greedy the single baseline. The "+-" is the sample standard deviation of
the finished tour over the maps (n = 2 here, so read it as a rough spread).
"""
import sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
sys.path.insert(0, str(Path(__file__).resolve().parent))
import tables as T

SURF, INK, SEC, GRID, AXIS = '#fcfcfb', '#191919', '#666666', '#e7e7e7', '#c2c2c2'
FAM = [('greedy + cleanup', [T.G], '#257de0'),
       ('boxes (25 settings)', [f'box:{c}:{m}:t10' for c in T.CURVES[:5] for m in T.MER], '#f47230'),
       ('hull rings (6 merges)', T.HULL, '#00b894')]
LV = ['raw', 'TO Rn']


def stats(kind):
    ms = [m for m in T.grp(kind, '100k') if m in T.HK]
    out = []
    for label, names, col in FAM:
        names = [n for n in names if all(T.level_len(n, m, l) for m in ms for l in LV)]
        per = {l: [T.mean([(T.level_len(n, m, l) / T.HK[m] - 1) * 100 for n in names]) for m in ms] for l in LV}
        out.append((label, col, T.mean(per['raw']), T.mean(per['TO Rn']), T.sd(per['TO Rn'])))
    return ms, out


def fmt(x):
    return f'{x:,.0f}' if x >= 1000 else f'{x:.1f}'


fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.0), sharex=True, sharey=True)
fig.patch.set_facecolor(SURF)
for ax, (kind, title) in zip(axes, [('E', 'uniform 100k'), ('C', 'clustered 100k')]):
    ms, rows = stats(kind)
    ax.set_facecolor(SURF)
    for i, (label, col, raw, fin, s) in enumerate(rows):
        y = len(rows) - 1 - i
        ax.plot([fin, raw], [y, y], color=col, lw=2, solid_capstyle='round', zorder=2)
        ax.scatter([raw], [y], s=70, facecolor=SURF, edgecolor=col, linewidth=2, zorder=3)
        ax.scatter([fin], [y], s=70, facecolor=col, edgecolor=SURF, linewidth=1.5, zorder=4)
        spread = '' if s is None else ' ± ' + ('<0.1' if s < 0.05 else f'{s:.1f}')
        ax.annotate(fmt(raw), (raw, y), xytext=(0, 9), textcoords='offset points', ha='center', fontsize=8.5, color=SEC)
        ax.annotate(fmt(fin) + spread, (fin, y), xytext=(0, -17), textcoords='offset points', ha='center',
                    fontsize=9, color=INK, fontweight='bold')
        if ax is axes[0]:
            ax.text(-0.02, y, label, transform=ax.get_yaxis_transform(), ha='right', va='center', fontsize=9.5, color=INK)
    ax.set_xscale('log')
    ax.set_xlim(2.5, 6000)
    ax.set_ylim(-0.6, len(rows) - 0.35)
    ax.set_yticks([])
    ax.set_xticks([3, 10, 30, 100, 300, 1000, 3000])
    ax.set_xticklabels(['3', '10', '30', '100', '300', '1,000', '3,000'])
    ax.minorticks_off()
    ax.tick_params(axis='x', colors=SEC, labelsize=8.5, length=0)
    ax.grid(axis='x', color=GRID, lw=0.8, zorder=0)
    for sp in ('top', 'right', 'left'):
        ax.spines[sp].set_visible(False)
    ax.spines['bottom'].set_color(AXIS)
    ax.set_title(f'{title}  (n = {len(ms)} maps)', loc='left', fontsize=10.5, color=INK, pad=8)
    ax.set_xlabel('% above the Held-Karp bound (log scale)', fontsize=8.5, color=SEC, labelpad=6)
handles = [Line2D([], [], marker='o', ls='', markersize=8, markerfacecolor=SURF, markeredgecolor=SEC, markeredgewidth=2,
                  label='raw tour (before cleanup)'),
           Line2D([], [], marker='o', ls='', markersize=8, markerfacecolor=SEC, markeredgecolor=SURF,
                  label='after 2-opt + Or-opt  (value ± sample SD over maps)')]
fig.legend(handles=handles, loc='lower center', ncol=2, frameon=False, fontsize=8.8, labelcolor=SEC,
           bbox_to_anchor=(0.55, -0.015))
fig.subplots_adjust(left=0.17, right=0.985, top=0.9, bottom=0.25, wspace=0.06)
out = Path(__file__).resolve().parents[1] / 'figures' / 'fig_hero.png'
fig.savefig(out, dpi=200, facecolor=SURF)
print('wrote', out)
