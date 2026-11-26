"""3D: best box setting (after cleanup) against greedy plus cleanup, by size. Reads report/data/table_3d_best_setting.md."""
from pathlib import Path
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
rows = {}
for line in open(ROOT / 'report/data/table_3d_best_setting.md'):
    p = [x.strip() for x in line.strip().strip('|').split('|')]
    if len(p) == 3 and re.fullmatch(r'[+-]\d+\.\d+', p[2]):
        rows[p[0]] = float(p[2])

size = {'1k': 1e3, '10k': 1e4, '100k': 1e5, '1M': 1e6}
series = []
for key, label, color, marker in [('uniform', 'uniform', '#257de0', 'o'), ('clustered', 'clustered', '#f47230', 'o'),
                                  ('tight clustered', 'tight clustered', '#00b894', 'o')]:
    pts = [(size[s], rows[f'{key} {s}']) for s in size if f'{key} {s}' in rows]
    series.append((label, color, marker, pts))
stars = {'star1k': 1e3, 'star10k': 1e4, 'kj37859': 37859, 'hyg109399': 109399, 'hyg119614': 119614,
         'star250k': 250000, 'gaia1': 2079471}
series.append(('star catalogues', '#555555', 's', [(stars[k], rows[k]) for k in stars if k in rows]))

fig, ax = plt.subplots(figsize=(7.2, 3.9), facecolor='#fcfcfb')
ax.set_facecolor('#fcfcfb')
ax.axhline(0, color='#222', lw=1.2)
ax.text(1.1e3, 0.08, 'greedy plus cleanup (0)', fontsize=9, color='#222', va='bottom')
offs = {'uniform': 0.0, 'clustered': -0.18, 'tight clustered': -0.38, 'star catalogues': 0.28}
for label, color, marker, pts in series:
    xs, ys = zip(*pts)
    ax.plot(xs, ys, color=color, lw=2.2, marker=marker, ms=6.5, mec='#fcfcfb', mew=1.2, ls='-' if label != 'star catalogues' else '--')
    ax.text(xs[-1] * 1.12, ys[-1] + offs[label], f'{label}  {ys[-1]:+.1f}', fontsize=9.5, color=color if label != 'star catalogues' else '#444',
            fontweight='bold', va='center')
ax.set_xscale('log')
ax.set_xlim(8e2, 2.2e7)
ax.set_ylim(-0.3, 4)
ax.set_xticks([1e3, 1e4, 1e5, 1e6])
ax.set_xticklabels(['1k', '10k', '100k', '1M'])
ax.set_xlabel('cities or stars')
ax.set_ylabel('% longer than greedy plus cleanup')
for s in ('top', 'right'):
    ax.spines[s].set_visible(False)
ax.grid(axis='y', color='#e4e4e0', lw=0.8)
ax.set_axisbelow(True)
fig.tight_layout()
out = ROOT / 'report/figures/fig_3d.png'
fig.savefig(out, dpi=200, facecolor='#fcfcfb')
print('wrote', out)
