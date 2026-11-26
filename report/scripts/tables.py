"""Rebuild the matrix tables of the project report from the run files.

    py report/scripts/tables.py            # writes report/data/*.md

Reads results/2d/runs.csv, results/2d/box2/runs.csv, results/2d/hull/runs.csv,
results/3d/runs.csv and refs/bounds. Nothing in the repo is changed except
the files it writes in report/data/. Every table is a mean over the maps of a
group; a cell is left as a dash when a run does not exist.
"""
import collections, csv, os, re, statistics as st
from pathlib import Path

ROOT = Path(os.environ.get('TSP_ROOT', Path(__file__).resolve().parents[2]))
OUT = ROOT / 'report' / 'data'
OUT.mkdir(parents=True, exist_ok=True)
mean = lambda v: st.mean(v) if v else None


def load(path):
    """(algorithm, instance) -> first OK row: length, extras dict, user time, peak MB, n."""
    d = {}
    for r in csv.DictReader(open(path)):
        k = (r['algorithm'], r['instance'])
        if k in d or r['status'] != 'OK':
            continue
        ex = dict(kv.split('=', 1) for kv in r['extras'].split(';') if '=' in kv)
        d[k] = dict(L=float(r['tour_length']), ex=ex, user=float(r['user_s']),
                    mb=float(r['peak_mb']), n=int(r['n']))
    return d


HK = {}
for line in open(ROOT / 'refs' / 'bounds'):
    p = line.split()
    if len(p) >= 2:
        HK[p[0]] = float(p[1])
HR = load(ROOT / 'results/2d/hull/runs.csv')      # all curves x merges, hull, greedy, R1 variants
BR = load(ROOT / 'results/2d/runs.csv')           # DIMACS bench: nn, greedy, spacefill, twoopt, ...
B2 = load(ROOT / 'results/2d/box2/runs.csv')      # bare curve orders for knopp, peano, morton, strip
R3 = load(ROOT / 'results/3d/runs.csv')

CURVES = ['hilbert', 'moore', 'sierpinski', 'knopp', 'peano', 'morton', 'strip']
MER = ['M1', 'M2', 'M3', 'M4', 'M5']
HULL = [f'hull:{m}' for m in ['H0'] + MER]
G = 'greedy+2opt+oropt'
maps2 = sorted({m for _, m in HR})
COLS = [('E', '1k'), ('E', '10k'), ('E', '100k'), ('E', '1M'), ('C', '1k'), ('C', '10k'), ('C', '100k')]
HEAD = ('| Method | uniform 1k | uniform 10k | uniform 100k | uniform 1M | clustered 1k | '
        'clustered 10k | clustered 100k |\n|---|---|---|---|---|---|---|---|')


def grp(kind, size):
    return [m for m in maps2 if re.fullmatch(rf'{kind}{size}\.\d', m)]


def length(D, name, m, level='full'):
    r = D.get((name, m))
    if not r:
        return None
    if level == 'full':
        return r['L']
    key = {'raw': 'start_length', 'T Rn': 'after_2opt_length'}[level]
    return float(r['ex'][key]) if key in r['ex'] else None


def hk_gap(D, name, ms, level='full'):
    return [(length(D, name, m, level) / HK[m] - 1) * 100 for m in ms if m in HK and length(D, name, m, level)]


def vs_greedy(D, name, ms, level='full'):
    return [(length(D, name, m, level) / HR[(G, m)]['L'] - 1) * 100 for m in ms
            if (G, m) in HR and length(D, name, m, level)]


def f1(x, big=1000):
    return '-' if x is None else (f'{x:.1f}' if x < big else f'{x:,.0f}')


def best_cell(D, names, ms, level, score):
    """merge/setting with the smallest score over ms; returns (value-as-text, label) or None"""
    cand = [n for n in names if all((n, m) in D for m in ms)]
    sc = {n: mean(score(D, n, ms, level)) for n in cand}
    sc = {n: v for n, v in sc.items() if v is not None}
    if not sc:
        return None
    b = min(sc, key=sc.get)
    return b, sc[b]


def save(name, text):
    (OUT / name).write_text(text + '\n', encoding='utf-8')
    print('wrote', name)


# ---- Matrix B: every solver with cleanup (2D, % above HK) -------------------------------------
def matrix_b():
    rows = []
    def bench(label, name):
        rows.append('| ' + label + ' | ' + ' | '.join(f1(mean(hk_gap(BR, name, grp(k, s)))) for k, s in COLS) + ' |')
    bench('greedy + 2-opt + Or-opt, Rn (baseline)', G)
    bench('greedy + 2-opt + Or-opt, R1', 'twoopt_or')
    bench('greedy + 2-opt, R1', 'twoopt')
    def boxes(label, names):
        vals = []
        for k, s in COLS:
            ms = grp(k, s)
            b = best_cell(HR, names, ms, 'full', vs_greedy)       # chosen by gap to greedy TO Rn
            vals.append('-' if b is None else f1(mean(hk_gap(HR, b[0], ms))) + ' (' + b[0].split(':')[-2 if b[0].startswith('box:') else -1] + ')')
        rows.append('| ' + label + ' | ' + ' | '.join(vals) + ' |')
    for c in CURVES:
        boxes(f'boxes, {c}', [f'box:{c}:{m}:t10' for m in MER])
    boxes('hull rings', HULL)
    save('matrix_B_with_cleanup.md', HEAD + '\n' + '\n'.join(rows))


# ---- Matrix C: every solver without cleanup (2D, % above HK) ----------------------------------
def matrix_c():
    rows = []
    def simple(label, D, name, level='full'):
        rows.append('| ' + label + ' | ' + ' | '.join(f1(mean(hk_gap(D, name, grp(k, s), level))) for k, s in COLS) + ' |')
    simple('nearest neighbour (nn)', BR, 'nn')
    simple('greedy (DIMACS greedy)', BR, 'greedy')
    simple('greedy start of our solver', HR, G, 'raw')
    simple('Hilbert order (spacefill)', BR, 'spacefill')
    simple('Moore order', BR, 'curve:moore')
    simple('Sierpinski order', BR, 'curve:sierpinski')
    for c in ['knopp', 'peano', 'morton', 'strip']:
        simple(f'{c.capitalize()} order', B2, f'curve:{c}')
    simple('hull rings in order', HR, 'curve:hull')
    def boxes(label, names):
        vals = []
        for k, s in COLS:
            ms = grp(k, s)
            b = best_cell(HR, names, ms, 'raw', hk_gap)           # shortest RAW tour
            vals.append('-' if b is None else f1(b[1]) + ' (' + b[0].split(':')[-2 if b[0].startswith('box:') else -1] + ')')
        rows.append('| ' + label + ' | ' + ' | '.join(vals) + ' |')
    for c in CURVES:
        boxes(f'boxes + merge, {c}', [f'box:{c}:{m}:t10' for m in MER])
    boxes('hull rings + merge', HULL)
    save('matrix_C_without_cleanup.md', HEAD + '\n' + '\n'.join(rows))


# ---- Matrix D: curve x merge, raw -> final (2D % above HK) ------------------------------------
def matrix_d():
    out = []
    for lab, (k, s) in [('uniform 10k', ('E', '10k')), ('clustered 10k', ('C', '10k')),
                        ('uniform 100k', ('E', '100k')), ('uniform 1M', ('E', '1M'))]:
        ms = [m for m in grp(k, s) if m in HK]
        out.append(f'### {lab} ({len(ms)} maps with a bound)\n\n| Curve | M1 | M2 | M3 | M4 | M5 |\n|---|---|---|---|---|---|')
        for c in CURVES:
            out.append(f'| {c} | ' + ' | '.join(
                f"{f1(mean(hk_gap(HR, f'box:{c}:{m}:t10', ms, 'raw')))} -> {f1(mean(hk_gap(HR, f'box:{c}:{m}:t10', ms)))}"
                for m in MER) + ' |')
        out.append('| hull (H0, M1..M5) | ' + ' ; '.join(
            f"{f1(mean(hk_gap(HR, h, ms, 'raw')))} -> {f1(mean(hk_gap(HR, h, ms)))}" for h in HULL) + ' | | | | |')
        out.append('')
    save('matrix_D_curve_by_merge_2d.md', '\n'.join(out))


# ---- Matrix A: cleanup levels (2D and 3D) -----------------------------------------------------
def level_len(name, m, level):
    """raw / T R1 / TO R1 / T Rn / TO Rn length of a 2D run, from the hull run file"""
    r = HR.get((name, m))
    if level == 'raw':
        return length(HR, name, m, 'raw')
    if level == 'TO Rn':
        return r['L'] if r else None
    if level == 'T Rn':
        return length(HR, name, m, 'T Rn')
    x = HR.get(('greedy:' + level if name == G else name + ':' + level, m))   # one-round (R1) runs
    return x['L'] if x else None


def matrix_a():
    LV = ['raw', 'T R1', 'TO R1', 'T Rn', 'TO Rn']
    fam = [('greedy + 2-opt + Or-opt', [G]),
           ('boxes, 5 main curves (25 settings)', [f'box:{c}:{m}:t10' for c in CURVES[:5] for m in MER]),
           ('boxes, Morton (5 settings)', [f'box:morton:{m}:t10' for m in MER]),
           ('hull rings (6 merges)', HULL)]
    rows = ['| Maps | Method | raw | T R1 | TO R1 | T Rn | TO Rn | share of raw excess removed |\n|---|---|---|---|---|---|---|---|']
    for kind, kn in [('E', 'uniform'), ('C', 'clustered')]:
        for size in ['10k', '100k']:
            ms = [m for m in grp(kind, size) if m in HK]
            for label, names in fam:
                names = [n for n in names if all(level_len(n, m, l) for m in ms for l in LV)]
                row = [mean([mean([(level_len(n, m, l) / HK[m] - 1) * 100 for m in ms]) for n in names]) for l in LV]
                rows.append(f'| {kn} {size} | {label} | ' + ' | '.join(f1(x) for x in row) +
                            f' | {(row[0] - row[4]) / row[0] * 100:.0f}% |')
    save('matrix_A_cleanup_levels_2d.md', '\n'.join(rows))


# ---- 3D: cleanup levels, raw starts, best merges ----------------------------------------------
maps3 = sorted({m for _, m in R3})
GROUPS3 = [('uniform 1k', r'E3d1k\.\d'), ('uniform 10k', r'E3d10k\.\d'), ('uniform 100k', r'E3d100k\.\d'),
           ('uniform 1M', r'E3d1M\.\d'), ('clustered 100k', r'C3d100k\.\d'),
           ('tight clustered 100k', r'CT3d100k\.\d'), ('star10k', r'star10k'), ('star250k', r'star250k'),
           ('gaia1', r'gaia1')]
CUR3 = ['hilbert3', 'moore3', 'sierpinski3', 'peano3', 'morton3', 'snake3']


def gm(p):
    return [m for m in maps3 if re.fullmatch(p, m)]


def vs3(name, ms, level):
    out = []
    for m in ms:
        r = R3.get((name, m))
        if not r:
            continue
        a = r['L'] if level == 'full' else (float(r['ex']['start_length']) if 'start_length' in r['ex'] else None)
        if a:
            out.append((a / R3[(G, m)]['L'] - 1) * 100)
    return out


def matrix_3d_raw():
    fmt = lambda x: '-' if x is None else (f'+{x:.0f}' if abs(x) < 10000 else f'+{x:,.0f}')
    t1 = ['| Maps | greedy start | ' + ' | '.join(f'{c} order' for c in CUR3) + ' |\n|---|' + '---|' * 7]
    t2 = ['| Maps | ' + ' | '.join(f'{c} + merge' for c in CUR3) + ' |\n|---|' + '---|' * 6]
    for lab, p in GROUPS3:
        ms = gm(p)
        g = [(R3[('greedy', m)]['L'] / R3[(G, m)]['L'] - 1) * 100 for m in ms if ('greedy', m) in R3]
        row = [fmt(mean(g))]
        for c in CUR3:
            v = vs3(f'curve:{c}', ms, 'full')
            row.append(fmt(mean(v)) if len(v) == len(ms) else '-')
        t1.append(f'| {lab} | ' + ' | '.join(row) + ' |')
        row = []
        for c in CUR3:
            best = None
            for m in MER:
                n = f'box:{c}:{m}:t10'
                v = vs3(n, ms, 'raw')
                if len(v) == len(ms) and (best is None or mean(v) < best[0]):
                    best = (mean(v), m)
            row.append('-' if best is None else f'{fmt(best[0])} ({best[1]})')
        t2.append(f'| {lab} | ' + ' | '.join(row) + ' |')
    save('matrix_C_3d.md', '\n'.join(t1) + '\n\n' + '\n'.join(t2))



# ---- B2: every solver with cleanup, % longer than greedy TO Rn (2D, includes C1M) -------------
def matrix_b2():
    cols = COLS[:4] + COLS[4:] + [('C', '1M')]
    head = ('| Method | uniform 1k | uniform 10k | uniform 100k | uniform 1M | clustered 1k | clustered 10k | '
            'clustered 100k | clustered 1M |\n|---|---|---|---|---|---|---|---|---|')
    rows = []
    for c in CURVES:
        vals = []
        for k, s_ in cols:
            ms = grp(k, s_)
            b = best_cell(HR, [f'box:{c}:{m}:t10' for m in MER], ms, 'full', vs_greedy)
            vals.append('-' if b is None else f1(b[1]))
        rows.append(f'| boxes, {c} | ' + ' | '.join(vals) + ' |')
    vals = []
    for k, s_ in cols:
        ms = grp(k, s_)
        b = best_cell(HR, HULL, ms, 'full', vs_greedy)
        vals.append('-' if b is None else f1(b[1]))
    rows.append('| hull rings | ' + ' | '.join(vals) + ' |')
    save('matrix_B2_vs_greedy_2d.md', head + '\n' + '\n'.join(rows))


# ---- A2: effort of the cleanup (moves per city, rounds, time, memory) -------------------------
def matrix_a_effort():
    def fx(n, ms, key):
        return mean([float(HR[(n, m)]['ex'][key]) for m in ms if (n, m) in HR and key in HR[(n, m)]['ex']])
    meths = [('greedy + 2-opt + Or-opt', G), ('boxes, Hilbert, M1', 'box:hilbert:M1:t10'),
             ('boxes, Hilbert, M4', 'box:hilbert:M4:t10'), ('boxes, Sierpinski, M5', 'box:sierpinski:M5:t10'),
             ('boxes, Morton, M2', 'box:morton:M2:t10'), ('boxes, Strip, M1', 'box:strip:M1:t10'),
             ('hull rings, M1', 'hull:M1')]
    rows = ['| Maps | Start | 2-opt moves per city | Or-opt moves per city | Rounds T / O | 2-opt time (s) | '
            'Or-opt stage time (s) | Total user time (s) | Peak memory (MB) |\n|---|---|---|---|---|---|---|---|---|']
    for kind, size, lab in [('E', '100k', 'uniform 100k'), ('E', '1M', 'uniform 1M')]:
        ms = grp(kind, size)
        n = HR[(G, ms[0])]['n']
        for label, name in meths:
            rows.append(f"| {lab} | {label} | {fx(name, ms, 'twoopt_moves') / n:.2f} | {fx(name, ms, 'oropt_moves') / n:.3f} | "
                        f"{fx(name, ms, 'twoopt_rounds'):.1f} / {fx(name, ms, 'oropt_rounds'):.1f} | "
                        f"{fx(name, ms, 'twoopt_time'):.2f} | {fx(name, ms, 'oropt_time'):.2f} | "
                        f"{mean([HR[(name, m)]['user'] for m in ms]):.1f} | {mean([HR[(name, m)]['mb'] for m in ms]):.0f} |")
    save('matrix_A_effort_2d.md', '\n'.join(rows))


# ---- 2D scaling: user time and peak memory, 1k to 1M -------------------------------------------
def scaling_2d():
    mm = [('greedy + cleanup', G), ('boxes, Hilbert, M1', 'box:hilbert:M1:t10'), ('boxes, Hilbert, M4', 'box:hilbert:M4:t10'),
          ('boxes, Sierpinski, M5', 'box:sierpinski:M5:t10'), ('boxes, Morton, M2', 'box:morton:M2:t10'),
          ('boxes, Strip, M1', 'box:strip:M1:t10'), ('hull rings, M1', 'hull:M1')]
    rows = ['| Method | 1k time | 10k time | 100k time | 1M time | 1k MB | 10k MB | 100k MB | 1M MB |\n|---|---|---|---|---|---|---|---|---|']
    for lab, n in mm:
        t, b = [], []
        for k, s_ in COLS[:4]:
            ms = grp(k, s_)
            t.append(mean([HR[(n, m)]['user'] for m in ms if (n, m) in HR]))
            b.append(mean([HR[(n, m)]['mb'] for m in ms if (n, m) in HR]))
        rows.append(f'| {lab} | ' + ' | '.join(f'{x:.2f}' if x < 10 else f'{x:.1f}' for x in t) + ' | ' +
                    ' | '.join(f'{x:.0f}' for x in b) + ' |')
    save('scaling_2d.md', '\n'.join(rows))


# ---- 3D: cleanup levels, curve x merge (raw -> final), time and memory --------------------------
def level3(name, m, level):
    if name == G:
        nm = {'raw': 'greedy', 'T R1': 'twoopt', 'TO R1': 'twoopt_or', 'T Rn': 'greedy+2opt', 'TO Rn': G}[level]
        r = R3.get((nm, m))
        return r['L'] if r else None
    r = R3.get((name, m))
    if level == 'raw':
        return float(r['ex']['start_length']) if r and 'start_length' in r['ex'] else None
    if level == 'TO Rn':
        return r['L'] if r else None
    if level == 'T Rn':
        return float(r['ex']['after_2opt_length']) if r and 'after_2opt_length' in r['ex'] else None
    x = R3.get((name + ':' + level, m))
    return x['L'] if x else None


def matrix_a3():
    LV = ['raw', 'T R1', 'TO R1', 'T Rn', 'TO Rn']
    fmt = lambda x: '-' if x is None else (f'{x:+.1f}' if abs(x) < 1000 else f'{x:+,.0f}')
    rows = ['| Maps | Start | raw | T R1 | TO R1 | T Rn | TO Rn |\n|---|---|---|---|---|---|---|']
    for lab, p in [('uniform 100k', r'E3d100k\.\d'), ('clustered 100k', r'C3d100k\.\d'), ('tight clustered 100k', r'CT3d100k\.\d')]:
        ms = gm(p)
        for ml, name in [('greedy + 2-opt + Or-opt', G), ('boxes, hilbert3, M1', 'box:hilbert3:M1:t10'),
                         ('boxes, hilbert3, M5', 'box:hilbert3:M5:t10'), ('boxes, snake3, M2', 'box:snake3:M2:t10')]:
            vals = []
            for l in LV:
                v = [(level3(name, m, l) / R3[(G, m)]['L'] - 1) * 100 for m in ms if level3(name, m, l)]
                vals.append(mean(v) if len(v) == len(ms) else None)
            rows.append(f'| {lab} | {ml} | ' + ' | '.join(fmt(x) for x in vals) + ' |')
    save('matrix_A_cleanup_levels_3d.md', '\n'.join(rows))


def matrix_b3_d3():
    allg = [('uniform 1k', r'E3d1k\.\d'), ('uniform 10k', r'E3d10k\.\d'), ('uniform 100k', r'E3d100k\.\d'),
            ('uniform 1M', r'E3d1M\.\d'), ('clustered 1k', r'C3d1k\.\d'), ('clustered 10k', r'C3d10k\.\d'),
            ('clustered 100k', r'C3d100k\.\d'), ('clustered 1M', r'C3d1M\.\d'), ('tight clustered 1k', r'CT3d1k\.\d'),
            ('tight clustered 10k', r'CT3d10k\.\d'), ('tight clustered 100k', r'CT3d100k\.\d'),
            ('tight clustered 1M', r'CT3d1M\.\d'), ('star1k', r'star1k'), ('star10k', r'star10k'), ('kj37859', r'kj37859'),
            ('hyg109399', r'hyg109399'), ('hyg119614', r'hyg119614'), ('star250k', r'star250k'), ('gaia1', r'gaia1')]
    rows = ['| Maps | ' + ' | '.join(CUR3) + ' |\n|---|' + '---|' * 6]
    for lab, p in allg:
        ms = gm(p)
        cells = []
        for c in CUR3:
            best = None
            for m in MER:
                v = vs3(f'box:{c}:{m}:t10', ms, 'full')
                if len(v) == len(ms) and (best is None or mean(v) < best):
                    best = mean(v)
            cells.append('-' if best is None else f'{best:.1f}')
        rows.append(f'| {lab} | ' + ' | '.join(cells) + ' |')
    save('matrix_B3_3d.md', '\n'.join(rows))
    out = []
    fm = lambda x: '-' if x is None else (f'{x:+.1f}' if abs(x) < 100 else f'{x:+,.0f}')
    for lab, p in [('uniform 100k', r'E3d100k\.\d'), ('clustered 100k', r'C3d100k\.\d'),
                   ('tight clustered 100k', r'CT3d100k\.\d'), ('star250k', r'star250k')]:
        ms = gm(p)
        out.append(f'### {lab}\n\n| Curve | M1 | M2 | M3 | M4 | M5 |\n|---|---|---|---|---|---|')
        for c in CUR3:
            cells = []
            for m in MER:
                n = f'box:{c}:{m}:t10'
                raw, fin = vs3(n, ms, 'raw'), vs3(n, ms, 'full')
                cells.append(f'{fm(mean(raw))} -> {fm(mean(fin))}' if len(raw) == len(ms) else '-')
            out.append(f'| {c} | ' + ' | '.join(cells) + ' |')
        out.append('')
    save('matrix_D_curve_by_merge_3d.md', '\n'.join(out))


def scaling_3d():
    sel = [('uniform 10k', r'E3d10k\.\d'), ('uniform 100k', r'E3d100k\.\d'), ('uniform 1M', r'E3d1M\.\d'),
           ('tight clustered 100k', r'CT3d100k\.\d'), ('star250k', r'star250k'), ('gaia1', r'gaia1')]
    rows = ['| Method | ' + ' | '.join(l for l, _ in sel) + ' |\n|---|' + '---|' * len(sel)]
    for what, key, f in [('time (s)', 'user', '{:.1f}'), ('peak MB', 'mb', '{:,.0f}')]:
        for lab, n in [('greedy + cleanup', G), ('boxes, hilbert3, M1', 'box:hilbert3:M1:t10'),
                       ('boxes, hilbert3, M4', 'box:hilbert3:M4:t10'), ('boxes, hilbert3, M5', 'box:hilbert3:M5:t10'),
                       ('boxes, snake3, M2', 'box:snake3:M2:t10')]:
            cells = []
            for _, p in sel:
                v = [R3[(n, m)][key] for m in gm(p) if (n, m) in R3]
                cells.append('-' if not v else f.format(mean(v)))
            rows.append(f'| {lab}, {what} | ' + ' | '.join(cells) + ' |')
    save('scaling_3d.md', '\n'.join(rows))


# ---- TSPLIB, by type: mean % above the proven optimum (third column of refs/bounds) ------------
def tsplib():
    opt = {}
    for line in open(ROOT / 'refs' / 'bounds'):
        p = line.split()
        if len(p) >= 3:
            opt[p[0]] = float(p[2])
    T = {'Cities and towns': 'usa13509 brd14051 d15112 d18512 fnl4461 nrw1379',
         'Circuit-board drilling': 'd1291 d1655 d2103 fl1400 fl1577 fl3795 pcb1173 pcb3038 u1060 u1432 u1817 u2152 u2319 rl5915',
         'Chip layout': 'pla7397 pla33810 pla85900',
         'Origin not stated': 'rl1304 rl1323 rl1889 rl5934 rl11849 pr1002 pr2392 vm1084 vm1748'}
    rows = ['| Method | ' + ' | '.join(f'{k} ({len(v.split())})' for k, v in T.items()) + ' |\n|---|---|---|---|---|']
    for lab, n in [('nearest neighbour', 'nn'), ('greedy', 'greedy'), ('Hilbert order', 'spacefill'),
                   ('Moore order', 'curve:moore'), ('Sierpinski order', 'curve:sierpinski'),
                   ('greedy + 2-opt, R1', 'twoopt'), ('greedy + 2-opt + Or-opt, R1', 'twoopt_or'),
                   ('greedy + 2-opt + Or-opt, Rn', G), ('boxes, Sierpinski, M5, t=4', 'box:sierpinski:M5:t4'),
                   ('boxes, Hilbert, M5, t=10', 'box:hilbert:M5:t10'), ('boxes, Sierpinski, M1, t=10', 'box:sierpinski:M1:t10'),
                   ('boxes, Hilbert, M1, t=10', 'box:hilbert:M1:t10'), ('boxes, Knopp, M1, t=10', 'box:knopp:M1:t10'),
                   ('boxes, Sierpinski, M2, t=12', 'box:sierpinski:M2:t12')]:
        cells = []
        for ms in T.values():
            v = [(BR[(n, m)]['L'] / opt[m] - 1) * 100 for m in ms.split() if (n, m) in BR]
            cells.append(f'{mean(v):.2f}' if len(v) == len(ms.split()) else '-')
        rows.append(f'| {lab} | ' + ' | '.join(cells) + ' |')
    save('tsplib_by_type.md', '\n'.join(rows))


# ---- Spread: the same headline tables with the sample standard deviation across maps --------
def sd(v):
    return st.stdev(v) if len(v) >= 2 else None


def ms_(x, s, big=1000):
    """'mean +- sd'; the sd is left out when only one map has a bound"""
    if x is None:
        return '-'
    if s is None:
        return f1(x, big)
    return f'{f1(x, big)} ± ' + ('<0.1' if s < 0.05 else f'{s:.1f}')


def sd_headline():
    nmaps = ' | '.join(str(len([m for m in grp(k, s) if m in HK])) for k, s in COLS)
    # B1 with sd
    rows = []
    def bench(label, name):
        cells = []
        for k, s in COLS:
            v = hk_gap(BR, name, grp(k, s))
            cells.append(ms_(mean(v), sd(v)))
        rows.append('| ' + label + ' | ' + ' | '.join(cells) + ' |')
    bench('greedy + 2-opt + Or-opt, Rn (baseline)', G)
    bench('greedy + 2-opt + Or-opt, R1', 'twoopt_or')
    bench('greedy + 2-opt, R1', 'twoopt')
    def boxes(label, names):
        cells = []
        for k, s in COLS:
            ms = grp(k, s)
            b = best_cell(HR, names, ms, 'full', vs_greedy)
            if b is None:
                cells.append('-')
                continue
            v = hk_gap(HR, b[0], ms)
            cells.append(ms_(mean(v), sd(v)) + ' (' + b[0].split(':')[-2 if b[0].startswith('box:') else -1] + ')')
        rows.append('| ' + label + ' | ' + ' | '.join(cells) + ' |')
    for c in CURVES:
        boxes(f'boxes, {c}', [f'box:{c}:{m}:t10' for m in MER])
    boxes('hull rings', HULL)
    rows.append('| maps with a bound (n) | ' + nmaps + ' |')
    save('matrix_B_with_cleanup_sd.md', HEAD + '\n' + '\n'.join(rows))
    # C1 with sd
    rows = []
    def simple(label, D, name, level='full'):
        cells = []
        for k, s in COLS:
            v = hk_gap(D, name, grp(k, s), level)
            cells.append(ms_(mean(v), sd(v)))
        rows.append('| ' + label + ' | ' + ' | '.join(cells) + ' |')
    simple('nearest neighbour (nn)', BR, 'nn')
    simple('greedy (DIMACS greedy)', BR, 'greedy')
    simple('greedy start of our solver', HR, G, 'raw')
    simple('Hilbert order (spacefill)', BR, 'spacefill')
    simple('Moore order', BR, 'curve:moore')
    simple('Sierpinski order', BR, 'curve:sierpinski')
    for c in ['knopp', 'peano', 'morton', 'strip']:
        simple(f'{c.capitalize()} order', B2, f'curve:{c}')
    simple('hull rings in order', HR, 'curve:hull')
    def boxes_raw(label, names):
        cells = []
        for k, s in COLS:
            ms = grp(k, s)
            b = best_cell(HR, names, ms, 'raw', hk_gap)
            if b is None:
                cells.append('-')
                continue
            v = hk_gap(HR, b[0], ms, 'raw')
            cells.append(ms_(mean(v), sd(v)) + ' (' + b[0].split(':')[-2 if b[0].startswith('box:') else -1] + ')')
        rows.append('| ' + label + ' | ' + ' | '.join(cells) + ' |')
    for c in CURVES:
        boxes_raw(f'boxes + merge, {c}', [f'box:{c}:{m}:t10' for m in MER])
    boxes_raw('hull rings + merge', HULL)
    rows.append('| maps with a bound (n) | ' + nmaps + ' |')
    save('matrix_C_without_cleanup_sd.md', HEAD + '\n' + '\n'.join(rows))
    # At a glance (A1 families), 100k maps: sd over maps of the family-mean gap
    LV = ['raw', 'TO Rn']
    fam = [('Greedy', [G]),
           ('Boxes (25 settings)', [f'box:{c}:{m}:t10' for c in CURVES[:5] for m in MER]),
           ('Hull rings (6 merges)', HULL)]
    out = ['| Maps | ' + ' | '.join(f for f, _ in fam) + ' |\n|---|---|---|---|']
    for kind, kn in [('E', 'uniform'), ('C', 'clustered')]:
        ms = [m for m in grp(kind, '100k') if m in HK]
        cells = []
        for label, names in fam:
            names = [n for n in names if all(level_len(n, m, l) for m in ms for l in LV)]
            per_map = {l: [mean([(level_len(n, m, l) / HK[m] - 1) * 100 for n in names]) for m in ms] for l in LV}
            cells.append(f"{f1(mean(per_map['raw']))} → {ms_(mean(per_map['TO Rn']), sd(per_map['TO Rn']))}")
        out.append(f'| {kn} 100k (n = {len(ms)}) | ' + ' | '.join(cells) + ' |')
    save('matrix_at_a_glance_sd.md', '\n'.join(out))


# ---- Checks added after review: 1M on three maps, and seam share against the t sweep ------------
def one_m_three_maps():
    """B4: box settings that ran on all three 1M maps of each type, % longer than greedy TO Rn (mean +- sd)."""
    names = [('Sierpinski, M1, t=10', 'box:sierpinski:M1:t10'), ('Knopp, M1, t=10', 'box:knopp:M1:t10'),
             ('Sierpinski, M2, t=12', 'box:sierpinski:M2:t12'), ('Hilbert, M1, t=10', 'box:hilbert:M1:t10'),
             ('Sierpinski, M5, t=4', 'box:sierpinski:M5:t4'), ('Hilbert, M5, t=10', 'box:hilbert:M5:t10')]
    rows = ['| Setting | uniform 1M (3 maps) | clustered 1M (3 maps) |\n|---|---|---|']
    for label, a in names:
        cells = []
        for kind in ('E1M', 'C1M'):
            v = [(B2[(a, f'{kind}.{k}')]['L'] / B2[(G, f'{kind}.{k}')]['L'] - 1) * 100
                 for k in range(3) if (a, f'{kind}.{k}') in B2 and (G, f'{kind}.{k}') in B2]
            cells.append(f'+{mean(v):.2f} ± {sd(v):.2f} (n={len(v)})')
        rows.append(f'| {label} | ' + ' | '.join(cells) + ' |')
    save('matrix_B4_1M_three_maps.md', '\n'.join(rows))


def seam_vs_t():
    """Seam share of the start tour (join_cost / start_length) and the finished gap, Hilbert M1, t = 4..12, 10k maps."""
    rows = ['| t | uniform 10k: seam share of start | uniform 10k: % above HK | clustered 10k: seam share of start | clustered 10k: % above HK |\n|---|---|---|---|---|']
    for t in (4, 6, 8, 10, 12):
        a = f'box:hilbert:M1:t{t}'
        cells = []
        for kind in ('E', 'C'):
            ms = grp(kind, '10k')
            js = [float(B2[(a, m)]['ex']['join_share']) for m in ms if (a, m) in B2]
            fg = [(B2[(a, m)]['L'] / HK[m] - 1) * 100 for m in ms if (a, m) in B2 and m in HK]
            cells += [f'{mean(js):.2f}', f'{mean(fg):.2f} ± {sd(fg):.2f}']
        rows.append(f'| {t} | ' + ' | '.join(cells) + ' |')
    save('table_seam_vs_t.md', '\n'.join(rows))


def seam_correlations():
    """Across every box setting (curve, merge, t) run on the 10k maps: does the seam share or the start gap predict the finished gap?"""
    import math
    data = collections.defaultdict(list)
    for D in (B2, HR):
        for (a, m), r in D.items():
            if not re.fullmatch(r'box:\w+:M\d:t\d+', a) or not re.fullmatch(r'[EC]10k\.\d', m) or m not in HK:
                continue
            if 'join_share' not in r['ex'] or 'start_length' not in r['ex'] or (a, m[0], m) in data:
                continue
            data[(a, m[0], m)] = (float(r['ex']['join_share']), (float(r['ex']['start_length']) / HK[m] - 1) * 100,
                                  (r['L'] / HK[m] - 1) * 100)
    per = collections.defaultdict(list)
    for (a, k, m), v in data.items():
        per[(a, k)].append(v)
    def corr(x, y):
        mx, my = mean(x), mean(y)
        return sum((p - mx) * (q - my) for p, q in zip(x, y)) / math.sqrt(sum((p - mx) ** 2 for p in x) * sum((q - my) ** 2 for q in y))
    out = ['| maps | settings | corr(seam share, finished gap) | corr(start gap, finished gap) |\n|---|---|---|---|']
    for k, kn in (('E', 'uniform 10k'), ('C', 'clustered 10k')):
        ks = [x for x in per if x[1] == k and len(per[x]) >= 2]
        j = [mean([v[0] for v in per[x]]) for x in ks]
        sg = [mean([v[1] for v in per[x]]) for x in ks]
        fg = [mean([v[2] for v in per[x]]) for x in ks]
        out.append(f'| {kn} | {len(ks)} | {corr(j, fg):+.2f} | {corr(sg, fg):+.2f} |')
    save('table_seam_correlations.md', '\n'.join(out))


def best3d():
    """Best box setting (curve, merge, t = 10) after cleanup for each 3D map group: mean % longer than greedy TO Rn."""
    groups = [('uniform 1k', r'E3d1k\.\d'), ('uniform 10k', r'E3d10k\.\d'), ('uniform 100k', r'E3d100k\.\d'),
              ('uniform 1M', r'E3d1M\.\d'), ('clustered 1k', r'C3d1k\.\d'), ('clustered 10k', r'C3d10k\.\d'),
              ('clustered 100k', r'C3d100k\.\d'), ('clustered 1M', r'C3d1M\.\d'),
              ('tight clustered 1k', r'CT3d1k\.\d'), ('tight clustered 10k', r'CT3d10k\.\d'),
              ('tight clustered 100k', r'CT3d100k\.\d'), ('tight clustered 1M', r'CT3d1M\.\d'),
              ('star1k', r'star1k'), ('star10k', r'star10k'), ('kj37859', r'kj37859'), ('hyg109399', r'hyg109399'),
              ('hyg119614', r'hyg119614'), ('star250k', r'star250k'), ('gaia1', r'gaia1')]
    rows = ['| Maps | Best setting | % longer than greedy TO Rn |\n|---|---|---|']
    for lab, p in groups:
        ms = gm(p)
        best = None
        for c in CUR3:
            for m in MER:
                v = vs3(f'box:{c}:{m}:t10', ms, 'full')
                if len(v) == len(ms) and (best is None or mean(v) < best[0]):
                    best = (mean(v), f'{c} {m}')
        rows.append(f'| {lab} | {best[1]} | {best[0]:+.2f} |')
    save('table_3d_best_setting.md', '\n'.join(rows))


if __name__ == '__main__':
    matrix_b(); matrix_b2(); matrix_c(); matrix_d(); matrix_a(); matrix_a_effort(); scaling_2d()
    matrix_3d_raw(); matrix_a3(); matrix_b3_d3(); scaling_3d(); tsplib(); sd_headline(); one_m_three_maps(); seam_vs_t(); seam_correlations(); best3d()
