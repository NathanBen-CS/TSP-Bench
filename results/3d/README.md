# results/3d

Results of the 3D box study. The write-up is [docs/box_tsp_3d.md](../../docs/box_tsp_3d.md).
Every run was made on one laptop, pinned to CPU core 0 unless noted, and every tour was
checked by `bin/length3d.exe`.

**Code names** (as in 2D): a bare name = that order only; B = boxes (exact loops +
merge); T = 2-opt; O = Or-opt; R1 = one round of cleanup; Rn = rounds until nothing
improves. In the CSV files the setting `box:hilbert3:M5:t10` is hilbert3 BTO Rn M5 t=10;
a suffix `:TO R1` or `:T R1` is the one-round cleanup, `:k24` a neighbour-list size other
than 16, and `kd:` instead of `box:` the kd-median boxes. `greedy+2opt+oropt` is greedy
TO Rn, `twoopt_or` greedy TO R1, `greedy+2opt` greedy T Rn, `twoopt` greedy T R1.

## Data

| file | what it holds |
|---|---|
| `runs.csv` | one row per timed run (6,360 runs, all OK): length, user / kernel / wall seconds, peak memory, and the program's own stage times and counts in `extras` |
| `runs_untimed.csv` | 24 runs at 1M cities and on gaia1 (morton3 and peano3 with M2–M4), run six at a time on other cores: **lengths only**, their times are not comparable and are not used |
| `run.log` | console log of all runs |
| `summary.md` | every table, made by `bench/box3_report.py` from the two CSV files |
| `summary.csv` | the same numbers, one row per map group and setting |
| `diagnostics.csv` | box sizes, neighbours outside the box and tour edges crossing box borders, 2D and 3D, made by `bench/box3_diag.py` |

Bounds used for the percentages are in `refs/`: `refs/star_bounds.csv` (published
Waterloo optimal or best tours) and `refs/hk_own.csv` (**our own** Held-Karp 1-tree
bound, not a published one).

## Charts

Made by `bench/box3_charts.py`. Unless a chart says otherwise, the measure is % longer
than greedy TO Rn on the same map, mean of 3 seeds (stars: one map each).

| file | the question it answers |
|---|---|
| `box3_curve_merge.png` | which curve and merge do best, at 100k and on the stars, after full cleanup and before any |
| `box3_overview.png` | how far every baseline and the best box settings are from the bound (100k) and from the best known tour (stars) |
| `box3_uniform_vs_clustered.png` | does the gap change with the kind of map |
| `box3_t_sweep.png` | does the box size t matter (hilbert3, 10k cities, t = 4 to 12) |
| `box3_t_sweep_sierpinski3.png` | the same for sierpinski3, whose boxes come close to t |
| `box3_k_sweep.png` | does the neighbour-list size K matter (100k cities and two star maps, K = 8 to 32; hilbert3 solid, sierpinski3 M5 dashed) |
| `box3_k_time.png` | what a larger K costs in time |
| `box3_scale.png` | does the gap change from 1k to 1M cities |
| `box3_1M.png` | every setting run at 1M cities and on gaia1, with user seconds |
| `box3_kd.png` | kd-median boxes (exactly t cities) against curve-cell boxes |
