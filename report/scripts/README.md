# Report scripts

`tables.py` rebuilds every matrix table of the project report from the run files:

```
py report/scripts/tables.py      # writes report/data/*.md
```

It reads `results/2d/runs.csv`, `results/2d/box2/runs.csv`, `results/2d/hull/runs.csv`,
`results/3d/runs.csv` and `refs/bounds`, and writes only into `report/data/`.

| Output | Report table |
|---|---|
| `matrix_A_cleanup_levels_2d.md`, `matrix_A_effort_2d.md`, `matrix_A_cleanup_levels_3d.md` | A1, A2, A3 (what cleanup does) |
| `matrix_B_with_cleanup.md`, `matrix_B2_vs_greedy_2d.md`, `matrix_B3_3d.md` | B1, B2, B3 (every solver with cleanup) |
| `matrix_C_without_cleanup.md`, `matrix_C_3d.md` | C1 and C1b, C3 and C4 (every solver without cleanup) |
| `matrix_D_curve_by_merge_2d.md`, `matrix_D_curve_by_merge_3d.md` | D1 to D7 (curve by merge, raw to final) |
| `scaling_2d.md`, `scaling_3d.md` | E1, E2, E4 (time and memory) |
| `tsplib_by_type.md` | F1 (real-world maps by type) |
| `matrix_B_with_cleanup_sd.md`, `matrix_C_without_cleanup_sd.md`, `matrix_at_a_glance_sd.md` | B1, C1 and the at-a-glance table with the sample standard deviation over maps (`±`; none when only one map has a bound) |

| `matrix_B4_1M_three_maps.md`, `table_seam_vs_t.md`, `table_seam_correlations.md` | B4 (six box settings on all three 1M maps per type), and the join-share tests in section 7 (Hilbert M1 over t = 4 to 12, and the correlations over 63 settings) |

Not rebuilt by it: C2 (the 2D raw table against greedy plus cleanup), E3 (the stage breakdown at
1M) and G1, G2 (the summary views), which were computed with one-off commands from the same files,
and the figures in `figures/`. The main-result figure `figures/fig_hero.png` is drawn by `fig_hero.py`
(`py report/scripts/fig_hero.py`), from the same run files.
