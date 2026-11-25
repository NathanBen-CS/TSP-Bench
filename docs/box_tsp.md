# Box-TSP: build notes and results

A Karp-style partitioning heuristic (boxes + exact solve per box + glue +
cleanup), added to the DIMACS-style Euclidean TSP bench and compared with the
bench's existing solvers. **It has no worst-case guarantee, and it is not a
PTAS.**

**Short answer:** the full method works and produces valid tours at every size.
But it is **worse than the same cleanup started from a greedy tour**, by about
3–4% on uniform data and about 3–8% on clustered data, and it is not faster. The
box size *t* has almost no effect after cleanup.

Full tables: [`results/2d/box/summary.md`](../results/2d/box/summary.md) ·
raw runs: `results/2d/box/runs.csv` · summary numbers: `results/2d/box/summary.csv`.

---

## 1. What the bench already had (Task 0)

1. **Language:** C solvers; Python stdlib runner (`bench/run.py`) on Windows, with gcc from MSYS2.
2. **Registering a solver:** an entry in a config JSON (`exe`, `max_n`, `timeout_s`, optional `geometric_only`, now also `args`).
3. **Contract:** `solver <instance> <tour_out> [seed]`; the tour is the city numbers 1..N, one per line; `key=value` lines on stdout become extra columns.
4. **Input format:** TSPLIB (`EUC_2D`, `CEIL_2D`, `UPPER_DIAG_ROW`), read by the DIMACS reader `codes/Instance.c`. Distances are rounded to integers (nint for EUC_2D).
5. **Validation:** `codes/length.exe` checks that the tour is a permutation and recomputes its length.
6. **Scoring:** user CPU time and peak memory of the whole process (Windows API), pinned to one CPU; one row per run.
7. **Reference values:** `refs/bounds` (DIMACS: Held-Karp bound, plus optimal/best-known where available).
8. **Existing solvers:** `greedy` (DIMACS), `nn`, `spacefill` (Hilbert order), `twoopt`, `twoopt_or` (greedy start, 2-opt / 2-opt + Or-opt).
9. **Shared pieces already there:** Hilbert sort, grid neighbour lists (K = 16), greedy start, 2-opt and Or-opt with a work queue, shorter-side reversal and a reversal cap.
10. **Generators:** `portgen` (uniform, integer, 10⁶ square) and `portcgen` (DIMACS clustered) already match the plan, so they were reused.

The new code follows the bench's conventions, not the prototype's: integer
distances, the bench's neighbour lists and cleanup, and K = 16.

## 2. What was built

| Piece | Where | Notes |
|---|---|---|
| Hilbert sort on a 2²⁰ grid | `solvers/d2/hilbert.c` (`hilbertSort`) | factored out of the existing code; `spacefill` gives byte-identical tours |
| Quadtree boxes | `solvers/common/boxes.c` (`boxSplit`) | boxes = contiguous runs of the sorted list, as in `proto.c`; depth cap 20 |
| Exact loop per box | `boxes.c` (`heldKarpLoop`) | Held-Karp, table 2^(t−1) × (t−1); t ≤ 12 |
| Glue | `boxes.c` (`boxTour`) | the `proto.c` method; records G and the sum of the box loops separately |
| Cleanup | `solvers/common/localsearch.c` (`localSearchToOptimum`) | the bench's 2-opt / Or-opt, repeated until a whole round with every city re-examined finds nothing |
| Solver and variants | `solvers/d2/boxtsp.c` | one program; `--start`, `--cleanup`, `--t`, `--k` |
| Tests | `tests/boxtest2d.c` | 33 tests, all pass |
| Data | `bench/box_instances.py`, `bench/check_clumpiness.py` | 38 instances; the clumpiness test passes |
| Experiment | `bench/box_config.json`, `bench/box_report.py` | 13 algorithms × 38 instances |
| Validator | `bench/run.py` | new `MISMATCH` status when the solver's own `tour_cost` ≠ `length.exe` |

**Deviations from `proto.c`:**
- **Distances:** the bench's integer distances, so lengths match `length.exe`.
- **Depth-capped boxes:** a box at the depth cap with more than t cities keeps its curve order. In `proto.c` it would overflow the DP table. It never happened on the experiment data, only in the duplicate-point tests.
- **Stronger 2-opt:** best move per city, 16 neighbours, rounds until clean.

## 3. Tests (Task 2): `bin/boxtest2d.exe`, 33 / 33 pass

- **Hilbert:** on a 16×16 grid, consecutive curve positions are neighbouring cells; the 2×2 order is (0,0), (0,1), (1,1), (1,0).
- **Boxes:** uniform, clustered and duplicate-heavy data, t = 6, 8, 10, 12. Every city is in exactly one box. Every box is ≤ t, except depth-capped boxes, which are checked to hold only cities in the same grid cell.
- **DP:** equals brute force on 300 random sets of 4–9 cities (some with duplicate points), and the returned loop has exactly that length.
- **Glue:** the tour is a permutation, and its length ≤ loops + G.
- **2-opt:** on 300 cities with complete neighbour lists, the tour never gets longer, and **brute force over all edge pairs finds no improving move**.
- **Or-opt:** never longer; no 2-opt move left; no improving Or-opt move left *within the gain criterion* the search uses. Brute force still finds 4–8 improving Or-opt moves *outside* that criterion, which the search deliberately never tries.
- **Whole pipeline:** 20,000 cities each of uniform, clustered and heavy duplicates (25 cities per spot). Valid tours, nothing crashes or loops.

The existing solvers were re-run after the shared-code changes. Their tours are
**byte-identical** to before (24 of 24 checked).

## 4. Data (Task 1)

38 instances in `C:\Users\Natha\tsp-instances-box`, kept apart from the DIMACS
testbed:

| | 1k | 10k | 100k | 1M |
|---|---|---|---|---|
| uniform `E` | 10 seeds | 3 | 3 | 3 |
| clustered `C` | 10 | 3 | 3 | 3 |

- **Sources:** 29 are the DIMACS testbed files, hard-linked. 9 were generated with the same generators and the README's seed pattern: E100k.2, C100k.2, E1M.1–2, C1M.0–2.
- **Clustered generator:** `portcgen` doesn't clip points to the square. The bench's version was kept.

**Clumpiness check** (`bench/check_clumpiness.py`, run automatically):

| size | empty cells, uniform → clustered | mean nearest-neighbour distance, uniform → clustered | ratio |
|---|---|---|---|
| 1k | 0.904 → 0.930 | 16,314 → 6,639 | 0.41 |
| 10k | 0.362 → 0.820 | 5,115 → 2,038 | 0.40 |
| 100k | 0.000 → 0.585 | 1,578 → 680 | 0.43 |
| 1M | 0.000 → 0.125 | 493 → 217 | 0.44 |

The test passes at every size (it needs a ratio < 0.7 and more empty cells). At 1M
the empty-cell difference is small, because 10,000 clusters cover most of the map.
The nearest-neighbour measure separates the two types clearly at every size.

## 5. Calibration

Full method on 100,000 uniform cities, gap vs 0.7124·√(N·A), mean of 3 seeds:

| t | box tour | + 2-opt (ours) | prototype (box + its 2-opt) | + Or-opt (ours) |
|---|---|---|---|---|
| 6 | +28.97% | +15.53% | +16.8% | +8.13% |
| 8 | +28.11% | +15.15% | +16.3% | +8.08% |
| 10 | +27.65% | +14.87% | +16.0% | +8.02% |
| 12 | +27.46% | +14.58% | +15.7% | +7.95% |

This is the same region as the prototype and follows the same trend. Ours is about
1 point better after 2-opt, which fits a stronger 2-opt. Before cleanup, the box
tour is about 28% over, matching the plan's "real gap near 28%".

## 6. Results (Task 4)

Mean ± sample standard deviation over seeds.

**References:**
- **Uniform:** 0.7124·√(N·A).
- **Clustered:** the Held-Karp bound, which is available for all 1k and 10k seeds and 2 of the 3 100k seeds. None exists for 1M.
- **"vs baseline":** % longer than `greedy+2opt+oropt` on the same instance, which works for every seed.

**Gap, uniform:**

| algorithm | 1k | 10k | 100k | 1M |
|---|---|---|---|---|
| greedy (DIMACS) | +20.14 ± 2.80 | +16.72 ± 0.55 | +14.20 ± 0.22 | +13.63 ± 0.09 |
| hilbert-only | +41.19 ± 1.27 | +39.01 ± 0.70 | +37.90 ± 0.07 | +37.52 ± 0.02 |
| boxes-only | +31.26 ± 1.59 | +29.63 ± 0.41 | +27.65 ± 0.18 | +27.24 ± 0.01 |
| boxes+2opt | +15.48 ± 1.50 | +15.29 ± 0.48 | +14.87 ± 0.16 | +15.39 ± 0.08 |
| **boxes+2opt+oropt** | +9.44 ± 1.40 | +8.59 ± 0.24 | +8.02 ± 0.28 | +8.11 ± 0.05 |
| **greedy+2opt+oropt** | **+6.00 ± 1.60** | **+4.48 ± 0.30** | **+3.65 ± 0.21** | **+4.46 ± 0.23** |
| twoopt_or (existing) | +6.01 ± 1.34 | +4.45 ± 0.33 | +3.69 ± 0.17 | +4.57 ± 0.26 |

At 1k the 0.7124 estimate is too low (cities near the border have fewer
neighbours), so all 1k uniform gaps look a few points worse than they are. The
comparisons between algorithms aren't affected.

**Gap, clustered (vs Held-Karp; number of seeds with a bound in brackets):**

| algorithm | 1k (10) | 10k (3) | 100k (2) |
|---|---|---|---|
| greedy (DIMACS) | +19.52 ± 3.32 | +23.09 ± 0.15 | +21.74 ± 0.58 |
| hilbert-only | +48.01 ± 3.83 | +55.55 ± 2.94 | +55.81 ± 0.03 |
| boxes-only | +39.32 ± 3.57 | +46.51 ± 2.60 | +46.05 ± 0.23 |
| boxes+2opt | +16.54 ± 3.73 | +24.31 ± 1.88 | +24.94 ± 0.18 |
| **boxes+2opt+oropt** | +10.45 ± 2.77 | +17.61 ± 1.72 | +17.59 ± 0.03 |
| **greedy+2opt+oropt** | **+7.45 ± 2.40** | **+9.16 ± 1.25** | **+9.51 ± 0.73** |

**Full method vs greedy baseline, paired (% longer; all seeds):**

| | 1k | 10k | 100k | 1M |
|---|---|---|---|---|
| uniform | +3.25 ± 0.49 | +3.94 ± 0.32 | +4.21 ± 0.45 | +3.50 ± 0.28 |
| clustered | +2.87 ± 4.35 | +7.76 ± 2.80 | +7.06 ± 0.77 | +6.79 ± 0.13 |

**Time (user seconds, uniform):** at 100k the full method takes 1.57 s and the
greedy baseline 1.10 s. At 1M they take 23.1 s and 20.2 s. On clustered 1M data,
26.6 s and 14.9 s.

**Where the full method's time goes (1M uniform):**

| stage | time |
|---|---|
| Hilbert sort | 0.71 s |
| boxes | 0.02 s |
| DP | 2.41 s |
| glue | 0.11 s |
| neighbour lists | 2.87 s |
| 2-opt | 6.95 s |
| Or-opt | 9.25 s |

**Glue share G / box tour:** 0.29–0.38 on uniform data and about 0.49 on clustered
data. Between a third and a half of the box tour's length is joining edges.

## 7. The three questions

**(1) After cleanup, does the box size t still matter? No, not in any practical
way.**
- **Before cleanup:** larger boxes help a little, about 1–1.5 points, e.g. uniform 100k: +28.97% at t = 6 vs +27.46% at t = 12.
- **After 2-opt:** the difference shrinks to about 1 point.
- **After Or-opt:** it shrinks to 0.1–0.3 points (uniform 100k: 8.13 / 8.08 / 8.02 / 7.95%), which is about the size of the seed-to-seed spread. On clustered data there is no consistent order at all.
- **Time:** t = 12 costs 1.4–1.8× the time of t = 10 at 100k–1M.

One reason the effect is small: t is only a *cap*. Because a quadtree cell
splits into four, the average box holds about **4–5 cities** whatever t is (with
t = 10: 4.0–5.4 cities per box).

**(2) Does the full method beat greedy + 2-opt + Or-opt? No.** It ends **3.3–4.2%
longer** on uniform data and **2.9–7.8% longer** on clustered data, at every size,
with the same cleanup code. It isn't faster either (1.4× the time at 100k, 1.1–1.8×
at 1M).
- **The box start itself is better than the plain Hilbert order** (about 10 points
  better before cleanup), so the exact per-box solving does improve the
  space-filling tour.
- **But it starts much worse than greedy** (+28% vs +14% on uniform 100k), and
  neighbour-list 2-opt and Or-opt can't close that gap. The long joining edges
  between boxes need fixes that connect cities outside each other's 16-nearest
  lists, so those moves are never tried.
- This matches the earlier finding in this bench that 2-opt from a space-filling
  start ends far worse than 2-opt from greedy.

**(3) Do uniform and clustered data behave differently? Yes.** Clustered data
hurts the box method more:
- **Before cleanup:** the box tour is about 46% over the Held-Karp bound, against about 28% on uniform data.
- **Glue share:** about 0.49, against about 0.3 on uniform data.
- **After full cleanup:** about 7% behind the greedy baseline at 10k–1M, against about 4% on uniform data.

The Hilbert order runs through empty space between clusters, so many joining
edges are long, inter-cluster seams. Those are exactly the edges that
nearest-neighbour cleanup can't reach, because every city's 16 nearest neighbours
are inside its own cluster. At 1k clustered the difference to the baseline (+2.9 ±
4.4%) isn't distinguishable from noise.

## 8. Known limitations and things not done

- **K:** the cleanup uses K = 16 neighbours (the bench's value), not 8–10. `--k` exists if you want to test this.
- **Cleanup comparison:** the box variants and `greedy+2opt+oropt` repeat cleanup rounds until clean. The existing `twoopt_or` runs one work-queue pass. The difference in tour length is ≤ 0.1%, but the repeated rounds cost time at 1M (20 s vs 11 s).
- **Missing references:** no Held-Karp bound exists for C100k.2 and C1M.*, so clustered gaps cover fewer seeds there. The paired comparison covers all of them.
- **Stage times:** these use `clock()` inside the solver. The run totals are user CPU time measured outside it.
- **Random shift:** the optional random-shift experiment (Arora-style shifted grids) was not run.
- **Bounds:** nothing here is presented as a bound. G is reported only as a measured share. Prior work already covers boxes plus exact solving (Karp 1977, Arora, Huang et al., DCAP). What this adds is the Hilbert-ordered glue, the cleanup, and this comparison.

## 9. Reproduce

```powershell
$env:PATH = "C:\msys64\ucrt64\bin;C:\msys64\usr\bin;" + $env:PATH
cd solvers; make; cd ..; .\bin\boxtest2d.exe             # build + 33 tests
py bench\box_instances.py                           # data + clumpiness check
py bench\run.py --config bench\box_config.json      # ~30 min, resumable
py bench\box_report.py                              # results/2d/box/summary.md
```
