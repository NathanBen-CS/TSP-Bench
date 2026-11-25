# Box-TSP in 3D

The 2D box study ([box_tsp.md](box_tsp.md), [box_tsp_2d.md](box_tsp_2d.md))
extended to 3D: the same curves → boxes → Held-Karp loops → merge → cleanup,
compared with **greedy + 2-opt + Or-opt**. It is a Karp-style heuristic with no
worst-case guarantee, and it is not a PTAS. (Arora's scheme in R^d has a running
time with an exponent of about (c·√d)^(d−1), so it is not a practical route in 3D
either.)

**Short answer: no.** No box setting beats **greedy TO Rn** on average on any
3D map group with 1,000 cities or more: uniform, clustered, tight clustered or
stars, from 1k to 1M cities, and gaia1 (2.08M stars). The margin is about as in 2D
on uniform maps (2D: about 3 points above the bound) and smaller on clustered maps
(2D: about 4.5), but it never reaches 0:
- **Best box setting after full cleanup**, % longer than greedy TO Rn (section 8):
  - **Uniform:** +2.0 to +3.2.
  - **Stars:** +2.4 to +3.4.
  - **Clustered and tight clustered:** +1.6 to +2.5 at 10k–1M, +0.5 to +0.8 at 1k.
- **Clustered helps a little, not enough:** the gap is smallest on tight
  clustered maps with merge M5 (+1.9 to +2.4% at 10k–1M), as in 2D.
- **Single wins only at 1k:** of 900 single runs at 1k (all settings, 3 seeds),
  26 beat greedy TO Rn: 7 on clustered and 19 on tight clustered maps. That is
  within seed noise there (sd 0.6–2.6 points); no setting's *mean* is below 0.
  At 10k and above, none of 4,801 single runs beats it.
- **star100 is the only exception:** a single map of 100 stars, where 28 of 90
  box settings beat greedy TO Rn (best −3.39%), because greedy TO Rn is
  unusually weak there (3.5% above optimal). 100 stars is too small to say anything.
- **Why:** in 3D a box holds on average only about 3 cities at t = 10 (2D: 4–5),
  82% of each city's 16 nearest neighbours lie outside its box (2D: 72%), and
  57% of the edges of a good tour cross a box border (2D: 38%) (section 7).
- **t and K:** t barely matters (section 11); a larger K helps M1–M4 on tight
  clustered maps a lot, but helps greedy TO Rn as well, and M5's gap does not close (section 12).
- **Cost:** box runs take about 1.5–2.7× the user time of greedy TO Rn at
  100k (section 14).

Code names are the 2D ones: a bare name = that order only; B = boxes (exact loops
+ merge); T = 2-opt; O = Or-opt; **R1** = one round of cleanup; **Rn** = rounds
until nothing improves. **greedy TO Rn** is the main baseline.

Raw runs: `results/3d/runs.csv` (6,360 timed runs, all OK) and
`results/3d/runs_untimed.csv` (24 runs, lengths only). Full tables:
[`results/3d/summary.md`](../results/3d/summary.md). Files:
[`results/3d/README.md`](../results/3d/README.md).

---

## 1. What was built

**One set of shared code, two dimensions.** `solvers/common/` holds everything
that only uses `distance()` and the neighbour lists (tour, local search, greedy
start, box split, Held-Karp, merges M1–M5). It is compiled **once** and linked
into both the 2D and the 3D programs. Only three things differ, and the Makefile
picks the implementation behind the same header:

| | 2D | 3D |
|---|---|---|
| reader | `codes/Instance.c` (unchanged) | `solvers/d3/instance3d.c` |
| neighbour lists (`neighbors.h`) | uniform grid | kd-tree, `solvers/d3/neighbors3d.c` |
| curves (`curves.h`) | `solvers/d2/curves2d.c` | `solvers/d3/curves3d.c` |

**No 2D result changed.** Before anything was moved, 244 tours and all 120
`boxtest` lines were recorded from a fresh build. `tests/regress_2d.py`
compares every tour byte for byte, and it passes.

**3D curves**

| Curve | Base | Digits | Closed | Note |
|---|---|---|---|---|
| hilbert3 | 8 | 20 | no | Skilling's transpose algorithm (2004) |
| morton3 | 8 | 20 | no | negative control: 50% of steps go to a face neighbour |
| peano3 | 27 | 13 | no | 3 × 3 × 3 split, Peano's digit-parity rule |
| snake3 | 8 | 20 | no | rows, layers, all alternating. The 3D version of the 2D "strip", at full grid resolution, so a layer holds about one city and the order is close to "sort by z" |
| moore3 (Stage F) | 8 | 20 | yes | eight Hilbert3 curves in a closed cycle of octants (section 10) |
| sierpinski3 (Stage F) | 2 | 60 | yes | tetrahedra cut in half (section 10) |

The box split needed no change for base 8 and 27: the powers (8²⁰ = 2⁶⁰,
27¹³ = 3³⁹ ≈ 4.05 × 10¹⁸) fit in `unsigned long long`.

## 2. Conventions kept and changed

**Kept from 2D:** K = 16 neighbours (`--k`), greedy candidates 10 in round 1 and 8
in the join rounds, Or-opt chains of 1–3, work-queue local search, reversal cap
100,000, rounded integer distances, deterministic runs, `CFLAGS = -O2` everywhere,
user CPU time of the whole process pinned to core 0, every tour checked by the
DIMACS `length` program (linked with the 3D reader: `bin/length3d.exe`).

**Changed, and why:**
- **Neighbour search: kd-tree instead of grid.** A uniform grid is hopeless on the
  stars (hyg119614 has 109,399 of its 119,614 stars within 1% of the box width).
  Lists have the same meaning: true distance, nearest first, deterministic ties.
- **greedy.** The DIMACS `greedy` program is 2D only. "greedy" here is the
  multi-fragment start that `twoopt`, `twoopt_or` and `boxtsp` use in both
  dimensions (`boxtsp3d --start greedy --cleanup none`).
- **R1 for box starts** (`boxtsp3d --rounds 1`) is new. In 2D only greedy had it.
- **Seeds:** 3 per size at 1k to 100k and 1 at 1M, as the task set (2D had 10 at 1k
  and 3 at 1M).
- **Stage times** use `clock()`, which on Windows is wall time, as in 2D.

## 3. Instances

**Uniform (E3d):** `codes/portgen3d.c`, integer coordinates uniform in the
1,000,000 cube. Names `E3d1k.0` and so on, with the DIMACS seeds (.0 is N, .k is
10·N + k).

**Clustered (C3d):** `codes/portcgen3d.c`, the direct analogue of `portcgen.c`.
N/100 centres uniform in the cube; each point is a random centre plus a
truncated normal deviate per axis. As in `portcgen.c`, points are **not clipped** to
the cube. The standard deviation is **10⁶ / N^(1/3)** instead of 10⁶ / N^(1/2): in
2D that is the side of a square holding one city on average, in 3D the side of
such a cube.

**C3d is much less clustered than the 2D C maps.** The scaling keeps the spread in
units of city spacing, but a 3D Gaussian cloud fills far more of its space.
Measured by `bench/check_clumpiness3d.py` (mean nearest-neighbour distance,
clustered / uniform; 1 = not clustered at all):

| size | 2D C / E | C3d / E3d | C3d cities outside the cube |
|---|---|---|---|
| 1k | 0.41 | 0.87 | 32% |
| 10k | 0.42 | 0.79 | 12% |
| 100k | 0.43 | 0.78 | 5% |
| 1M | 0.43 | 0.78 | 2% |

At 1k, where the spread is 10% of the cube's side and a third of the points lie
outside it, C3d is close to a uniform cloud.

**Tight clustered (CT3d), added.** Because of that, a second clustered set uses
half the standard deviation (`portcgen3d N seed 0.5`). On 10k it has a
nearest-neighbour ratio of 0.44, the same as the 2D C maps. It was not in the
original plan; C3d is kept exactly as specified.

**Stars** (Waterloo, `bench/get_star.py`): star100, star1k, star10k, kj37859,
hyg109399, hyg119614, star250k and gaia1 (2,079,471 stars). **star10m is not
used:** its data page offers the file only by email request
([star10m_data.html](https://math.uwaterloo.ca/tsp/star/star10m_data.html)). Mean
nearest-neighbour distance against an even spread over the bounding box (an
indication only, a box is not the stars' shape): 0.77–0.81 for star100 to star250k,
0.61 for kj37859, 0.45 for hyg109399, 0.19 for gaia1, 0.07 for hyg119614.

## 4. Validation

**Waterloo tours** (`tests/validate_star.py`): the published tour files, measured
with our reader and distance function. All eight lengths are reproduced exactly
(units of 1/10 parsec, so the published parsecs × 10):

| map | measured | published | status | source |
|---|---|---|---|---|
| star100 | 1,795 | 179.5 pc | optimal | [tour page](https://math.uwaterloo.ca/tsp/star/star100_tour.html) |
| star1k | 22,227 | 2,222.7 pc | optimal | [tour page](https://math.uwaterloo.ca/tsp/star/star1k_tour.html) |
| star10k | 276,750 | 27,675.0 pc | optimal | [tour page](https://math.uwaterloo.ca/tsp/star/star10k_tour.html) |
| kj37859 | 28,235,453 | 2,823,545.3 pc | optimal | [page](https://math.uwaterloo.ca/tsp/star/kj37859.html) |
| hyg109399 | 13,750,874 | 1,375,087.4 pc | optimal | [tour page](https://math.uwaterloo.ca/tsp/star/hyg_tour.html) |
| hyg119614 | 253,077,734 | 25,307,773.4 pc | optimal | [page](https://math.uwaterloo.ca/tsp/star/hyg119614.html) |
| star250k | 9,794,191 | 979,419.1 pc | bounds 979,414.0 (LP), 979,414.9 (branch and bound) | [page](https://math.uwaterloo.ca/tsp/star/star250k.html) |
| gaia1 | 288,843,524 | 28,884,352.4 pc | bound 28,884,139.0 | [page](https://math.uwaterloo.ca/tsp/star/gaia1.html) |

**Our own Held-Karp bound** (`solvers/common/hkbound.c`, `refs/hk_own.csv`). There
are no published bounds for E3d and C3d, so we compute a 1-tree bound: subgradient
ascent on a sparse graph, then the final numbers evaluated on the **complete**
graph, which makes the printed value a proven lower bound. It can only be too low.
Checked first on 2D maps with a published Held-Karp bound, and on stars with a
known optimum:

| map | ours vs published |
|---|---|
| E1k.0, E1k.1 | −0.004%, −0.008% vs HK |
| C1k.0, C1k.1 | −0.003%, −0.004% vs HK |
| E10k.0, C10k.0 | −0.003%, −0.090% vs HK |
| E100k.0, C100k.0 | −0.015%, −0.274% vs HK |
| star1k, star10k, kj37859 | −0.27%, −0.21%, −0.19% vs the **optimal tour** |

It is never above the published value. On uniform maps it is essentially the
Held-Karp bound; on clustered maps it can be up to about 0.3% too low at 100k, so
gaps over it look that much too big. **Every "own HK" number below is our own
bound.**

## 5. Tests

- `tests/unit3d.c` (`bin/unit3d.exe`): the reader (real and negative coordinates,
  any city order, no EOF line); neighbour lists against brute force on 300 cities
  (uniform, clustered, duplicates, one plane, dense core + far points; K = 1, 8, 16, 32;
  and a subset); every cube curve is a bijection on small grids, steps to a face
  neighbour (Morton: 50%, noted), and its cells at every level are contiguous key
  ranges; Moore3 closes; keys rise along the sorted order and fit the digits.
- `tests/boxtest3d.c` (`bin/boxtest3d.exe`): boxes are whole curve cells for every
  curve, kind of data and t; depth-capped boxes keep curve order; kd boxes
  hold exactly t; Held-Karp equals brute force on 400 sets of 1–8 points; two
  unit cubes two apart (M2–M5 reach the optimum merge 18, M1 gives 19, because
  its first cut is fixed, as in 2D); every curve × merge × data
  kind gives one valid tour whose added cost equals the recount; tiny inputs
  (2–30 cities); M2 equals the brute-force best snake; greedy and cleanup leave no
  improving 2-opt move on 300 cities.
- `tests/regress_2d.py`, `tests/validate_star.py`, `tests/check_links.py`.

All pass.

## 6. Baseline (Stage C)

Before any box run: greedy alone and with every cleanup level, K = 16, on every map
(190 runs). **On the stars**, % above the published optimal or best tour (Waterloo,
section 4); one map each, so no spread:

| map | cities | greedy | greedy T R1 | greedy TO R1 | greedy T Rn | **greedy TO Rn** | TO Rn user s |
|---|---|---|---|---|---|---|---|
| star100 | 100 | +12.87 | +4.46 | +4.01 | +3.57 | **+3.51** | 0.00 |
| star1k | 1,000 | +7.87 | +2.99 | +2.33 | +2.67 | **+2.29** | 0.02 |
| star10k | 10,000 | +7.02 | +2.89 | +2.49 | +2.63 | **+2.32** | 0.12 |
| kj37859 | 37,859 | +7.05 | +2.85 | +2.52 | +2.55 | **+2.21** | 0.66 |
| hyg109399 | 109,399 | +6.57 | +2.79 | +2.40 | +2.47 | **+2.16** | 2.81 |
| hyg119614 | 119,614 | +12.58 | +4.70 | +4.13 | +4.56 | **+3.89** | 2.94 |
| star250k | 250,000 | +6.45 | +3.03 | +2.64 | +2.50 | **+2.18** | 6.95 |
| gaia1 | 2,079,471 | +6.49 | +4.56 | +4.27 | +2.87 | **+2.62** | 169.69 |

**On the synthetic maps**, % above **our own** Held-Karp bound (section 4), mean ± sd over 3 seeds:

| maps | size | greedy | greedy T R1 | greedy TO R1 | greedy T Rn | **greedy TO Rn** |
|---|---|---|---|---|---|---|
| uniform | 1k | +8.80 ± 0.36 | +3.05 ± 0.17 | +2.42 ± 0.19 | +2.77 ± 0.10 | **+2.40 ± 0.15** |
| uniform | 10k | +7.45 ± 0.14 | +3.08 ± 0.09 | +2.78 ± 0.14 | +2.83 ± 0.07 | **+2.49 ± 0.09** |
| uniform | 100k | +6.93 ± 0.09 | +3.03 ± 0.04 | +2.62 ± 0.02 | +2.72 ± 0.02 | **+2.39 ± 0.02** |
| clustered | 1k | +11.03 ± 0.71 | +4.59 ± 1.12 | +3.86 ± 0.82 | +4.25 ± 0.85 | **+3.65 ± 0.83** |
| clustered | 10k | +9.40 ± 0.41 | +4.51 ± 0.31 | +3.78 ± 0.29 | +4.08 ± 0.32 | **+3.56 ± 0.41** |
| clustered | 100k | +8.46 ± 0.12 | +4.11 ± 0.04 | +3.53 ± 0.07 | +3.69 ± 0.02 | **+3.15 ± 0.02** |
| tight clustered | 1k | +12.68 ± 0.98 | +5.72 ± 0.81 | +5.37 ± 0.97 | +5.47 ± 0.84 | **+5.02 ± 0.85** |
| tight clustered | 10k | +12.15 ± 0.65 | +5.84 ± 0.70 | +5.10 ± 0.69 | +5.52 ± 0.80 | **+4.81 ± 0.85** |
| tight clustered | 100k | +11.92 ± 0.07 | +5.96 ± 0.10 | +5.13 ± 0.02 | +5.66 ± 0.07 | **+4.88 ± 0.05** |

**greedy TO Rn is a strong baseline in 3D**, as the task expected: 2.2–2.6% above the
optimal or best tour on 6 of the 8 star maps (3.5% on star100, 3.9% on hyg119614, where
109,399 of the 119,614 stars sit in a tiny core), and 2.4–2.5% above our own bound on
uniform maps. On clustered maps the gap to our bound is larger (C3d 3.2–3.7%, CT3d
4.8–5.0%). We cannot say how much of that is the tour and how much is our bound being too low.
On the 2D C maps, which are about as clustered as CT3d, our bound was up to 0.27% below the
published one (section 4); for CT3d there is no published bound to check against.
Every box result below is measured against greedy TO Rn on the same map, not against the bound.

**R1 against Rn at 1M.** One round of cleanup is close to rounds-until-clean up to 250k
(greedy TO R1 0.0–0.5 points behind TO Rn) but not on uniform and clustered 1M and on gaia1
(1.4–1.6 points behind; on gaia1 +4.27% against +2.62% above the best tour). Tight
clustered 1M is again close (0.22). Section 9 has the box runs.

## 7. Why boxes do less in 3D (diagnostics)

`bench/box3_diag.py` (`results/3d/diagnostics.csv`) builds the boxes exactly as the
solver does and measures them against the **greedy TO Rn** tour of the same map,
with K = 16. The 2D rows use the 2D Hilbert curve on the 2D maps. One map each.

| map | avg cities per box (t = 10) | of a city's 16 neighbours, outside its box | city whose nearest neighbour is outside its box | tour edges crossing a box border | their share of the tour's length |
|---|---|---|---|---|---|
| 2D E10k, hilbert | 4.48 | 73% | 28% | 39% | 48% |
| 2D E100k, hilbert | 5.38 | 69% | 26% | 36% | 46% |
| 2D C100k, hilbert | 4.32 | 73% | 28% | 39% | 57% |
| E3d10k, hilbert3 | 2.70 | 85% | 48% | 59% | 66% |
| E3d100k, hilbert3 | 3.19 | 82% | 47% | 57% | 63% |
| E3d1M, hilbert3 | 3.86 | 79% | 45% | 54% | 61% |
| C3d100k, hilbert3 | 3.10 | 82% | 47% | 57% | 64% |
| CT3d100k, hilbert3 | 3.13 | 81% | 46% | 56% | 65% |
| star250k, hilbert3 | 2.86 | 83% | 49% | 59% | 65% |
| gaia1, hilbert3 | 3.28 | 81% | 47% | 57% | 63% |

No box hit the depth cap on any map.

- **Boxes are small.** A cube cell splits into 8, so a cell that holds 11 cities becomes
  8 cells holding about 1.4 each. t = 10 gives 2.7–3.9 cities per box, and t = 12
  hardly more (E3d100k: 3.19 at t = 10, 3.21 at t = 12). A box of 3 cities has only one possible
  loop, so Held-Karp has nothing to optimise inside it.
- **The boxes cut the good tour in many places.** More than half of the edges of
  greedy TO Rn cross a box border, and these carry about two thirds of its
  length. Whatever the merge does with them, the box tour has to be rebuilt
  between boxes by the cleanup, which is what greedy + cleanup already does well.
- **Other curves:** peano3 splits a cell into 27, so its box size jumps with N (E3d:
  1.44 at 10k, 4.52 at 100k). sierpinski3 halves a tetrahedron at every step, so
  its boxes come close to t (6.9 at t = 10, like the 2D Sierpinski curve). But even
  its boxes are crossed by 54% of the tour's edges (E3d100k), against 36% for 2D Hilbert.

## 8. Main result: curve × merge

Every curve (hilbert3, moore3, peano3, sierpinski3, morton3, snake3) with every merge
(M1–M5) at t = 10, K = 16, full cleanup (BTO Rn): 30 settings on every map up to 250k
(at 1M and gaia1: 18, section 13). % longer than greedy TO Rn on the same map, mean ± sd
over 3 seeds (stars and 1M: one map, no spread). "Single runs < 0" counts every
setting × seed that came out shorter than greedy TO Rn.

| maps | best setting | best | hilbert3 M1 | hilbert3 M5 | single runs < 0 |
|---|---|---|---|---|---|
| uniform 1k | peano3 M2 | +2.00 ± 0.15 | +2.27 ± 0.20 | +3.06 ± 0.36 | 0 of 90 |
| uniform 10k | snake3 M5 | +2.59 ± 0.11 | +2.91 ± 0.11 | +2.94 ± 0.16 | 0 of 90 |
| uniform 100k | snake3 M2 | +2.77 ± 0.02 | +3.31 ± 0.08 | +3.20 ± 0.02 | 0 of 90 |
| uniform 1M | sierpinski3 M5 | +3.17 | +3.41 | +3.29 | 0 of 18 |
| clustered 1k | sierpinski3 M4 | +0.84 ± 1.01 | +2.00 ± 0.59 | +1.78 ± 1.01 | 4 of 90 |
| clustered 10k | moore3 M5 | +1.84 ± 0.61 | +2.34 ± 0.61 | +1.88 ± 0.67 | 0 of 90 |
| clustered 100k | morton3 M5 | +2.30 ± 0.02 | +2.87 ± 0.11 | +2.31 ± 0.08 | 0 of 90 |
| clustered 1M | peano3 M5 | +2.45 | +2.91 | +2.71 | 0 of 18 |
| tight clustered 1k | morton3 M4 | +0.66 ± 1.13 | +0.88 ± 0.81 | +2.46 ± 1.92 | 8 of 90 |
| tight clustered 10k | moore3 M5 | +2.22 ± 0.97 | +3.72 ± 1.35 | +2.35 ± 0.60 | 0 of 90 |
| tight clustered 100k | peano3 M5 | +1.94 ± 0.08 | +4.42 ± 0.23 | +1.99 ± 0.11 | 0 of 90 |
| tight clustered 1M | moore3 M5 | +2.01 | +5.32 | +2.04 | 0 of 18 |
| star100 | moore3 M4 | **−3.39** | −1.24 | +0.22 | 12 of 30 |
| star1k | snake3 M2 | +2.42 | +2.81 | +3.01 | 0 of 30 |
| star10k | snake3 M2 | +2.75 | +3.19 | +2.99 | 0 of 30 |
| kj37859 | snake3 M2 | +2.65 | +3.12 | +3.19 | 0 of 30 |
| hyg109399 | snake3 M5 | +2.80 | +3.37 | +3.29 | 0 of 30 |
| hyg119614 | peano3 M2 | +2.93 | +3.14 | +4.35 | 0 of 30 |
| star250k | snake3 M5 | +2.94 | +3.42 | +3.33 | 0 of 30 |
| gaia1 | sierpinski3 M5 | +3.40 | +3.59 | +3.73 | 0 of 18 |

**Counting every setting we ran** (other t, other K, kd boxes, R1 and T cleanups:
up to 305 settings per map group, `results/3d/summary.md` section 3), the picture is
the same:
- **No mean below 0** except star100.
- **The best means at 10k and above:** clustered 10k hilbert3 M5 t = 7, +1.63 ± 0.39;
  clustered 100k sierpinski3 M5 K = 32, +2.22 ± 0.08; tight clustered 100k hilbert3 M5 K = 32, +1.79 ± 0.12.
- **Single runs below 0, only at 1k:** 0 of 300 on uniform, 7 of 300 on clustered and
  19 of 300 on tight clustered maps. The best mean there is kd M2 t = 12, +0.47 ± 0.70.
- **None at 10k and above:** 0 of 4,801 single runs.

How to read it:
- **The 1k wins are noise.** On clustered maps at 1k the seed sd is 0.6 to 2.6 points, larger than
  the mean gap of the best settings (+0.5 to +0.9). Some seeds go below 0, and no
  setting is below 0 on average. With 3 seeds per setting and 30–100 settings, a few
  runs below 0 are to be expected when the mean is that close.
- **star100** is one map of 100 stars. greedy TO Rn is 3.51% above optimal there, weak by
  its own standard (2.2–2.6% on the other stars), so the boxes get 12 of 30 wins. It shows
  that boxes can win on a tiny map, nothing more.
- **Clustered vs uniform:** the gap is smaller on clustered maps (best +1.6 to +2.5 at
  10k–1M, against +2.6 to +3.2 on uniform), as in 2D, and **M5 is the merge that does it**:
  on tight clustered 100k hilbert3 M5 is +1.99 against +4.42 for M1. On uniform maps and stars the
  merge hardly matters after cleanup.
- **snake3** (sort by z, layer by layer) is the best curve on uniform maps and most stars
  after cleanup, as "strip" was on uniform maps in 2D, and the worst on tight clustered maps
  (+18 to +22% at 100k). Its start is 10–50 times as long as greedy TO Rn and takes 9–16×
  the time (section 14).
- **Before cleanup** M4 gives the shortest box tours (+10 to +14% above greedy TO Rn at 100k),
  M1 and M2 the longest (+27 to +64%). Chart `box3_curve_merge.png` shows both levels.

Charts: `box3_curve_merge.png` (every curve and merge, 100k and stars),
`box3_uniform_vs_clustered.png` (by map type), `box3_overview.png` (against the
bound and the best tours), `box3_scale.png` (1k to 1M).

## 9. Less cleanup: R1 and T

Before cleanup, box starts (other than snake3) are 10–65% longer than greedy TO Rn,
while the greedy start is only 4–7% longer (uniform 100k: +4.4%, tight clustered 100k:
+6.7%). So box starts depend more on the cleanup. With
one round (R1) or 2-opt alone (T), the gap gets bigger, not smaller. Best box setting
of all curves and merges at t = 10, % longer than greedy TO Rn:

| maps | greedy TO R1 | best box TO R1 | greedy T R1 | best box T R1 | best box TO Rn |
|---|---|---|---|---|---|
| uniform 10k | +0.28 ± 0.05 | +3.07 ± 0.19 | +0.58 ± 0.01 | +5.22 ± 0.04 | +2.59 ± 0.11 |
| uniform 100k | +0.23 ± 0.03 | +3.23 ± 0.09 | +0.62 ± 0.02 | +5.50 ± 0.09 | +2.77 ± 0.02 |
| uniform 1M | +1.49 | +5.35 | +1.82 | +7.95 | +3.17 |
| clustered 10k | +0.21 ± 0.20 | +2.21 ± 0.49 | +0.91 ± 0.18 | +5.05 ± 0.66 | +1.84 ± 0.61 |
| clustered 100k | +0.37 ± 0.05 | +2.76 ± 0.03 | +0.93 ± 0.02 | +4.92 ± 0.06 | +2.30 ± 0.02 |
| clustered 1M | +1.38 | +3.47 | +1.84 | +5.94 | +2.45 |
| tight clustered 10k | +0.28 ± 0.16 | +2.57 ± 0.68 | +0.98 ± 0.18 | +4.86 ± 0.88 | +2.22 ± 0.97 |
| tight clustered 100k | +0.24 ± 0.04 | +2.26 ± 0.10 | +1.03 ± 0.08 | +4.55 ± 0.11 | +1.94 ± 0.08 |
| tight clustered 1M | +0.22 | +2.38 | +0.97 | +4.66 | +2.01 |
| star10k | +0.17 | +3.24 | +0.56 | +5.77 | +2.75 |
| star250k | +0.45 | +4.15 | +0.83 | +6.32 | +2.94 |
| gaia1 | +1.61 | +4.96 | +1.89 | +7.32 | +3.40 |

At 1M only hilbert3 was run with R1, so the 1M and gaia1 R1 columns are hilbert3's best merge.
At 1M and on gaia1 **one round of greedy cleanup leaves much more undone** than at
100k (uniform: +1.49 against +0.23). We did not isolate why; tight clustered 1M does
not show it (+0.22).

## 10. Stage F: Moore3, Sierpinski3 and kd-median boxes

**Moore3.** Haverkort's inventory of 3D Hilbert curves
([arXiv:1109.2323](https://arxiv.org/abs/1109.2323)) shows there are many 3D
Hilbert-type curves and closed variants, with no single standard "3D Moore". We
built the direct analogue of the 2D construction (Moore = four Hilbert curves in a
closed cycle): eight hilbert3 curves, one per octant, visited in the Gray-code cycle
x y x z x y x z, each turned so that it starts where the last one ended. That makes
it closed and keeps hilbert3's locality inside each octant. It is our construction;
we did not check it against a specific curve of that paper.

**Sierpinski3.** The 2D Sierpinski curve halves triangles. The task pointed to
[arXiv:1711.04473](https://arxiv.org/abs/1711.04473) for a Sierpinski-type 3D curve; we
did not reproduce that paper's curve but built one on the same idea. We cut the cube into
the 6 Kuhn tetrahedra and halve them by Maubach's bisection. For
the walk order we searched a 30-state rule table by hill climbing so that
consecutive tetrahedra always share at least one corner (81% share a face); the
curve is closed. Its plain order is worse than hilbert3's (+62% against +46% above
greedy TO Rn on uniform 100k), and it is our own table, not a published curve.

After full cleanup, t = 10 (section 8 has all merges):

| maps | hilbert3 M1 | hilbert3 M5 | moore3 M1 | moore3 M5 | sierpinski3 M1 | sierpinski3 M5 |
|---|---|---|---|---|---|---|
| uniform 100k | +3.31 ± 0.08 | +3.20 ± 0.02 | +3.32 ± 0.04 | +3.12 ± 0.02 | +3.07 ± 0.04 | +3.01 ± 0.04 |
| uniform 1M | +3.41 | +3.29 | +3.41 | +3.27 | – | +3.17 |
| clustered 100k | +2.87 ± 0.11 | +2.31 ± 0.08 | +2.75 ± 0.06 | +2.31 ± 0.04 | +2.64 ± 0.04 | +2.31 ± 0.08 |
| clustered 1M | +2.91 | +2.71 | +2.83 | +2.73 | – | +2.50 |
| tight clustered 100k | +4.42 ± 0.23 | +1.99 ± 0.11 | +4.23 ± 0.15 | +2.01 ± 0.09 | +4.59 ± 0.11 | +2.65 ± 0.12 |
| tight clustered 1M | +5.32 | +2.04 | +5.35 | +2.01 | – | +3.04 |
| star10k | +3.19 | +2.99 | +2.90 | +2.91 | +3.03 | +3.13 |
| star250k | +3.42 | +3.33 | +3.45 | +3.33 | +3.31 | +3.29 |
| gaia1 | +3.59 | +3.73 | +3.69 | +3.75 | – | +3.40 |

**sierpinski3 with every merge** (t = 10, BTO Rn, % longer than greedy TO Rn; at 1M and on
gaia1 only M5 was run). "Rank" is the place of its best merge among the curves run on those maps.

| maps | M1 | M2 | M3 | M4 | M5 | rank | best other curve |
|---|---|---|---|---|---|---|---|
| uniform 1k | +2.43 ± 0.77 | +2.55 ± 0.23 | +2.35 ± 0.30 | +2.30 ± 0.20 | +2.37 ± 0.19 | 5 of 6 | peano3 M2 +2.00 |
| uniform 10k | +2.85 ± 0.25 | +2.90 ± 0.19 | +2.72 ± 0.20 | +2.93 ± 0.18 | +2.80 ± 0.14 | 2 of 6 | snake3 M5 +2.59 |
| uniform 100k | +3.07 ± 0.04 | +3.07 ± 0.07 | +3.12 ± 0.05 | +3.12 ± 0.00 | +3.01 ± 0.04 | 2 of 6 | snake3 M2 +2.77 |
| uniform 1M | – | – | – | – | +3.17 | 1 of 5 | morton3 M1 +3.20 |
| clustered 1k | +1.82 ± 0.91 | +1.54 ± 0.71 | +1.59 ± 0.65 | +0.84 ± 1.01 | +1.08 ± 1.23 | 1 of 6 | peano3 M4 +0.88 |
| clustered 10k | +2.25 ± 0.58 | +2.47 ± 0.52 | +2.31 ± 0.48 | +2.09 ± 0.61 | +2.09 ± 0.45 | 6 of 6 | moore3 M5 +1.84 |
| clustered 100k | +2.64 ± 0.04 | +2.64 ± 0.04 | +2.57 ± 0.08 | +2.44 ± 0.09 | +2.31 ± 0.08 | 4 of 6 | morton3 M5 +2.30 |
| clustered 1M | – | – | – | – | +2.50 | 2 of 5 | peano3 M5 +2.45 |
| tight clustered 1k | +1.46 ± 2.67 | +1.34 ± 0.68 | +2.32 ± 0.27 | +2.66 ± 2.46 | +1.04 ± 2.02 | 4 of 6 | morton3 M4 +0.66 |
| tight clustered 10k | +4.48 ± 1.05 | +4.33 ± 0.20 | +5.44 ± 1.63 | +3.29 ± 0.39 | +2.59 ± 0.37 | 5 of 6 | moore3 M5 +2.22 |
| tight clustered 100k | +4.59 ± 0.11 | +4.93 ± 0.19 | +4.42 ± 0.33 | +3.30 ± 0.16 | +2.65 ± 0.12 | 5 of 6 | peano3 M5 +1.94 |
| tight clustered 1M | – | – | – | – | +3.04 | 5 of 5 | moore3 M5 +2.01 |
| star100 | −0.59 | +0.54 | +1.18 | −0.70 | −0.48 | 3 of 6 | moore3 M4 −3.39 |
| star1k | +3.10 | +3.28 | +2.82 | +3.33 | +3.42 | 5 of 6 | snake3 M2 +2.42 |
| star10k | +3.03 | +3.07 | +2.93 | +3.12 | +3.13 | 3 of 6 | snake3 M2 +2.75 |
| kj37859 | +2.99 | +3.10 | +2.91 | +3.14 | +3.07 | 2 of 6 | snake3 M2 +2.65 |
| hyg109399 | +3.13 | +3.22 | +3.15 | +3.20 | +3.15 | 2 of 6 | snake3 M5 +2.80 |
| hyg119614 | +7.12 | +8.73 | +6.06 | +3.40 | +3.61 | 5 of 6 | peano3 M2 +2.93 |
| star250k | +3.31 | +3.28 | +3.31 | +3.28 | +3.29 | 3 of 6 | snake3 M5 +2.94 |
| gaia1 | – | – | – | – | +3.40 | 1 of 5 | morton3 M1 +3.50 |

Its user time is 1.2–1.8× that of hilbert3 with the same merge (uniform 100k M5: 7.8 s
against 4.9 s), because the cleanup has a longer start to repair.

- **Moore3 equals hilbert3** within seed noise everywhere. Closing the curve
  changes one join out of thousands.
- **Sierpinski3** makes the longest box starts (M5 before cleanup: +31.7% on uniform 100k
  against +20.3% for hilbert3), but after cleanup it is 0.1–0.3 points better on uniform maps
  and gaia1, and 0.6–1.0 worse on tight clustered maps. It never comes near 0.

**kd-median boxes** (`--split kd`, `solvers/d3/kdboxes3d.c`): the cities are cut at
the median of the widest axis until each piece holds t cities (one may hold fewer),
and the pieces are visited in hilbert3 order of their centres. So a box really holds t
cities (curve cells hold about 3). Run with every merge at t = 4, 6, 8, 10, 12 on
the 10k maps, star10k and kj37859; at t = 10 and 12 on the 1k and 100k maps, hyg109399 and
star250k; and M5 t = 12 at 1M and on gaia1 (479 runs).

| maps | kd M5 t=12 | cells M5 t=10 | kd M4 t=10 | cells M4 t=10 | kd M1 t=12 | cells M1 t=10 |
|---|---|---|---|---|---|---|
| uniform 100k | +3.26 ± 0.05 | +3.20 ± 0.02 | +3.30 ± 0.02 | +3.52 ± 0.04 | +3.19 ± 0.03 | +3.31 ± 0.08 |
| clustered 100k | +2.71 ± 0.12 | +2.31 ± 0.08 | +2.71 ± 0.06 | +2.60 ± 0.03 | +2.71 ± 0.07 | +2.87 ± 0.11 |
| tight clustered 100k | +4.05 ± 0.25 | +1.99 ± 0.11 | +3.93 ± 0.15 | +3.41 ± 0.08 | +4.77 ± 0.21 | +4.42 ± 0.23 |
| tight clustered 1M | +6.45 | +2.04 | – | – | – | – |
| star250k | +3.54 | +3.33 | +3.57 | +3.67 | +3.42 | +3.42 |

- **Bigger boxes do not help.** On uniform maps kd boxes are within 0.25 points of curve cells.
- **On tight clustered maps they are much worse with M5** (+4.05 against +1.99 at 100k,
  +6.45 against +2.04 at 1M). A median cut goes through the middle of a cluster, where a
  cube cell more often follows the empty space between clusters.
- **Before cleanup kd boxes are worse** with every merge (for example M4 on uniform 100k:
  +18.9% against +12.1%), so the extra cities per box do not make up for the worse borders.
- The only box setting with the best mean on a whole group is kd at 1k (tight
  clustered 1k, kd M2 t = 12: +0.47 ± 0.70), within noise of the others there.

## 11. Box size t

hilbert3, peano3, morton3 and sierpinski3 with every merge at t = 4 to 12, on the 10k maps and star10k
(1,750 runs). Full table: `results/3d/summary.md`, section 4; charts `box3_t_sweep.png` (hilbert3) and
`box3_t_sweep_sierpinski3.png`.

- **After cleanup t does not matter** beyond noise. On uniform 10k all 60 values lie
  between +2.58 and +3.28. The seed sd is about 0.1–0.2, and no curve shows a trend.
- **On tight clustered 10k** the merge matters, t does not: M5 stays at +2.0 to +2.7 for every t and
  curve, the others at +3.3 to +6.9.
- **Before cleanup** a larger t helps M1 and M2 a little (hilbert3 M1 on clustered 10k:
  +42% at t = 4, +34% at t = 12), and M3–M5 not at all.
- **Why so flat:** for hilbert3 the boxes hardly change with t (section 7: 1.93 cities at
  t = 4, 3.19 at t = 10, 3.21 at t = 12 on E3d100k), because a cell splits into 8.
- **sierpinski3 is the test of that**, because its boxes do grow with t (2.8 cities at t = 4,
  7.0 at t = 10, 8.3 at t = 12 on E3d100k). Its box tours get *longer* before cleanup as t grows
  (M4 on uniform 10k: +14% at t = 4, +23% at t = 12). After cleanup there is still no trend:
  - **Uniform 10k:** M5 +2.58 ± 0.15 at t = 4 (the best mean of any setting there) and +2.80 at t = 10.
  - **Clustered 10k:** M5 between +1.69 and +2.09 for every t.
  - **Tight clustered 10k:** M5 between +2.21 and +2.70.

  So larger boxes do not help either. The differences between t values stay within the seed sd.

## 12. Neighbour-list size K

K = 8, 16, 24, 32 at 10k and 100k (3 seeds each) and on star10k, kj37859, hyg109399,
for greedy TO Rn, and for hilbert3 and sierpinski3 with every merge (819 runs; sierpinski3 added
later, 315 runs, 0 of them shorter than greedy TO Rn at the same K). Charts `box3_k_sweep.png`,
`box3_k_time.png`.

**What K controls** (`--k`, one list per city, nearest first): the candidate moves of
2-opt and Or-opt, the joins tried by merges M3, M4 and M5, and greedy's first round
of candidates (min(10, K)). **It does not change** the boxes, Held-Karp, merges M1 and
M2, or greedy's later join rounds (8 candidates).

% longer than greedy TO Rn **with K = 16**, 100k, mean ± sd over 3 seeds:

| maps | setting | K = 8 | K = 16 | K = 24 | K = 32 |
|---|---|---|---|---|---|
| uniform | greedy TO Rn | +0.06 ± 0.03 | 0 | −0.00 ± 0.02 | −0.03 ± 0.03 |
| uniform | hilbert3 M1 | +3.39 ± 0.07 | +3.31 ± 0.08 | +3.29 ± 0.05 | +3.27 ± 0.03 |
| uniform | hilbert3 M5 | +3.21 ± 0.06 | +3.20 ± 0.02 | +3.18 ± 0.03 | +3.14 ± 0.03 |
| clustered | greedy TO Rn | +0.19 ± 0.05 | 0 | −0.06 ± 0.01 | −0.12 ± 0.07 |
| clustered | hilbert3 M1 | +3.10 ± 0.06 | +2.87 ± 0.11 | +2.63 ± 0.04 | +2.59 ± 0.05 |
| clustered | hilbert3 M5 | +2.47 ± 0.06 | +2.31 ± 0.08 | +2.25 ± 0.05 | +2.27 ± 0.05 |
| tight clustered | greedy TO Rn | +0.30 ± 0.03 | 0 | −0.14 ± 0.05 | −0.20 ± 0.04 |
| tight clustered | hilbert3 M1 | +5.76 ± 0.17 | +4.42 ± 0.23 | +3.69 ± 0.20 | +3.28 ± 0.07 |
| tight clustered | hilbert3 M4 | +4.19 ± 0.08 | +3.41 ± 0.08 | +3.13 ± 0.15 | +2.82 ± 0.11 |
| tight clustered | hilbert3 M5 | +2.37 ± 0.08 | +1.99 ± 0.11 | +1.86 ± 0.07 | +1.79 ± 0.12 |
| uniform | sierpinski3 M5 | +3.12 ± 0.03 | +3.01 ± 0.04 | +3.05 ± 0.04 | +3.03 ± 0.06 |
| clustered | sierpinski3 M1 | +3.13 ± 0.09 | +2.64 ± 0.04 | +2.51 ± 0.04 | +2.47 ± 0.09 |
| clustered | sierpinski3 M5 | +2.59 ± 0.05 | +2.31 ± 0.08 | +2.27 ± 0.10 | +2.22 ± 0.08 |
| tight clustered | sierpinski3 M1 | +6.26 ± 0.39 | +4.59 ± 0.11 | +3.94 ± 0.13 | +3.58 ± 0.09 |
| tight clustered | sierpinski3 M4 | +4.01 ± 0.17 | +3.30 ± 0.16 | +2.94 ± 0.27 | +2.72 ± 0.23 |
| tight clustered | sierpinski3 M5 | +3.30 ± 0.08 | +2.65 ± 0.12 | +2.40 ± 0.12 | +2.34 ± 0.10 |

- **Uniform maps and stars:** K changes little (uniform: within 0.16 points, about the seed sd;
  stars, one map each: within 0.33).
- **Tight clustered maps:** a larger K helps the box runs a lot (M1: −2.5 points from
  K = 8 to 32), because a city's neighbours on the far side of a cluster gap only enter longer
  lists. But greedy TO Rn also gets better with K, so **against greedy TO Rn at the same K**
  M5 is +2.07 at K = 8 and +1.99 at K = 32. M5 stays about 2 points behind.
- **Time** (100k, user seconds): greedy TO Rn 2.2 → 2.9 s from K = 8 to 32 on uniform maps,
  hilbert3 M1 5.6 → 6.3 s, M4 4.0 → 8.1 s (its join search grows with K), M5 4.2 → 6.2 s.
- **sierpinski3 behaves like hilbert3:** flat on uniform maps and stars (within 0.25 points), better with K on
  tight clustered maps (M1 −2.7 points, M5 −1.0 from K = 8 to 32), but greedy TO Rn improves too, so
  against greedy TO Rn at the same K M5 is +3.00 at K = 8 and +2.54 at K = 32 on tight clustered 100k
  (hilbert3 M5: +2.07 and +1.99). On clustered 100k its M5 at K = 32 (+2.22 ± 0.08) is the best mean of
  any setting there, level with hilbert3 M5 (+2.27 ± 0.05) within seed noise; against greedy TO Rn at
  K = 32 it is +2.34. At 10k the same pattern holds with larger seed sd (0.4–0.6 on clustered maps).
  Time at 100k: sierpinski3 M5 7.3 → 9.1 s on uniform maps from K = 8 to 32, 3.2 → 4.6 s on tight clustered.
- **K = 16 is the right default** for the comparison: no setting closes the gap at any K.

## 13. 1 million cities and gaia1

One map each: E3d1M.0, C3d1M.0, CT3d1M.0 and gaia1 (2,079,471 stars). 18 box settings
after full cleanup (hilbert3 with every merge; morton3 and peano3 with every merge;
moore3 M1 and M5; sierpinski3 M5), kd M5 t = 12, and hilbert3 with R1 cleanups. The
morton3 and peano3 M2–M4 runs are lengths only (section 16). % longer than greedy TO Rn,
user seconds in brackets:

| map | greedy TO Rn | greedy TO R1 | greedy T Rn | hilbert3 M1 | hilbert3 M5 | best box setting |
|---|---|---|---|---|---|---|
| uniform 1M | 0 (85 s) | +1.49 (21 s) | +0.32 (51 s) | +3.41 (235 s) | +3.29 (166 s) | sierpinski3 M5 +3.17 (258 s) |
| clustered 1M | 0 (68 s) | +1.38 (18 s) | +0.45 (44 s) | +2.91 (165 s) | +2.71 (60 s) | peano3 M5 +2.45 (50 s) |
| tight clustered 1M | 0 (23 s) | +0.22 (13 s) | +0.70 (16 s) | +5.32 (86 s) | +2.04 (34 s) | moore3 M5 +2.01 (35 s) |
| gaia1 | 0 (170 s) | +1.61 (31 s) | +0.25 (121 s) | +3.59 (433 s) | +3.73 (218 s) | sierpinski3 M5 +3.40 (444 s) |

- **No box setting gets closer at 1M.** The best gap is +2.0 to +3.4, as at 100k.
- **gaia1 against the best known tour** (Waterloo): greedy TO Rn +2.62%, greedy TO R1 +4.27%,
  the best box setting +6.11% (sierpinski3 M5), hilbert3 M5 +6.45%.
- **Memory:** M4 needs much more than the others (625 MB at 1M, 1.25 GB on gaia1, against
  139–330 MB at 1M for M1 and M5, and 170–180 MB for greedy TO Rn).
- No seeds: with one map each, take differences under about 0.2 points as noise (the
  seed sd at 100k is 0.02–0.12).

Chart `box3_1M.png` shows every setting on all four maps.

## 14. Time

User CPU seconds of the whole program (read, neighbour lists, start, cleanup), mean
over seeds, one core, same laptop:

| maps | greedy TO Rn | greedy TO R1 | hilbert3 M1 | hilbert3 M4 | hilbert3 M5 | peano3 M5 | snake3 M2 |
|---|---|---|---|---|---|---|---|
| uniform 100k | 2.5 | 1.8 | 5.9 | 5.2 | 4.9 | 5.6 | 25.5 |
| clustered 100k | 2.4 | 1.6 | 5.4 | 4.3 | 3.5 | 3.1 | 23.4 |
| tight clustered 100k | 1.5 | 1.1 | 3.3 | 3.6 | 2.7 | 2.6 | 18.5 |
| hyg109399 | 2.8 | 2.0 | 7.3 | 5.8 | 6.2 | 5.6 | 32.5 |
| star250k | 7.0 | 4.8 | 19.0 | 13.2 | 14.2 | 10.2 | 111.8 |
| uniform 1M | 85.0 | 20.6 | 235.2 | 159.2 | 166.0 | 94.2 | – |
| clustered 1M | 68.3 | 18.2 | 165.1 | 71.9 | 60.5 | 50.4 | – |
| tight clustered 1M | 23.1 | 12.6 | 86.0 | 50.7 | 34.3 | 32.9 | – |
| gaia1 | 169.7 | 31.0 | 432.9 | 207.4 | 218.4 | 176.6 | – |

- **Building the start is cheap; the cleanup is the cost.** On uniform 100k, the curve
  sort, boxes and Held-Karp of hilbert3 take 0.2 s, M5's merge 0.8 s; greedy itself takes 0.35 s.
  The rest of a box run is the cleanup of a longer start.
- **snake3 is the best curve after cleanup on uniform maps and stars, and by far the
  slowest:** its start is 10–50 times as long as greedy TO Rn, so the cleanup does a huge
  amount of work (25.5 s against 2.5 s at uniform 100k). It was not run at 1M.
- Timing: user time of the whole process, pinned to core 0, laptop plugged in. The laptop
  went into Modern Standby with the screen off for several periods; we checked the stage
  times of the same map in standby and awake (932 and 2,178 runs): the same within 0.01 s.

## 15. What the results show, and what they don't

**They show:**
- On every 3D map group of 1,000 cities or more that we ran (uniform, clustered,
  tight clustered, 8 star maps, 1k to 2.08M cities), **greedy TO Rn is shorter than
  every box setting on average**. At 10k and above the best box setting is 1.6 to 3.4% longer,
  more than the seed noise (sd 0.02–1.0).
- The same holds with less cleanup (R1, T) and at every t and K we tried.
- The merge matters on tight clustered maps (M5 best, as in 2D), the curve matters
  little after cleanup, t not at all.
- The reason we can measure: boxes of about 3 cities that cut more than half of a good tour's edges.

**They don't show:**
- **That boxing can never win in 3D.** We tried one family (curve cells or kd cells,
  t ≤ 12, five merges, 2-opt + Or-opt cleanup). Stronger cleanup (Or-3opt, Lin-Kernighan)
  was not tried on either start; it could narrow or widen the gap.
- **Anything about 100-city maps** beyond star100: one map, where boxes do win.
- **Anything about clustered data in general:** C3d at 1k is hardly clustered (section 3),
  CT3d is one choice of tighter clusters, and the stars are what they are.
- **Exact gaps to optimal on synthetic maps:** those use our own Held-Karp bound,
  which can be up to about 0.3% too low (section 4).
- **Worst-case behaviour:** the method has no guarantee; it is not a PTAS.

## 16. Not done, and limits

- **star10m** not used: the file is only sent on request by email (section 3).
- **snake3 at 1M and on gaia1** not run: at star250k it already took 112 s against 7 s for
  greedy TO Rn. It is the best curve after cleanup on uniform maps and stars, so the 1M
  "best box" in section 13 is the best of the curves that were run there.
- **At 1M and on gaia1** only part of the matrix was timed (hilbert3 with every merge, morton3,
  peano3, moore3 with M1 and M5, sierpinski3 M5, kd M5 t = 12). morton3 and peano3 with M2–M4
  were run for length only, six at a time (`results/3d/runs_untimed.csv`). R1 at 1M only for hilbert3.
- **T R1 and the extra curves** were run at t = 10 only; the t sweep covers hilbert3, peano3, morton3 and sierpinski3.
- **The DIMACS `greedy` program is 2D only**; "greedy" is our multi-fragment start (section 2).
- **nn** (nearest neighbour) only up to 40k cities: it is quadratic and only a reference.
- **No own Held-Karp bound at 1M**, so 1M results are only against greedy TO Rn.
- **Moore3 and Sierpinski3** are our own constructions (section 10), not checked against a published curve.

## 17. Reproduce

```
cd solvers && make                         # bin/*2d.exe and bin/*3d.exe (make 2d / make 3d)
py tests/regress_2d.py                     # 2D results unchanged
bin/unit3d.exe && bin/boxtest3d.exe        # 3D tests
py tests/validate_star.py                  # Waterloo tours
py bench/get_star.py                       # download the star maps (data3d/star)
py bench/gen_instances3d.py                # E3d, C3d, CT3d (data3d/instances)
py bench/hk_bounds3d.py                    # our own Held-Karp bounds (refs/hk_own.csv)
py bench/box3_plan.py                      # writes bench/box3_*.json
py bench/run.py --config bench/box3_baseline.json    # then main, tsweep, ksweep, r1, r1t,
                                           # orders, big, bigr1, extra, kd, kdscale, r1extra, bigF
py bench/run.py --config bench/box3_bigrest.json     # lengths only
py bench/box3_diag.py
py bench/box3_report.py && py bench/box3_charts.py
```

Every `run.py` call is resumable and never repeats a finished run.
