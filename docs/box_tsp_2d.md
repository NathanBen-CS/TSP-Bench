# Box-TSP 2D extension: curves, merges, t

This extends the box solver (exact Held-Karp loops inside boxes, the loops
joined into one tour, then 2-opt and Or-opt) with **7 space-filling curves**,
**5 merge methods** and **box sizes t = 4–12**. It's a Karp-style heuristic
with no worst-case guarantee, and it is not a PTAS.

**Short answer:**
- **Before cleanup,** the merge method matters most: the greedy-style merges
  M4 and M5 make start tours 5–12% shorter.
- **After full cleanup,** most of that disappears:
  - **Uniform data:** every sensible setting lands within ±1% of the reference.
  - **Clustered data:** the curve and merge matter, by about 3% (closed and
    triangle curves help; M3 and M4 hurt).
- **t** barely matters after cleanup.
- **The greedy start** with the same cleanup still beats every box setting:
  by about 3 points on uniform data and 4.5 on clustered (% above the
  Held-Karp bound), at both 10k and 100k.

Raw runs: `results/box2/runs.csv` (522 runs, all OK). Full tables:
[`results/box2/summary.md`](../results/box2/summary.md). Charts:
`results/box2/box2_*.png`.

---

## 1. What was built

**Curves** (`solvers/curves.c`): every curve gives each city a key, and a box is
one curve cell, meaning all keys that share a prefix in the curve's base. So a
box is always one unbroken run of the sorted list.

| Curve | Status | Base | Cells | Closed | Note |
|---|---|---|---|---|---|
| Hilbert | built (reference) | 4 | squares | no | existing code, results unchanged |
| Moore | built | 4 | squares | yes | four Hilbert curves, rotated |
| Sierpinski | built | 2 | triangles | yes | Platzman–Bartholdi: two triangles over the square |
| Knopp | built | 2 | triangles | no | one triangle covering the square (twice its area) |
| Peano | built | 9 | squares | no | 3 × 3 split, grid 3¹³ |
| Morton | built (negative control) | 4 | squares | no | only 50.2% of steps go to a neighbouring cell (measured) |
| Strip | built (baseline) | 4 | rows | no | snake rows at full grid resolution (see §4) |
| Beta-Omega | **excluded** | | | | only an index → xy description was found (Math::PlanePath, half-plane variant); no verified xy → index for the square, closed curve |
| H-index | **not built** | | | | optional in the plan |
| Gosper, Dragon | **excluded** (per plan) | | | | can't tile a square / not one-to-one |

**Merges** (`solvers/merges.c`). The shared move: remove one edge from each of
two loops, add two edges between them, and take the cheaper of the two ways to
reconnect.

| ID | Method |
|---|---|
| M1 | snake glue (original): one cut per loop, chosen box by box with a one-box look-ahead; first box cut at a fixed spot |
| M2 | snake with every cut and direction chosen jointly by dynamic programming, closing edge included; exact for the snake shape |
| M3 | best edge pair in curve order: running loop + next box; running-loop edges near the new box only |
| M4 | greedy global cheapest: candidate merges from neighbour lists, cheapest first, until one loop is left |
| M5 | loop MST: Kruskal over loops by closest cities, best-pair merge along each tree edge |
| M6, stretch | not built (optional) |

The defaults (Hilbert, M1) reproduce the original solver exactly: 25 of 25 saved
tours were byte-identical after the refactor. The other bench solvers are
unchanged (24 of 24).

## 2. Tests (`solvers/boxtest.exe`, all pass)

- **Curve cells:** every cell gets its own position (16×16 grid, 27×27 for Peano).
- **Steps:** consecutive cells are neighbours for Hilbert, Moore, Peano and Strip
  (100%). Morton is recorded at 50.2%, as expected.
- **Closing:** Moore and Sierpinski really close.
- **Triangle curves:** each of 1,024 triangles maps back to its own position, and
  consecutive triangles share an edge.
- **Boxes, every curve, on uniform, clustered and duplicate-heavy data:** each box
  is exactly one curve cell with no other city in it, and only depth-capped boxes
  go over t.
- **Plan's merge test:** two unit squares 2 apart. M2, M3, M4 and M5 give exactly
  **10**, the optimum. M1 gives 12, because its first box's cut is fixed.
- **Added cost:** each method's reported added cost equals the independently
  recounted (tour − loops) in all 45 combinations (3 curves × 3 data types × 5
  merges), and every tour is a valid permutation.
- **M2:** equals the brute-force best snake on 16 small instances, and is never
  worse than M1.

## 3. Runs

| Stage | What | Runs |
|---|---|---|
| A (screen) | N = 10k, t = 10, 7 curves × 5 merges, + plain curve orders + every baseline; E and C × 3 seeds | 288 |
| B (t sweep) | N = 10k, t = 4–12, best 3 curves (Sierpinski, Moore, Knopp) × best 2 merges (M5, M2) + reference | 168 new |
| C (scale) | N = 100k, 4 finalists + reference + every baseline; E and C × 3 seeds | 66 |

- **Cleanup levels:** every box run used full cleanup and also recorded the tour
  before cleanup and after 2-opt. So "none", "2-opt" and "full" come from the same
  run.
- **Comparisons:** each setting is compared with the reference (Hilbert, M1,
  t = 10) **on the same instance**, as mean ± sample standard deviation over
  seeds.
- **Gaps:** gaps are to the Held-Karp bound, which exists for all 10k seeds and 2
  of 3 at 100k.

## 4. Results

**Curve × merge at 10k, t = 10, % vs reference** (negative = shorter; mean of 3
seeds; the `box2_curve_merge.png` heatmap shows the same):

| | before cleanup, uniform | after full cleanup, uniform | before cleanup, clustered | after full cleanup, clustered |
|---|---|---|---|---|
| best curves (Sierpinski, Moore, Knopp) with M1/M2 | −2.0 to +6.6 | −0.5 to −0.0 | −3.2 to +2.0 | **−2.5 to −3.0** |
| any curve but Strip, with M4 | −7.1 to −10.4 | −0.7 to +0.1 | −5.5 to −11.8 | +0.3 to +6.1 |
| any curve but Strip, with M5 | −2.2 to −5.1 | −0.2 to +1.0 | −9.4 to −12.1 | −1.3 to −3.3 |
| any curve but Strip, with M3 | −2.5 to −5.0 | +0.2 to +1.0 | −2.1 to +4.0 | **+2.4 to +7.0** |
| Peano | −9.6 to +3.2 | +0.1 to +0.8 | −11.1 to +4.2 | −2.9 to +4.3 |
| Morton (negative control) | −10.4 to +29.8 | −0.7 to +1.0 | −9.8 to +27.7 | −1.3 to +7.0 |
| Strip | +762 to +1979 | −0.6 to −1.2 | +683 to +2726 | +10.5 to +15.6 |

**Best settings vs the baselines, % above the Held-Karp bound, after full cleanup:**

| | 10k uniform | 10k clustered | 100k uniform | 100k clustered |
|---|---|---|---|---|
| reference: Hilbert, M1, t = 10 | 8.3 | 17.6 | 8.3 | 17.6 |
| Sierpinski, M1, t = 10 | 7.8 | 14.2 | 8.0 | 16.3 |
| Knopp, M1, t = 10 | 7.8 | 14.1 | 8.0 | 16.3 |
| Sierpinski, M2, t = 12 | 7.8 | 14.5 | 8.2 | 15.8 |
| Sierpinski, M5, t = 4 | 8.9 | **13.9** | 8.9 | **14.3** |
| **greedy + 2-opt + Or-opt** | **4.2** | **9.2** | **4.2** | **9.5** |
| twoopt_or (existing) | 4.2 | 9.5 | 4.2 | 9.6 |

## 5. Answers

**Does the merge (glue) method matter?**
- **Before cleanup: yes, a lot.** M4 and M5 cut the box tour by up to 10–12%.
  They join loops wherever it's cheapest instead of in curve order, so the long
  seams get fewer.
- **After cleanup on uniform data: no.** All methods end within about 1% of each
  other.
- **After cleanup on clustered data: yes, by several points** (Hilbert: from
  −3.3 with M5 to +3.4 with M3). M5 is the most reliable (−1.3 to −3.3 vs reference for every curve). M1 and M2 do well
  with good curves. **M3 is the worst, and M4 isn't far behind:** their much
  better start tours end *longer* after cleanup. A shorter start tour doesn't
  mean a shorter final tour.

**Does the curve matter?**
- **On uniform data: hardly** (±1%).
- **On clustered data: yes.** The closed curves (Moore, Sierpinski) and the
  triangle curves (Sierpinski, Knopp) beat Hilbert by about 3% with M1/M2.
- **Morton and Peano** are worse unless M5 repairs their joins.
- **Strip** is far the worst on clustered data.

**Does t matter?**
- **After cleanup: no.** Across t = 4–12, each curve and merge setting moves by
  less than about 1 point, without a consistent direction. That's within the
  seed-to-seed spread (`box2_t_sweep.png`).
- **Time:** t only costs time. At 100k, Sierpinski M2 with t = 12 takes 4.4–4.8 s
  against 1.7 s for the reference, because M2's dynamic programming grows with
  box size.

**Uniform vs clustered:**
- **Every choice studied** matters more on clustered data.
- **Every solver** loses about 5–9 points on clustered data; the plain curve
  orders lose 15–19 (Morton 10) (`box2_clustered_vs_uniform.png`).
- **The best box settings** lose about 5–7 points, close to greedy + cleanup (5.0).
  The reference loses 9.3.

**With vs without cleanup:**
- **Cleanup shrinks almost every difference.** Before cleanup the settings span
  hundreds of percent; after it, a few percent.

## 6. What didn't help (stated plainly)

- **No box setting beats greedy + 2-opt + Or-opt,** on either data type, at 10k
  or 100k. The gap to the best box setting is about 3 points on uniform data and
  4.5–5 on clustered (% above Held-Karp).
- **M3 (best edge pair in curve order),** the plan's requested method, is the
  worst merge after cleanup on clustered data (+2.4 to +7.0 vs reference). Its
  start tours are good, but they leave cleanup stuck in worse local optima.
- **M4's large improvement before cleanup** mostly vanishes after it, and on
  clustered data it ends behind M1.
- **t** had no reliable effect after cleanup.
- **Strip:** at full grid resolution each row holds about one city, so the strip
  order is essentially "sort by y", and each box is a full-width band. That makes
  it catastrophic before cleanup. This is the plan's definition, not the classic
  √N-strip heuristic. Oddly, on uniform data it ends 1.2% *better* than the
  reference after cleanup; on clustered data it ends 10–16% worse.

## 7. Caveats

- **Few seeds:** 3 per type, and only 2 with a Held-Karp bound at 100k. The
  top finalists differ by less than the seed noise, so their exact order isn't
  meaningful.
- **Stage B's choices:** the "best 2 merges" were chosen by their mean over all
  curves (M5, M2). M1 was strong with good curves and reached the finals through
  Stage A.
- **Repeated runs:** the existing solvers, the greedy baseline and the reference
  were re-run in this experiment, although the earlier box experiment had already
  run them on the same instances. All 84 repeats gave identical tour lengths. Only
  their times are new.
- **Not run:** 3D and n-dimensional versions (plan §5) were out of scope.

## 8. Reproduce

```powershell
$env:PATH = "C:\msys64\ucrt64\bin;C:\msys64\usr\bin;" + $env:PATH
cd solvers; make; .\boxtest.exe; cd ..
py bench\box2_plan.py A;  py bench\run.py --config bench\box2_stageA.json
py bench\box2_plan.py B --curves sierpinski,moore,knopp --merges M5,M2
py bench\run.py --config bench\box2_stageB.json
py bench\box2_plan.py C --settings knopp:M1:10,sierpinski:M1:10,sierpinski:M2:12,sierpinski:M5:4
py bench\run.py --config bench\box2_stageC.json
py bench\box2_report.py;  py bench\box2_charts.py
```
