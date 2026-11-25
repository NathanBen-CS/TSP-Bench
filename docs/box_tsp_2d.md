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

**Code names** used in every table and chart: a bare name = that order only;
B = boxes (exact loops + merge); T = 2-opt; O = Or-opt; **R1** = one round of
cleanup, **Rn** = rounds repeated until nothing improves (section 2, "R1 vs Rn").
For example `greedy TO R1` (the bench's `twoopt_or`), `greedy TO Rn`, `sierpinski BTO Rn M5 t=4`.

Raw runs: `results/2d/box2/runs.csv` (878 runs after Stage D, all OK). Full tables:
[`results/2d/box2/summary.md`](../results/2d/box2/summary.md). Charts:
`results/2d/box2/box2_*.png`.

---

## 1. What was built

**Curves** (`solvers/d2/curves2d.c`): every curve gives each city a key, and a box is
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
| Strip | built (baseline) | 4 | rows | no | snake rows at full grid resolution (see §7) |
| Beta-Omega | **excluded** | | | | only an index → xy description was found (Math::PlanePath, half-plane variant); no verified xy → index for the square, closed curve |
| H-index | **not built** | | | | optional in the plan |
| Gosper, Dragon | **excluded** (per plan) | | | | can't tile a square / not one-to-one |

**Merges** (`solvers/common/merges.c`). The shared move: remove one edge from each of
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

## 2. How the cleanup works: T (2-opt) and O (Or-opt)

This is exactly what `boxtsp --cleanup oropt` does after the start tour is built,
for **every** start: B, greedy or a plain curve order. The code is in
`solvers/common/localsearch.c`. Every BT and BTO number in this report comes from it.

### What both moves share

- **Distances** are the bench's rounded integers (DIMACS rounding), the same ones
  `length.exe` uses. A move is applied only if it shortens the tour by at least 1,
  so the tour length can never go up, and the cleanup always stops.
- **Neighbour lists:** each city's **16 nearest cities** (`--k 16`, the bench's
  value), found with a grid. A move must create an edge from a city to one of
  its 16 nearest. Moves needing any other new edge are never tried.
- **Tour storage:** an array of cities in tour order, plus each city's position.
  A move that reverses part of the tour reverses **the shorter side** of the
  loop.
- **Reversal cap:** moves that would reverse more than 100,000 cities are
  skipped. The shorter side is at most half the tour, so this **never triggers**
  at 10k or 100k cities.
- **Deterministic:** no randomness, so the same input always gives the same
  tour.

### T: one 2-opt move, step by step

Take a city **a**. For each of its two tour edges, toward the next city and
toward the previous one, call the other end **b**:

1. Walk through a's 16 neighbours **c**, nearest first.
2. Gain so far: `g = d(a,b) − d(a,c)` (swap edge a–b for edge a–c).
   **As soon as `g ≤ 0`, stop.** The remaining neighbours are farther away, so
   they can't do better. This is the standard *gain criterion*.
3. Let **d** be c's tour neighbour on the same side as b is to a. The full
   saving is `g + d(c,d) − d(b,d)`: remove a–b and c–d, add a–c and b–d.
4. Keep the **best** saving found over both edges and all neighbours. This is
   "best improvement" for city a, not "first improvement".

If the best saving is positive, apply it: reverse the stretch between b and c
(or the other side of the loop, whichever is shorter).

### O: one Or-opt move, step by step

Take a city **s1**. For chains of **1, 2 and 3** consecutive cities starting at
s1, going both ways round the tour, call the chain **s1 … s2**. Call the city
just before it **a** and just after it **b**:

1. Saving from cutting the chain out and closing the gap:
   `g0 = d(a,s1) + d(s2,b) − d(a,b)`. If `g0 ≤ 0`, this chain isn't worth moving.
2. For each chain end (s1, and s2 if the chain is longer than 1), walk through
   its 16 neighbours **c**, nearest first, **stopping once `d(end,c) ≥ g0`**. Each
   c gives two possible insertion spots, between c and the city after it, or
   between the city before it and c. Call that edge **u–v**. Spots touching the
   chain itself are skipped.
3. Price both orientations of the chain:
   - same way round, `u s1 … s2 v`: `g0 + d(u,v) − d(u,s1) − d(s2,v)`
   - reversed, `u s2 … s1 v`: `g0 + d(u,v) − d(u,s2) − d(s1,v)`
4. Keep the best over all chains, ends, spots and orientations.

If it's positive, apply it. The move is carried out as 1–3 2-opt reversals, so
the tour stays one valid loop throughout.

### The order cities are checked in (work queue and rounds)

- **One round:** every city is put in a first-in, first-out queue, in tour
  order. Take a city from the front:
  - In the **T** stage, try a 2-opt move.
  - In the **O** stage, try a 2-opt move first, and only if there's none, an
    Or-opt move.
  - **If a move was made,** the cities whose edges changed (4 for 2-opt, up to 6
    for Or-opt) and the city itself go back in the queue.
  - **If not,** the city leaves the queue.
  The round ends when the queue is empty.
- **Rounds until clean:** a city is only re-checked when *its own* edges change,
  so a round can miss a move created by a neighbour's change. So after each
  round, **every city is queued again** and another round runs. The stage stops
  after the first round that makes **no move at all**.

### The two stages of a BTO run

| Stage | What runs | Recorded as |
|---|---|---|
| B | the start tour (boxes + merge) | `start_length` |
| T | 2-opt only, rounds until a round makes no move | `after_2opt_length` (the BT result) |
| O | 2-opt and Or-opt together, rounds until a round makes neither move | final tour (the BTO result) |

So BTO is "2-opt until stuck, then 2-opt + Or-opt until stuck", in that order.
The BT result is the same run measured part way through.

### What "stuck" means

When a stage ends, **no improving move of its kind is left that passes the
gain criterion and uses one of each city's 16 nearest neighbours.** Moves
outside that are never tried, so the tour is not optimal.
- **2-opt:** with complete neighbour lists on 300 cities, a brute-force check of
  every pair of edges finds no improving move.
- **Or-opt:** a brute-force check still finds up to 8 improving moves that the
  gain criterion deliberately skips (`boxtest`).

### How much work it does (means over 3 seeds)

| Start, data | 2-opt moves (rounds) | then Or-opt moves + 2-opt moves (rounds) |
|---|---|---|
| hilbert B M1 t=10, E10k | 2,067 (3.0) | 892 + 392 (3.7) |
| hilbert B M1 t=10, C10k | 2,614 (4.0) | 845 + 414 (3.3) |
| greedy, E10k | 1,382 (3.0) | 205 + 54 (3.0) |
| greedy, C10k | 1,420 (3.3) | 259 + 90 (3.0) |
| hilbert B M1 t=10, E100k | 19,870 (4.0) | 9,094 + 3,970 (4.3) |
| greedy, E100k | 12,976 (5.7) | 1,996 + 546 (3.0) |

- **Rounds** include the final round that confirmed nothing was left, so "3"
  means two rounds made moves.
- **Starting from boxes takes far more moves than starting from greedy,** about
  3–5× as many Or-opt moves. That fits the explanation of why the box start ends
  worse: its flaws are many and spread out.

### R1 vs Rn: one round, or rounds until nothing improves

Every name with a cleanup (T or O) now says how long the cleanup ran:

| Tag | Meaning | Used by |
|---|---|---|
| **R1** | **one round** of the work queue, then stop | the bench's `twoopt` = **greedy T R1** and `twoopt_or` = **greedy TO R1** (`localSearch()`) |
| **Rn** | **rounds repeated until a whole round makes no move** (usually 3–6 rounds; see the table above) | every `boxtsp` cleanup: **greedy TO Rn**, every box setting (**BT Rn**, **BTO Rn**) (`localSearchToOptimum()`) |

So **greedy TO R1** and **greedy TO Rn** start from the same greedy tour and use
the same move code, 16-nearest-neighbour lists and reversal cap. They differ in
two ways:
1. **Rounds.** R1 stops when its one round's queue is empty. Because a city is
   only re-checked when its own edges change (the blind spot described above),
   moves can still be left. Rn queues every city again and keeps going until a
   full round finds nothing.
2. **Stage order.** R1 (`twoopt_or`) mixes 2-opt and Or-opt from the first move.
   Rn runs a 2-opt-only stage first (T), then 2-opt + Or-opt (O).

**What that does to the results:**
- **Quality: close, usually a little better with Rn.** Within about 0.3% in
  Stages A and C, within 0.1% at 1M, and 6.3% against 6.4% above Held-Karp on
  TSPLIB.
- **Time: Rn costs more.** At 1M, greedy TO Rn takes 20.2 s on uniform maps
  against 11.4 s for greedy TO R1.
- **Rn is not guaranteed to be better.** The moves are done in a different order,
  so the two can end in different local optima. On dsj1000, greedy TO R1 ends
  4.2% above optimal and greedy TO Rn 8.2%. On the TSPLIB city maps R1 is ahead
  by 0.02 points (`docs/real_world.md`).
- **Neither is optimal.** Both stop when no improving move is left among each
  city's 16 nearest neighbours ("What stuck means", above).

The box settings have no R1 version: `boxtsp` always runs its cleanup to Rn. So
a box setting against **greedy TO Rn** is a like-for-like comparison of start
tours, and against **greedy TO R1** it is not.

## 3. Tests (`bin/boxtest2d.exe`, all pass)

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

## 4. Runs

| Stage | What | Runs |
|---|---|---|
| A (screen) | N = 10k, t = 10, 7 curves × 5 merges, + plain curve orders + every baseline; E and C × 3 seeds | 288 |
| B (t sweep) | N = 10k, t = 4–12, best 3 curves (Sierpinski, Moore, Knopp) × best 2 merges (M5, M2) + reference | 168 new |
| C (scale) | N = 100k, 4 finalists + reference + every baseline; E and C × 3 seeds | 66 |
| D (1 million) | N = 1k–1M, 4 finalists + reference + Hilbert M5 t = 10 + 2 plain curves + every baseline except nn; all 38 maps of the set | 200 new (section 9) |

- **Cleanup levels:** every box run used full cleanup and also recorded the tour
  before cleanup and after 2-opt. So "none", "2-opt" and "full" come from the same
  run.
- **Comparisons:** each setting is compared with the reference (Hilbert, M1,
  t = 10) **on the same instance**, as mean ± sample standard deviation over
  seeds.
- **Gaps:** gaps are to the Held-Karp bound, which exists for all 10k seeds and 2
  of 3 at 100k.

## 5. Results

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

**The same Stage A matrix as % above the Held-Karp bound, after full cleanup**
(10k, t = 10, mean of 3 maps; lower is better):

| curve | M1 | M2 | M3 | M4 | M5 | | M1 | M2 | M3 | M4 | M5 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| | **uniform** | | | | | | **clustered** | | | | |
| hilbert | 8.3 | 8.2 | 9.3 | 8.2 | 9.4 | | 17.6 | 18.0 | 21.6 | 19.0 | **13.7** |
| moore | 8.3 | 8.3 | 9.1 | 8.2 | 9.3 | | 14.4 | 14.5 | 20.5 | 18.2 | 14.6 |
| sierpinski | 7.8 | 7.8 | 8.6 | **7.6** | 8.1 | | 14.2 | 14.7 | 21.7 | 18.0 | 14.6 |
| knopp | 7.8 | 8.0 | 8.6 | **7.6** | 8.2 | | 14.1 | 14.7 | 22.2 | 20.2 | 14.5 |
| peano | 8.7 | 8.7 | 9.0 | 8.4 | 9.2 | | 20.1 | 19.7 | 22.7 | 20.1 | 14.1 |
| morton | 7.6 | 7.6 | 9.3 | 8.3 | 9.4 | | 24.3 | 22.7 | 25.9 | 24.8 | 16.0 |
| strip | 7.0 | 7.0 | 7.6 | 7.5 | 7.3 | | 35.9 | 32.8 | 30.1 | 29.9 | 30.6 |
| **mean of the 5 main curves** | 8.2 | 8.2 | 8.9 | 8.0 | 8.8 | | 16.1 | 16.3 | 21.7 | 19.1 | **14.3** |

The "merges" ranking printed by `bench/box2_report.py` (M5, M2, M1, M4, M3) is
a mean over **all** curves and both data types. It is driven by the clustered
maps, where the merges differ most, so it doesn't mean M5 is best on uniform maps.

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

## 6. Answers

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

**So which merge is best?** There is no single answer; it depends on the map
type (Stage A matrix in section 5):
- **Clustered maps: M5, clearly.** Mean 14.3% above Held-Karp over the 5 main
  curves, against 16.1% for M1, 16.3% for M2, 19.1% for M4 and 21.7% for M3.
  This is where the merges differ most, and it's why M5 tops the report's
  pooled ranking.
- **Uniform maps: all within about 1 point.** M4, M1 and M2 are slightly ahead
  (8.0–8.2%), and M3 and M5 slightly behind (8.8–8.9%).
- **With a good curve on clustered maps** (Sierpinski, Knopp, Moore), M1 and M5
  nearly tie (14.1–14.4% against 14.5–14.6%). With a weak curve (Hilbert, Peano,
  Morton), M5 is far ahead.
- **In short: M5 is the most reliable merge. M1 only matches it when the curve
  is good,** and before cleanup M1 is the weakest of all (it has no freedom to
  choose where loops are joined).

**Why the original method (Hilbert, M1) does poorly.** The weak part is the
combination, not M1 on its own:
- **M1 joins the box loops in the curve's order,** so it inherits the curve's
  mistakes. With Hilbert on clustered maps that gives 17.6%, against 14.1–14.4%
  for the same M1 with Knopp, Sierpinski or Moore (only Peano is worse, 20.1%).
- **M5 picks its joins by cost (a minimum spanning tree over the loops), not by
  curve order,** so the curve matters much less. It rescues Hilbert (17.6% →
  13.7%, the best clustered result in Stage A) and also Peano and Morton.
- **On uniform maps the original is only slightly behind** (8.3% against 7.6–7.8%
  for the best settings), because there the curve and merge hardly matter after
  cleanup.
- **Separately, every box setting loses to greedy + 2-opt + Or-opt** (9.2% on
  clustered 10k). Even the best box result, 13.7%, is 4.5 points behind. That's a
  limit of the box approach, not of one merge (see section 7).

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

**Uniform vs clustered with no cleanup at all** (`box2_raw_clustered_vs_uniform.png`:
raw greedy, each plain curve order, and each curve's boxes B with every merge at
t = 10; 10k, % above HK):
- **Raw greedy:** 16.4% uniform, 23.1% clustered (+6.7).
- **M5 is the merge that holds up on clustered data:** B M5 loses only 5.4–9.5
  points, about as little as greedy. Every other merge loses 10–29 points.
- **M4 gives the shortest raw box tours on uniform data:** 15.9–20.2%. With
  hilbert, moore and morton (all 15.9%) that is just under raw greedy's 16.4%.
  But M4 loses 12–22 points on clustered data, so there it ends 29–38%.
- **On clustered data no raw box tour beats raw greedy:** the best is hilbert
  B M5 at 28.8% vs 23.1%.
- **Boxes always help a curve's plain order,** except sierpinski and knopp with
  M1, which are slightly worse on uniform data (36.8 vs 34.5, 37.8 vs 35.6).
- **Morton's plain order** (103/113%) is the worst; its M3–M5 box tours land with
  the others.

**With vs without cleanup:**
- **Cleanup shrinks almost every difference.** Before cleanup the settings span
  hundreds of percent; after it, a few percent.

## 7. What didn't help (stated plainly)

- **No box setting beats greedy + 2-opt + Or-opt,** on either data type, at any
  size from 1k to 1M. At 10k and 100k the gap to the best box setting is about 3
  points on uniform data and 4.5–5 on clustered (% above Held-Karp); at 1M the
  best box tours are 2.9% (uniform) and 4.2% (clustered) longer (section 9).
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

## 8. Caveats

- **Few seeds:** 3 per type, and only 2 with a Held-Karp bound at 100k. The
  top finalists differ by less than the seed noise, so their exact order isn't
  meaningful.
- **Stage B's choices:** the "best 2 merges" were chosen by their mean over all
  curves (M5, M2). M1 was strong with good curves and reached the finals through
  Stage A.
- **A setting the finals missed:** the best clustered result of Stage A,
  Hilbert M5 t = 10 (13.7%), was not a Stage C finalist, because Stage B only
  swept t for the 3 best curves (Sierpinski, Moore, Knopp), and the M5 finalist
  came from there (Sierpinski M5 t = 4). It was added to Stage D afterwards
  (section 9), so it has no 100k-only Stage C result.
- **Repeated runs:** the existing solvers, the greedy baseline and the reference
  were re-run in this experiment, although the earlier box experiment had already
  run them on the same instances. All 84 repeats gave identical tour lengths. Only
  their times are new.
- **Not run:** 3D and n-dimensional versions (plan §5) were out of scope.

## 9. Scaling to 1 million cities (Stage D)

**What was run.** The key settings on maps of 1k to 1M cities, on both map sets:

- **Settings:** the 4 finalists (sierpinski M1/M2/M5, knopp M1), the reference
  (hilbert M1 t=10), hilbert M5 t=10 (the best clustered result of Stage A,
  added in a second pass; see section 8), the two best plain curve orders
  (sierpinski, moore), and every existing solver except nn. nn is O(N²) and does
  not reach 1M.
- **Uniform + clustered set** (`results/2d/box2/runs.csv`): all 38 maps, i.e. 10
  maps of each type at 1k and 3 at 10k, 100k and 1M.
- **DIMACS testbed** (`results/2d/runs.csv`): the same settings on all 80
  geometric maps up to 1M, including the 33 geometric TSPLIB maps. The 8 matrix
  maps are skipped, because the box and curve solvers need coordinates, and so
  are E3M and E10M (above 1M).
- **Reuse:** runs already made with the same program, options and map file
  were copied instead of re-run (`bench/import_runs.py`): 162 into the box set,
  and 279 into the main file. Only 200 + 441 runs were new.

**Measure.** "% longer than greedy TO Rn" on the same map, where greedy TO Rn =
greedy + 2-opt + Or-opt, the best solver here. The 1M maps generated for this
set have no Held-Karp bound, so this paired measure is used at every size. Mean
over the maps of each size:

| | uniform 1k | 10k | 100k | 1M | clustered 1k | 10k | 100k | 1M |
|---|---|---|---|---|---|---|---|---|
| greedy TO R1 | +0.0 | −0.0 | +0.0 | +0.1 | −0.3 | +0.3 | +0.1 | +0.1 |
| greedy T R1 | +1.1 | +1.0 | +1.1 | +1.0 | +1.2 | +1.6 | +1.4 | +1.4 |
| sierpinski BTO Rn M1 t=10 | +3.3 | +3.4 | +3.7 | **+2.9** | +1.5 | +4.6 | +5.8 | +6.0 |
| knopp BTO Rn M1 t=10 | +3.2 | +3.5 | +3.7 | **+2.9** | +1.7 | +4.5 | +5.9 | +6.0 |
| sierpinski BTO Rn M2 t=12 | +3.3 | +3.4 | +3.9 | +3.1 | +1.2 | +4.9 | +5.6 | +6.1 |
| sierpinski BTO Rn M5 t=4 | +4.1 | +4.5 | +4.5 | +3.9 | +3.3 | +4.3 | +3.9 | **+4.2** |
| hilbert BTO Rn M5 t=10 | +4.7 | +5.0 | +5.2 | +4.4 | +4.2 | +4.1 | +4.1 | +4.3 |
| hilbert BTO Rn M1 t=10 (reference) | +3.2 | +3.9 | +4.2 | +3.5 | +2.9 | +7.8 | +7.1 | +6.8 |
| greedy | +13.3 | +11.7 | +10.2 | +8.8 | +11.2 | +12.8 | +10.9 | +10.2 |
| hilbert B M1 t=10 (box tour, no cleanup) | +23.8 | +24.1 | +23.1 | +21.8 | +29.7 | +34.2 | +33.2 | +32.9 |
| sierpinski (curve order only) | +27.0 | +29.1 | +29.5 | +28.6 | +30.7 | +37.4 | +38.4 | +38.6 |
| moore (curve order only) | +29.8 | +31.6 | +32.5 | +31.5 | +35.7 | +39.5 | +41.2 | +41.7 |
| hilbert (spacefill) | +33.2 | +33.0 | +33.0 | +31.6 | +37.8 | +42.5 | +42.1 | +41.8 |

**Time and memory at 1M** (user CPU seconds, mean of 3 maps; peak memory on E1M.0):

| | uniform | clustered | peak memory |
|---|---|---|---|
| greedy TO Rn | 20.2 s | 14.9 s | 161 MB |
| greedy TO R1 | 11.4 s | 10.1 s | 160 MB |
| greedy T R1 | 9.8 s | 9.2 s | 160 MB |
| hilbert BTO Rn M1 t=10 | 23.1 s | 26.6 s | 107 MB |
| sierpinski / knopp BTO Rn M1 t=10 | 28.7–29.1 s | 28.8–31.2 s | 107 MB |
| sierpinski BTO Rn M5 t=4 | 26.9 s | 27.2 s | 331 MB |
| hilbert BTO Rn M5 t=10 | 26.7 s | 29.0 s | 300 MB |
| sierpinski BTO Rn M2 t=12 | 58.1 s | 56.7 s | 162 MB |
| greedy | 9.2 s | 9.1 s | 69 MB |
| plain curve orders | 1.6–2.3 s | 1.7–2.3 s | 38–42 MB |

**Held-Karp gaps on the one DIMACS 1M map (E1M.0):** greedy TO Rn 4.8%, twoopt_or
4.9%, twoopt 5.9%, knopp / sierpinski BTO Rn M1 8.1%, sierpinski BTO Rn M2 8.3%,
hilbert BTO Rn M1 8.7%, sierpinski BTO Rn M5 9.1%, hilbert BTO Rn M5 9.8%, greedy 14.3%, sierpinski curve
35.1%, moore 38.1%, spacefill 38.3%.

**TSPLIB maps** (mean % above Held-Karp over the 33 geometric maps): greedy TO Rn
6.3, twoopt_or 6.4, twoopt 7.8, sierpinski BTO Rn M5 9.3, hilbert BTO Rn M5 9.5, knopp / sierpinski BTO Rn M1
10.4, sierpinski BTO Rn M2 10.7, hilbert BTO Rn M1 11.6, greedy 18.3; the box tour and
curve orders alone 38–49.

**What it shows:**
- **Every solver keeps a steady quality on uniform maps** from 1k to 1M; each
  moves by less than about 1 point. No box setting catches up with greedy TO Rn at
  any size.
- **On clustered maps the M1/M2 box settings fall further behind as the maps
  grow:** from +1.2 to +1.7% at 1k to about +6% at 1M. Both M5 settings stay at
  about +4% at every size, so from 10k up they are the best box settings on
  clustered maps, while they are the worst on uniform maps. This matches Stages
  A–C: M5 is the most reliable merge on clustered data.
- **With M5 the curve hardly matters on clustered maps:** hilbert M5 t=10 and
  sierpinski M5 t=4 are within 0.2 points of each other at every size from 10k
  to 1M (+4.3% and +4.2% at 1M). With M1 the curve does matter: hilbert M1 is
  +6.8% at 1M against +6.0% for sierpinski and knopp M1. On uniform maps
  hilbert M5 is the weakest box setting (+4.4 to +5.2%).
- **M2 with t = 12 costs twice the time** of the others at 1M for no gain, and
  M5 uses about 3× the memory of M1 (300–331 MB compared with 107 MB).
- **The best box settings take about 1.5–2× as long as greedy TO Rn** at 1M and end
  3–6% longer. greedy TO R1, the bench's one-round `twoopt_or`, is within
  0.1% of greedy TO Rn in about 60% of the time.
- **The plain curve orders are the fastest** (about 2 s at 1M) but end 29–42%
  longer than greedy TO Rn. The box tour alone (hilbert B M1) is 22–33% longer.

Charts: `results/2d/scale_1M.png` (every solver at 1M), `results/2d/scale_growth.png`
(each cleanup solver from 1k to 1M), `results/2d/tsplib_summary.png`; made by
`bench/plot_scaling.py`.

## 10. Reproduce

```powershell
$env:PATH = "C:\msys64\ucrt64\bin;C:\msys64\usr\bin;" + $env:PATH
cd solvers; make; cd ..; .\bin\boxtest2d.exe
py bench\box2_plan.py A;  py bench\run.py --config bench\box2_stageA.json
py bench\box2_plan.py B --curves sierpinski,moore,knopp --merges M5,M2
py bench\run.py --config bench\box2_stageB.json
py bench\box2_plan.py C --settings knopp:M1:10,sierpinski:M1:10,sierpinski:M2:12,sierpinski:M5:4
py bench\run.py --config bench\box2_stageC.json
py bench\box2_report.py;  py bench\box2_charts.py
# Stage D (1 million): reuse earlier runs, then run only what is missing
py bench\box2_plan.py D
py bench\import_runs.py --from results\2d\box\runs.csv --to-config bench\box2_stageD.json --alg greedy --alg spacefill --alg twoopt --alg twoopt_or --alg greedy+2opt+oropt --alg "boxes+2opt+oropt=box:hilbert:M1:t10"
py bench\run.py --config bench\box2_stageD.json
py bench\run.py --alg greedy+2opt+oropt,box:hilbert:M1:t10,box:knopp:M1:t10,box:sierpinski:M1:t10,box:sierpinski:M2:t12,box:sierpinski:M5:t4,box:hilbert:M5:t10,curve:sierpinski,curve:moore
py bench\report.py;  py bench\plot_scaling.py
```
