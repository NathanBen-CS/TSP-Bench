# Hull rings ("onion peeling") for 2D TSP

An empirical study of a convex-layer start tour, compared with the box methods
and greedy on the same 2D bench. **2D only. No worst-case guarantee is claimed or
tested.**

*Sections 2 and 3 (plan, hypotheses) were written on 2026-11-24 before any ring code
or run and are left as written; what changed afterwards is listed at the end of
section 2. Everything from section 4 on was written after the runs.*

**Short answer.** Hull rings never beat greedy TO Rn: 0 of 210 single runs on any
map, at any size, in any family. On uniform maps from 1k to 100k, cleaned hull tours
are level with the box methods: +7.5–7.6% above the Held-Karp bound against +6.7 to
+8.6% for the 7 curves at their best merge, and +4.2% for greedy TO Rn. On clustered
maps they are far behind: +22.6% at 10k and +30.7% at 100k, against about +14% for the
curves and +9.5% for greedy TO Rn. Before cleanup the hull tour is several times the
optimum and gets worse with N, from 4× the bound at 1k to 14× at 1M, as predicted
(H2). At 1M the benchmark's cleaned hull tour falls apart (+19.8%). A test run without the 2-opt reversal cap gets +7.5%, in line with 100k, so the collapse comes from that cap and not from the rings (section 11). The one
family where hull comes close to greedy TO Rn is concentric circles: +2.5 ± 2.1%,
the best mean of any start there, but within seed noise of the best box setting.

Data: [`results/2d/hull/summary.md`](../results/2d/hull/summary.md) (every table),
`results/2d/hull/runs.csv` (2,288 runs, all OK), charts `results/2d/hull_clean_vs_raw.png`
and `results/2d/hull_shapes.png`; pictures of hull order next to the curves in
`results/2d/curves/hull.png` and `all_curves.png`. Code names as in [`box_tsp_2d.md`](box_tsp_2d.md): a
bare name is that order only; B = boxes (here rings) + merge; T = 2-opt; O = Or-opt;
R1 = one cleanup round; Rn = rounds until nothing improves. **greedy TO Rn** is the
baseline.

---

## 1. The method (as given)

1. The convex hull of all cities is the outermost ring.
2. Remove those cities, take the hull of the rest: the next ring. Repeat until no
   city is left (convex layers, "onion peeling").
3. Each ring is a loop. For points in convex position the hull polygon is the
   optimal tour of those points, so a ring is the same kind of object as a box: an
   exact loop through a group of cities, with no Held-Karp needed.
4. Merge inward: remove one edge from the running loop and one from the next ring,
   and reconnect the four loose ends the cheaper of the two ways. Continue from
   the outside in until one tour is left.
5. Clean up with the same 2-opt + Or-opt (Rn) as every other method.

## 2. Plan (written before the runs)

**Route: rings plug into the box machinery.** The merges in
`solvers/common/merges.c` only see a `BoxSet`: a list of groups, each with its
loop in order. They never look at box geometry, so a ring can stand in for a box
unchanged. `boxtsp2d` gets a `--split rings` option. With it, the boxes are the
rings, outermost first, and each ring's loop is its hull order (no t, no
Held-Karp). Everything after that is the existing code: merges M1–M5, cleanup,
the stats the program prints, `bench/run.py` and its tour check.

**Files touched:**

| file | change |
|---|---|
| `solvers/d2/rings2d.c`, `rings2d.h` (new) | convex layers by repeated Andrew monotone chain over the cities still alive, presorted once by (x, y); collinear boundary points and duplicates kept in the ring; the hull-order reference tour; debug checks |
| `solvers/d2/boxtsp.c` | `--split rings`, `--start hull` (hull order alone), `--rounds 1` (R1 cleanup, as in 3D), `--check` (debug assertions). Nothing printed changes for the existing options |
| `solvers/common/merges.c`, `merges.h` | merge **H0** added with an id after M5, outside the M1–M5 range the tests loop over, so `boxtest2d` and the 3D program are unchanged |
| `solvers/Makefile` | `rings2d.o` in the 2D library; new test `ringtest2d` |
| `tests/ringtest2d.c` (new) | the ring checks below |
| `bench/gen_shapes2d.py` (new) | annulus, concentric circles, Gaussian blob, thin rectangle |
| `bench/hull_plan.py`, `hull_report.py`, `hull_charts.py` (new) | run configs, tables, charts; results in `results/2d/hull/` |

**Merge settings for rings:**
- **H0** (the literal version): merge ring k into the running loop using one
  surviving running-loop edge with both ends in ring k−1 and one edge of ring k,
  the cheapest of all such pairs. If ring k−1 has no surviving edge (a 1-city
  ring), fall back to the M3 candidates: running-loop edges at cities that are
  among the K nearest neighbours of ring k.
- **H-M1 … H-M5:** the existing merges on the rings, outermost first.
- **Hull order alone** (reference, like a plain curve order): ring 0 from its
  first city counter-clockwise, then each next ring entered at its city nearest
  to where the previous ring ended, once round counter-clockwise. No merge
  optimisation.

**Data and seeds:** the existing 2D maps, uniform (E) and clustered (C), 3 seeds
(.0 .1 .2) at 1k, 10k and 100k, and one map (.0) at 1M. Shape families at 10k,
3 seeds each.

**Bound:** the published DIMACS Held-Karp bound (`refs/bounds`). It exists for
every 1k and 10k map used, for seeds .0 and .1 at 100k, and for E1M.0. It does
not exist for E100k.2, C100k.2 or C1M.0, so the "% above Held-Karp" charts use
2 maps at 100k and no clustered point at 1M. Our own bound
(`bin/hkbound2d.exe`) is O(N²) and is 0.27% below the published one on C100k.0,
so it is not mixed in. Comparisons "% longer than greedy TO Rn" use every map.

**What changed after the plan:**
- **Strip at 1M:** only M1 was run (the user's call). Strip's tour starts so badly
  that one 1M run took 12–19 minutes, and strip was already known to be the worst
  curve on clustered data. So the charts have no strip point at 1M, and strip has
  no R1 runs.
- **Reused runs:** 292 box and greedy runs were copied from `results/2d/box2/runs.csv`
  (same program, options and map file; the regression test shows the tours are
  identical) instead of being run again. They keep their original timestamps.
- **Reversal-cap test (section 11):** not in the plan. It was added after the 1M result
  and is one untimed run of a separately built program.

## 3. Hypotheses (written before any run)

- **H1 (the user's):** hull rings are very good on uniform data and very bad on
  clustered data.
- **H2 (the counter-prediction):** on uniform data the number of layers grows like
  n^(2/3), so a ring holds about n^(1/3) cities. Neighbouring cities on one ring are
  then about n^(−1/3) apart, while a city's nearest neighbour is about n^(−1/2)
  away. If so, the raw (no cleanup) hull tour is several times optimal, its gap
  grows with n, and it may be the worst raw method of all. After cleanup it may
  still trail the box methods and greedy TO Rn.
- **H3:** clustered data hurts badly, because outer rings span the empty gaps
  between clusters again and again.
- **H4:** the method does best on data that is already ring-shaped (annulus,
  concentric circles, a single Gaussian blob) and at small n.

Each is answered with numbers in section 9.

---

## 4. How it was built

**Route: rings as boxes, as planned.** It needed no forcing: the merges only read a
`BoxSet` (each group's cities and its loop), so `ringBuild` fills one with the rings,
outermost first, and the rest of `boxtsp2d` runs unchanged.

| piece | where | notes |
|---|---|---|
| ring construction | `solvers/d2/rings2d.c` (`ringBuild`) | presort by (x, y) once; equal coordinates form one point whose cities stay next to each other; one Andrew monotone-chain pass per ring over the points still alive, popping only on a strict right turn, so collinear boundary points stay in the ring in boundary order; 1–2 points left, or all on one line, make the last ring |
| ring order | same | counter-clockwise, starting at the lowest of the leftmost points; the ring's loop is that order (no Held-Karp, no t) |
| hull order alone | `hullOrderTour` | ring 0 from its first city; each next ring entered at its city nearest to where the ring before ended, once round |
| merge H0 | `solvers/common/merges.c` (`mergeH0`) | the literal rule: cheapest pair of a running-loop edge with both ends in ring k−1 and an edge of ring k, both reconnections; if ring k−1 has no such edge left (a 1-city ring), M3's candidates |
| debug checks | `ringCheck`, `boxtsp2d --check` | every city once; every ring convex and counter-clockwise (no right turn, turning angles add up to one full turn); up to 20,000 cities, every city of a later ring strictly inside |
| options | `solvers/d2/boxtsp.c` | `--split rings`, `--start hull`, `--merge H0`, `--rounds 1`, `--check`; the old options print exactly what they did before |
| tests | `tests/ringtest2d.c` | section 5 |

H0's id comes after M5, outside the M1–M5 range the existing tests loop over, so it
changed none of their output. Merges M1–M5 needed no change for rings. Rings of 1 or 2
cities were already handled, since boxes can hold 1 or 2 cities: a 1-city loop points
to itself and a merge inserts it.

## 5. Checks

- **2D regression:** `py tests/regress_2d.py` gives 244 tours compared, 0 different, `boxtest` 120 PASS:
  every earlier 2D tour is byte-identical. The 3D tests (`boxtest3d`, `unit3d`) pass, and a
  3D run gives the length recorded in `results/3d/runs.csv`.
- **`bin/ringtest2d.exe`, 16 of 16 pass:**
  - **Convex position:** 300 sets of 3–12 cities, one per sector of a circle. Each is one
    ring, and the ring is exactly as long as the Held-Karp optimum.
  - **Random sets:** 300 sets of 4–12 cities, some with duplicates. Every ring of 4+ cities
    is an optimal loop of its cities. M1–M5, H0 and hull order all give permutations, and
    their length equals the ring loops plus the merge cost the method computed.
  - **Grids:** a 10×10 grid gives rings of 36, 28, 20, 12 and 4 cities; a 9×9 grid gives
    32, 24, 16, 8 and 1. So collinear boundary points stay in their ring.
  - **Small and degenerate cases:** points on one line (with a duplicate) make one ring that
    goes out and back. 1, 2 and 3 cities work, including 3 on a line.
  - **Larger maps:** 5,000 uniform, 5,000 clustered and 3,000 duplicate-heavy cities pass
    every check, including strict containment.
- **Every benchmark map** (E and C 1k–100k, E1M.0, C1M.0, the 12 shape maps) passes
  `--check`, with containment up to 20,000 cities: 32 maps, 0 failures.
- **Every run's tour** was checked by `codes/length.exe` (a permutation, and the program's
  own length equal to the recomputed one): 2,288 runs, all OK.
- **First test failure, and its cause:** the convex-position test first used fully random
  angles. Two cities can then be so close in angle that rounding to integers pulls one
  inside, so the set is not in convex position (one set gave 2 rings). The test now puts
  one city in each sector.

## 6. Rings on the maps

Picture: `results/2d/curves/hull.png` draws hull order (no cleanup) on a 16 × 16 grid,
on E1k.0 and on C1k.0, in the same style as the curve pictures; `all_curves.png` puts it
next to the seven curves on the grid. On C1k.0 the inner rings run between clusters.

From the H0 runs. "Ring edge / NN" = mean length of a ring edge divided by the mean
distance from a city to its nearest neighbour. The ring time is the program's own clock.

| maps | n maps | rings | cities per ring | rings / N^(2/3) | ring edge / NN | ring time |
|---|---|---|---|---|---|---|
| uniform 1k | 3 | 49 ± 2 | 20.6 | 0.487 | 6.5 | <0.01 s |
| uniform 10k | 3 | 226 ± 1 | 44.2 | 0.487 | 9.6 | 0.05 s |
| uniform 100k | 3 | 1,048 ± 2 | 95.4 | 0.486 | 14.1 | 2.6 s |
| uniform 1M | 1 | 4,860 | 205.8 | 0.486 | 20.7 | 126 s |
| clustered 1k | 3 | 57 ± 6 | 17.5 | 0.573 | 17.3 | <0.01 s |
| clustered 10k | 3 | 299 ± 6 | 33.5 | 0.644 | 29.4 | 0.07 s |
| clustered 100k | 3 | 1,532 ± 19 | 65.3 | 0.711 | 48.2 | 3.8 s |
| clustered 1M | 1 | 7,498 | 133.4 | 0.750 | 75.2 | 185 s |

- **Uniform maps follow N^(2/3) almost exactly:** 0.486–0.487 · N^(2/3) layers at every
  size. A ring holds 2.06 · N^(1/3) cities.
- **The ring-edge ratio grows like N^(1/6):** from 6.5 to 20.7, about ×1.47 per tenfold N,
  which is exactly 10^(1/6). This is the n^(−1/3) against n^(−1/2) of H2.
- **Clustered maps have more layers than N^(2/3) predicts, and the excess grows with N**
  (0.57 → 0.75 · N^(2/3)). Their ring edges are 3–4× longer relative to the nearest-neighbour
  distance than on uniform maps.
- **Construction time grows like N^(5/3):** about ×50 from 100k to 1M; 10^(5/3) is 46.
  At 1M that is 126 s (uniform) and 185 s (clustered), against 18 s and 16 s for the whole
  greedy TO Rn run. It finished, so 1M was run with the planned algorithm.

## 7. The two graphs

![Four line charts in one row, % above the Held-Karp bound against the number of cities
(1k to 1M, log scale). Left pair, after 2-opt + Or-opt: on uniform maps the orange hull line
sits among the grey curve lines at +6 to +8% up to 100k, then jumps to +19.8% at 1M (a hollow point marks +7.5% for one run without the
2-opt reversal cap), while
black greedy TO Rn stays near +4 to +5%. On clustered maps hull rises from +9.7% at 1k to
+30.7% at 100k, above every curve except strip, and greedy TO Rn stays at +5 to +9.5%.
Right pair, the same runs before cleanup, log scale: the hull line rises from about 330% to
1,350% above the bound on uniform maps and to 2,555% on clustered maps, far above every curve
except strip; greedy alone stays at +14 to +21%.](../results/2d/hull_clean_vs_raw.png)

*`results/2d/hull_clean_vs_raw.png`. Left: start + 2-opt + Or-opt, rounds until clean
(TO Rn). Right: the start tours of the same runs, no cleanup, log scale. Each curve
(grey) and hull (orange, best of H0 and M1–M5) at its best merge per size, chosen on the
cleaned mean, at t = 10 and K = 16. The raw panels use the same runs. Black: greedy TO Rn
(left) and greedy alone (right). Dotted: plain Hilbert order and hull order alone. Hollow
point: hull M1 at 1M without the 2-opt reversal cap (one run, section 11). Mean
± 1 sd over the maps with a published bound: 3 at 1k and 10k, 2 at 100k, 1 at 1M (uniform
only). Codes: B = boxes (for hull, rings) + merge, T = 2-opt, O = Or-opt, Rn = rounds until
nothing improves. Data: `results/2d/hull/hull_clean_vs_raw.csv`.*

**Key numbers** (% above Held-Karp, cleaned / raw):

| maps | hull (best merge) | best curve | greedy TO Rn / greedy alone |
|---|---|---|---|
| uniform 1k | +6.20 / +329 (M1) | +6.15 (morton M2) | +4.05 / +17.2 |
| uniform 10k | +7.52 / +559 (M1) | +7.04 (strip M2) | +4.20 / +16.4 |
| uniform 100k | +7.63 / +883 (M1) | +6.69 (strip M1); next morton M1 +7.68 | +4.18 / +15.1 |
| uniform 1M | +19.78 / +1,351 (M1); no reversal cap: +7.53 | +7.85 (sierpinski / knopp M4) | +4.76 / +14.3 |
| clustered 1k | +9.65 / +665 (M2) | +6.91 (knopp M2) | +5.35 / +17.3 |
| clustered 10k | +22.65 / +1,255 (M4) | +13.68 (hilbert M5) | +9.16 / +21.1 |
| clustered 100k | +30.69 / +2,555 (M3) | +14.20 (knopp M5) | +9.51 / +21.1 |

- **Merging matters little for hull.** Raw hull with the best merge is only 2–15% shorter
  than hull order alone. After cleanup the merges are within about 0.4 points of each other
  on uniform maps from 10k up. M1 (the snake) is the best hull merge on uniform maps at every
  size. H0, the literal version, is not the best anywhere on uniform maps, but is within seed
  noise of M1 up to 100k.
- **Every merged hull start is worse than every curve's start except strip,** and the gap
  grows with N. Strip boxes are the only raw starts worse than hull: +2,371% at 10k and
  +8,042% at 100k on uniform maps.

## 8. Where it thrives and fails: map families at 10k

% longer than greedy TO Rn on the same map, 3 maps per family. Best hull merge and best
box setting (of 7 curves × M1–M5) are chosen on the TO Rn mean; raw = the start tours of
those runs. Chart: `results/2d/hull_shapes.png`.

| maps | hull TO Rn | hull raw | best box TO Rn | best box raw | greedy alone |
|---|---|---|---|---|---|
| annulus | +3.45 ± 1.26 (M2) | +539 | +3.04 ± 0.19 (knopp M4) | +15 | +11.1 |
| concentric circles | **+2.51 ± 2.14** (M3) | +61 | +3.33 ± 2.65 (morton M5) | +44 | +2.3 |
| Gaussian blob | +3.20 ± 0.61 (M2) | +396 | +2.55 ± 0.56 (strip M1) | +1,788 | +11.1 |
| thin rectangle 10:1 | +2.88 ± 0.21 (M2) | +1,061 | +1.96 ± 0.14 (strip M1) | +7,605 | +12.5 |
| uniform (E) | +3.18 ± 0.29 (M1) | +532 | +2.73 ± 0.37 (strip M2) | +2,272 | +11.7 |
| clustered (C) | +12.36 ± 2.24 (M4) | +1,141 | +4.14 ± 1.08 (hilbert M5) | +18 | +11.0 |

**Diagnostics** (from the H0 runs):

| maps | rings | cities per ring | rings / N^(2/3) | ring edge / NN |
|---|---|---|---|---|
| annulus | 160 | 62.6 | 0.344 | 9.8 |
| concentric circles | 40 | 252.1 | 0.085 | 3.6 |
| Gaussian blob | 254 | 39.4 | 0.547 | 7.7 |
| thin rectangle | 226 | 44.2 | 0.487 | 18.4 |
| uniform | 226 | 44.2 | 0.487 | 9.6 |
| clustered | 299 | 33.5 | 0.644 | 29.4 |

What it shows:
- **Concentric circles are the one family where rings are nearly the right tour pieces.**
  The 20 circles become 40 rings: rounding to integers splits most circles into two layers.
  A ring edge is only 3.6 nearest-neighbour distances long (9.6 on uniform maps), and the
  raw hull tour is +61% above greedy TO Rn instead of +500% or more. After cleanup hull has
  the best mean of any start there, but within seed noise of the best box (sd 2.1 and 2.6).
  Greedy alone is already within 2.3% of greedy TO Rn on these maps, so there is little
  left to win.
- **Annulus and blob: no advantage.** Their ring edges are as long relative to the
  nearest-neighbour distance as on uniform maps (9.8 and 7.7), and hull ends as far behind as
  on uniform maps.
- **Thin rectangle:** the ring count matches the uniform square (226), but each ring is
  stretched along the long side. The ring-edge ratio doubles (18.4), and the raw tour is
  twice as bad as on the square.
- **Across these families, the raw hull gap follows the ring-edge ratio:** lowest for circles
  (3.6), then blob (7.7), uniform and annulus (9.6–9.8), rectangle (18.4) and clustered
  (29.4). That is a correlation over six families, not a tested cause. After cleanup the
  order mostly washes out, except on clustered maps.

## 9. The hypotheses, after the runs

- **H1 (amazing on uniform, terrible on clustered): partly supported.**
  - **Clustered half, confirmed:** after cleanup hull is +22.6% above Held-Karp at 10k and
    +30.7% at 100k. That is about twice the curves (+13.7 to +16.0) and three times greedy TO
    Rn (+9.2 to +9.5). Against greedy TO Rn it is +12.4% at 10k, where the best box setting
    is +4.1%.
  - **Uniform half, not supported:** hull is not amazing there. From 1k to 100k it is level
    with the box methods (+7.6% at 100k, against +6.7 to +8.6% for the curves), never beats
    greedy TO Rn (best mean +2.1%, 0 of 54 single runs). At 1M the benchmark run collapses to
    +19.8%, but that is the 2-opt reversal cap (section 11); without it, one run gives +7.5%.
- **H2 (N^(2/3) layers, raw tour far from optimal and getting worse): confirmed, except
  "worst raw method".**
  - **Layers:** 0.486 · N^(2/3) from 10k to 1M, with 2.06 · N^(1/3) cities per ring.
  - **Ring edges:** 6.5 to 20.7 nearest-neighbour distances, growing like N^(1/6).
  - **Raw tour:** 4.3× the bound at 1k, 6.6× at 10k, 9.8× at 100k and 14.5× at 1M (hull M1).
  - **Worst raw method:** no. Strip boxes are worse (+8,042% at 100k), but hull is worse than
    every other curve and than greedy alone at every size.
  - **"Trails boxes and greedy after cleanup":** it trails greedy TO Rn everywhere. Against the
    boxes it is level from 1k to 100k (better than 5 of the 7 curves' best settings at 100k).
    At 1M it is behind with the bench's reversal cap and level without it (one run).
- **H3 (clustered hurts, because rings span the gaps): confirmed that it hurts; the mechanism
  is supported but not directly tested.**
  - **It hurts:** clustered maps have more and smaller rings (0.64–0.75 · N^(2/3)), raw hull
    tours twice as bad as on uniform maps, and the largest cleaned gaps in the study.
  - **The mechanism:** ring edges are 3–4× longer relative to the nearest-neighbour distance
    (29.4 against 9.6 at 10k, 48.2 against 14.1 at 100k). That fits rings that cross the gaps
    between clusters, but no count of gap-crossing edges was made.
- **H4 (best on ring-shaped data and at small N): partly supported.**
  - **Small N:** yes. Against greedy TO Rn the best hull mean is +2.1% at 1k, +3.2% at 10k and
    +3.4% at 100k on uniform maps. On clustered maps it is +4.1% at 1k and +12–20% above that.
  - **Concentric circles:** yes, closest to greedy TO Rn there (+2.5 ± 2.1%).
  - **Annulus and Gaussian blob:** no. They do no better than the uniform square (+3.5 and
    +3.2, against +3.2), and both trail the best box setting.

## 10. Time, and less cleanup

**Time** (user CPU seconds, whole program; mean over maps; 1M is one map):

| maps | greedy TO Rn | ring construction | hull M1 (whole run) | hull M2 | hull order |
|---|---|---|---|---|---|
| uniform 10k | 0.08 | 0.05 | 0.15 | 0.63 | 0.06 |
| uniform 100k | 1.08 | 2.57 | 10.2 | 22.2 | 2.65 |
| uniform 1M | 18.1 | 126 | 527 | 866 | 125 |
| clustered 100k | 1.08 | 3.84 | 10.2 | 16.0 | 3.85 |
| clustered 1M | 15.6 | 185 | 561 | 691 | 193 |

- **Hull costs about 10× greedy TO Rn at 100k and 30× at 1M.** Most of that is 2-opt
  repairing the bad start (hull M1 on E1M.0: 407 s of its 527 s), then the ring
  construction.
- **M2 is the slowest merge on rings:** its dynamic programme grows with the square of the
  ring size, and rings are much larger than boxes.
- **The greedy and box times were partly measured in October** (reused runs), on the same
  laptop, pinned to core 0.

**Less cleanup** (% longer than greedy TO Rn; best hull merge and best box setting; strip
has no R1 runs; full table in `summary.md` section 7):

| maps | greedy TO R1 | hull TO R1 | box TO R1 | greedy T Rn | hull T Rn | box T Rn |
|---|---|---|---|---|---|---|
| uniform 10k | −0.03 ± 0.19 | +3.09 ± 0.45 | +3.04 ± 0.32 | +0.96 | +7.36 | +6.07 (strip) |
| uniform 100k | +0.03 ± 0.05 | +3.36 ± 0.16 | +3.50 ± 0.21 | +0.89 | +7.42 | +5.43 (strip) |
| clustered 10k | +0.32 ± 0.13 | +15.83 ± 2.27 | +4.12 ± 0.99 | +1.44 | +21.00 | +8.51 |
| clustered 100k | +0.07 ± 0.09 | +21.76 ± 1.63 | +4.16 ± 0.73 | +1.29 | +26.94 | +7.98 |

- **Less cleanup does not change the picture.** On uniform maps, hull with one TO round is
  level with the best box setting other than strip (within one sd). With 2-opt only, hull
  trails strip but ties or beats the other curves. On clustered maps hull stays 4–5× further
  behind greedy TO Rn than the boxes, at every cleanup level.
- **One round of TO is as good as Rn for greedy** (within 0.3%), as in the box studies.

## 11. Why hull collapses at 1M

At 1M on the uniform map, the cleaned hull tour is +19.8% to +26.7% above the bound,
depending on the merge, against +7.6–7.9% at 100k. Strip M1, the other very bad start,
jumps the same way (+6.7% at 100k → +31.3% at 1M); every other curve stays at +8 to +9%.

**The guess:** the bench's 2-opt skips any move that would reverse more than 100,000
cities (`MAX_REVERSAL` in `solvers/common/localsearch.c`, kept to stop O(N²) behaviour). At
100k a move never reverses more than 50,000 cities, so the cap never applies. At 1M a move
can need up to 500,000. A start that is 14 times too long needs many long reversals to
repair, so the cap would hurt it far more than a good start.

**The test:** one run of hull M1 on E1M.0 with a separately built `boxtsp2d`
(`-DMAX_REVERSAL=100000000`, i.e. no cap), everything else the same. It was not timed
against the benchmark: it ran on core 2 while other runs used core 0. The tour was checked by
`codes/length.exe`. Data: `results/2d/hull/reversal_cap_test.csv`.

| E1M.0, hull M1 | start | after 2-opt | 2-opt rounds | 2-opt time | final | vs Held-Karp |
|---|---|---|---|---|---|---|
| cap 100,000 (the benchmark run) | 10,284,714,366 | 881,762,368 | 16 | 407 s | 848,859,173 | +19.78% |
| no cap | 10,284,714,366 | 791,376,078 | 8 | 1,882 s | 762,075,030 | **+7.53%** |

- **The cap explains the collapse.** Without it, hull M1 at 1M ends at +7.53%, against +7.63%
  at 100k. It is shown as the hollow point in the chart.
- **The price is time.** 2-opt took 1,882 s instead of 407 s, about 100× greedy TO Rn's whole
  run (18 s).
- **This was one map and one merge, so it doesn't compare hull with the curves at 1M.**
  Their 1M runs used the cap too. Their good starts need far fewer long reversals, but whether
  the cap costs them anything was not tested. The 1M strip jump (+31.3%) fits the same
  explanation, but that was not tested either.

## 12. Limits

- **Sample sizes:**
  - **Maps:** 3 per group at 1k, 10k and 100k, and 1 at 1M.
  - **Bound:** 2 maps at 100k and 1 at 1M (uniform only) have the published Held-Karp bound.
  - **Clustered 1M:** no bound, so only the "vs greedy TO Rn" tables have it.
  - **Families:** 3 maps each, at 10k only.
- **Best-of selection:** picking the best of 5 or 6 merges on the same maps favours every
  method a little. It is the same rule for hull and for the curves.
- **The bench's integer distances:** rings are computed on the exact coordinates, but every
  length uses DIMACS rounding. The convex-position test shows no difference on the sets
  tested.
- **Circles in rounded coordinates** are not exact circles (40 rings for 20 circles). Exact
  circles might do better or worse; not tested.
- **Not tested as causes:**
  - the "rings span cluster gaps" mechanism (H3);
  - the link between the ring-edge ratio and the raw gap (section 8).
  The reversal-cap explanation was tested on one map with one merge (section 11); whether
  the cap also costs the curves anything at 1M was not.
- **Times** are single runs on one laptop, pinned to one core. The 1M numbers are one map
  each.

## 13. What is not done

- **3D hull shells:** out of scope. In 3D the convex hull of a layer is a surface, not a
  loop. Its cities have no free exact order, so a shell is not an exact tour piece as a ring
  is, and each shell would need its own tour (a different approach). The number of shells in
  3D should grow more slowly than N^(2/3).
- **Strip at 1M** with M2–M5, and strip with R1 (section 2).
- **A 2D tight clustered family:** 2D has none. The 3D study's CT3d was not ported to 2D.
- **Other ring merges:** for example merging each ring into its outer neighbour with
  Or-opt-style insertion, or a cleanup that repairs long reversals better (Or-3opt, Lin-Kernighan).
- **Worst-case analysis:** none. This is an empirical study.

## 14. Reproduce

```
cd solvers && make 2d            # boxtsp2d, ringtest2d
bin/ringtest2d.exe               # 16 tests
py tests/regress_2d.py           # every earlier 2D tour unchanged
py bench/gen_shapes2d.py         # annulus, circles, blob, rectangle at 10k
py bench/hull_plan.py            # bench/hull_{main,big,r1,shapes}.json
py bench/import_runs.py --from results/2d/box2/runs.csv --to-config bench/hull_main.json --alg ...
py bench/run.py --config bench/hull_main.json     # then big, r1, shapes
py bench/hull_report.py          # results/2d/hull/summary.md, summary.csv, hull_clean_vs_raw.csv
py bench/hull_charts.py          # results/2d/hull_clean_vs_raw.png, hull_shapes.png
py bench/plot_curves.py          # results/2d/curves/*.png, including hull.png
```

Every `run.py` call is resumable and skips finished runs. One hull run by hand:
`bin/boxtsp2d.exe MAP TOUR --split rings --merge H0 --check`.
