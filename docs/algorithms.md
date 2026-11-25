# Algorithm descriptions

DIMACS rules: at most two pages in total. One section per algorithm that appears
in `results/2d/runs.csv`. For each one, say what it does, which data structures it
uses, its parameters and seeds, and which **speed tricks it uses and which it
skips**, because readers compare running times.

All solvers are compiled with the same flags (`solvers/Makefile`, `CFLAGS = -O2`)
and read instances with the DIMACS reader (`codes/Instance.c`).

---

## greedy (DIMACS baseline, `codes/greedy.exe`)

The DIMACS reference implementation of greedy (multi-fragment) construction.
Candidate edges sit in a min-heap. It repeatedly adds the shortest edge that keeps
every city at degree ≤ 2 and closes no cycle early (fragment ends are tracked in
`tail[]`). Nearest-neighbour queries use a k-d tree for geometric instances and a
linear scan for matrix instances. Deterministic. Included unchanged as a baseline.

## nn: nearest neighbour (`solvers/common/nn.c`)

- **Idea:** start at city 1 and repeatedly move to the closest unvisited city.
  The last city connects back to city 1.
- **Data structures:** a `visited[]` flag array and the tour array.
- **Complexity:** O(N²) distance evaluations.
- **Parameters / seeds:** none. Deterministic, and the start city is always 1.
- **Speed tricks used:** none.
- **Speed tricks skipped:** no k-d tree or neighbour lists for the nearest-unvisited
  search (these would give roughly O(N log N) on geometric instances), and no
  caching of distances.
- **Extra output:** `read_time`, `solve_time` (seconds, from `clock()`).

## spacefill: space-filling curve (`solvers/d2/spacefill.c`, `hilbert.c`)

- **Idea:** scale the cities onto a 65536 × 65536 grid and visit them in the order
  a Hilbert curve passes through their cells (Bartholdi & Platzman's method, using
  the Hilbert curve instead of their Sierpinski curve).
- **Data structures:** one (curve index, city) pair per city, sorted with `qsort`.
- **Complexity:** O(N log N). It makes no distance evaluations at all.
- **Parameters / seeds:** curve order 16. Deterministic.
- **Restrictions:** geometric instances only. Matrix instances have no coordinates
  and are skipped (`geometric_only` in `bench/config.json`).
- **Speed tricks skipped:** none relevant. The sort dominates.

## twoopt: 2-opt local search (`solvers/common/twoopt.c`, `localsearch.c`)

- **Idea:** repeatedly remove two tour edges and reconnect the tour the other way
  (reversing the path between them) whenever that shortens it. Stop at a local
  optimum.
- **Starting tour:** greedy (multi-fragment), built in rounds from candidate edges
  (`greedystart.c`). In round 1 the candidates are each city's 10 nearest
  neighbours. In each later round they are each fragment end's 8 nearest other
  fragment ends, recomputed on a grid. Candidates are sorted by length and kept
  when degrees stay ≤ 2 and no early cycle forms (union-find). Whatever fragments
  remain are chained nearest-end-first. The result is within ±2% of the DIMACS
  `greedy` tour length on the testbed (dsj1000: +6.6%). A space-filling-curve start
  was tried and rejected: 2-opt then ended about 18% over HK on E10k.0, versus
  about 5% from greedy.
- **Data structures:** tour array plus position array; 16-nearest-neighbour lists
  built with a uniform grid (geometric) or a row scan (matrix).
- **Speed tricks used:** neighbour lists (a move must add an edge to one of the 16
  nearest neighbours, using the standard gain criterion); a FIFO work queue (the
  equivalent of don't-look bits), where a city is re-examined only after its tour
  edges change; reversal of the shorter side of the cycle; best improving move per
  city; a reversal cap that skips moves reversing more than 100,000 cities. This
  only affects instances above 200k cities, and keeps array reversals from making
  million-city searches roughly O(N²). Measured on E1M.0: with the cap, search
  takes 3.4 s (2-opt) and 5.9 s (+Or-opt) at 5.87% / 4.86% over HK; without it,
  18.4 s and 33.5 s at 5.01% / 4.06%.
- **Speed tricks skipped:** 2-level doubly-linked lists or splay trees for the tour
  (array reversal costs O(N) in the worst case); quadrant neighbours; no final full
  sweep to catch moves the work queue missed.
- **Parameters / seeds:** 16 neighbours, 10 greedy candidates. Deterministic.
- **Extra output:** `read_time`, `neighbor_time`, `start_time`, `start_length`,
  `search_time`, `moves_2opt`.

## twoopt_or: 2-opt + Or-opt (`solvers/common/twoopt_or.c`, `localsearch.c`)

- **Idea:** same as `twoopt`, but when a city has no improving 2-opt move, also try
  Or-opt. Or-opt cuts out a chain of 1–3 consecutive cities starting at that city
  (in either tour direction), closes the gap, and re-inserts the chain, either way
  round, next to one of a chain end's 16 nearest neighbours. Each Or-opt move is
  applied as 1–3 2-opt exchanges.
- **Starting tour, data structures, speed tricks, parameters:** as `twoopt`.
- **Extra output:** as `twoopt`, plus `moves_oropt`.

## boxtsp: box start tour + cleanup (`solvers/d2/boxtsp.c`, `boxes.c`)

A Karp-style partitioning heuristic (Karp 1977), with no worst-case guarantee. One
program provides all variants of the box experiment, so they run identical code.

- **Idea:** sort the cities along a Hilbert curve on a 2²⁰ × 2²⁰ grid. Cut the
  sorted list into quadtree boxes of at most *t* cities (each box is a contiguous
  run of the sorted list). Solve each box exactly with Held-Karp dynamic
  programming. Glue the box loops together in curve order: each loop is opened at
  the edge that minimises (edge in − loop edge cut + look-ahead to the next box).
  Then clean up.
- **Cleanup:** `2opt` runs 2-opt rounds until a whole round with every city
  re-examined makes no move. `oropt` does that, then rounds of 2-opt + Or-opt
  until neither improves. Both use the same neighbour lists, gain criterion and
  reversal cap as `twoopt_or`.
- **Data structures:** sorted city list plus box start positions (no explicit
  tree); DP table of 2^(t−1) × (t−1); the bench tour array and neighbour lists.
- **Complexity:** sort O(N log N); boxes O(N · depth); DP O(N · 2^t · t); glue
  O(N · t).
- **Depth cap:** a box still over *t* cities at depth 20 (only possible when many
  cities share a grid cell: duplicates) keeps its curve order instead of being
  solved exactly.
- **Parameters:** `--t` (default 10, max 12), `--k` neighbours (default 16, the
  bench's value), `--start boxes|greedy|hilbert|curve`, `--cleanup none|2opt|oropt`,
  `--curve hilbert|moore|sierpinski|knopp|peano|morton|strip` (`curves2d.c`) and
  `--merge M1..M5` (`merges.c`; M1 = the original snake glue). The defaults
  (Hilbert, M1) give exactly the original tours. Deterministic. Comparison of
  the curves and merges: `docs/box_tsp_2d.md`.
- **Speed tricks skipped:** as `twoopt`. The DP table is not reused across boxes of
  different sizes. No random shift of the quadtree.
- **Extra output:** per-stage times (`sort_time`, `split_time`, `dp_time`,
  `glue_time`, `neighbor_time`, `twoopt_time`, `oropt_time`), `boxes`,
  `capped_boxes`, `largest_box`, `loops_sum`, `glue_G`, `glue_share`,
  `start_length`, `after_2opt_length`, move and round counts.
