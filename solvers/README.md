# How the solvers work

This folder holds the TSP solvers used by the benchmark. This README explains
**what each algorithm does, how it does it, and why each detail is there**.
Every section starts with the simple idea and then goes deeper, so you can stop
at whatever level you need.

You don't need to know C to follow it. Code is quoted only where it helps.

**Where the code is:** `common/` holds everything that only needs `distance()` and the
neighbour lists; it is compiled once and linked into both the 2D and the 3D
programs. `d2/` and `d3/` hold the parts that need coordinates (reader, neighbour
search, curves) and the programs that differ. `make` in this folder builds all of
it into `bin/`, with a `2d` or `3d` suffix (`twoopt` below is `bin/twoopt2d.exe`,
and so on). This README describes the 2D programs; the 3D differences are in
[`docs/box_tsp_3d.md`](../docs/box_tsp_3d.md).

---

## Contents

1. [The problem in one minute](#1-the-problem-in-one-minute)
2. [How a tour is stored](#2-how-a-tour-is-stored)
3. [The pipeline: what each solver does, in order](#3-the-pipeline-what-each-solver-does-in-order)
4. [Space-filling curve (`spacefill`)](#4-space-filling-curve-spacefill)
5. [Neighbour lists: only look nearby](#5-neighbour-lists-only-look-nearby)
6. [Greedy start tour](#6-greedy-start-tour)
7. [2-opt (`twoopt`)](#7-2-opt-twoopt)
8. [Or-opt (`twoopt_or`)](#8-or-opt-twoopt_or)
9. [Making it fast on 10 million cities](#9-making-it-fast-on-10-million-cities)
10. [Results, and why they look the way they do](#10-results-and-why-they-look-the-way-they-do)
11. [Settings you can change](#11-settings-you-can-change)
12. [Glossary](#12-glossary)

---

## 1. The problem in one minute

You have N cities. Visit each exactly once and come back to the start, using
the **shortest total distance**. That round trip is called a **tour**.

The number of possible tours explodes: 10 cities have about 180,000 tours, and
1,000 cities have more than there are atoms in the universe. So nobody tries
them all. Instead we use **heuristics**: methods that find a *good* tour fast,
without promising the *best* one.

Every heuristic here makes one of two kinds of moves:

- **Construction** builds a tour from nothing (the space-filling curve, greedy).
- **Improvement** takes a finished tour and keeps making it shorter (2-opt, Or-opt).

The strong solvers do both: build a decent tour quickly, then polish it.

---

## 2. How a tour is stored

> Files: `common/tour.h`, `common/tour.c`

### The simple version

A tour is a list of city numbers in visiting order:

```
tour[] = [ 4, 1, 7, 3, 9, ... ]     → visit 4, then 1, then 7, ...
```

The trip back from the last city to the first is implied. Cities are numbered
**1 to N** (the DIMACS convention), but array slots start at **0**.

### One step deeper: the position array

The improvement algorithms constantly ask *"what comes after city 7?"*. With
only `tour[]` you'd have to search the whole list to find 7. So we keep a
second array that answers *"where is city c?"*:

```
pos[c] = the slot in tour[] where city c sits

tour[] = [ 4, 1, 7, 3, 9 ]
pos[4] = 0, pos[1] = 1, pos[7] = 2, pos[3] = 3, pos[9] = 4
```

Now "what comes after 7" is instant: `tour[pos[7] + 1]`, wrapping back to
slot 0 at the end. That is `succ(c)` (successor). `pred(c)` (predecessor) is the
same idea one slot back. The two arrays must always agree, so every function
that moves cities updates both.

### Deeper still: reversing a stretch of the tour

The key operation in this whole folder is **reversing a section of the tour**
(you'll see why in section 7). `reversePath(from, to)` swaps the two ends of the
section and walks inwards:

```
before:  1 [2 3 4 5] 6        reversePath(2, 5)
after:   1 [5 4 3 2] 6
```

One detail that matters a lot: a tour is a **circle**. Reversing the stretch
`2..5` gives exactly the same circle as reversing *everything else* (`6..1`); it
just goes round the other way. So the code always reverses **whichever side is
shorter**. On a 1,000,000-city tour, a reversal that would touch 900,000 cities
touches only 100,000 instead.

---

## 3. The pipeline: what each solver does, in order

| Step | `spacefill` | `twoopt` | `twoopt_or` |
|---|---|---|---|
| 1. Read the instance file | ✓ | ✓ | ✓ |
| 2. Build neighbour lists (16 nearest cities each) | | ✓ | ✓ |
| 3. Build a start tour | Hilbert curve **is** the answer | greedy | greedy |
| 4. Improve it | | 2-opt | 2-opt + Or-opt |
| 5. Write the tour file | ✓ | ✓ | ✓ |

Each solver prints its timings for each step (`read_time=`, `neighbor_time=`,
`start_time=`, `search_time=`), and the benchmark report shows them as extra columns.

Where the time goes on the biggest instance, E10M.0 (10 million cities), with
`twoopt_or`:

```
read file         ~3–10 s
neighbour lists     31 s
greedy start        51 s
2-opt + Or-opt      35 s
```

Construction costs more than improvement, which surprises most people. Section
9 explains why.

---

## 4. Space-filling curve (`spacefill`)

> Files: `d2/spacefill.c`, `d2/hilbert.c`

### The idea

Imagine drawing a single wiggly line that passes through **every** spot on a
map, never lifting the pen, in a way that keeps nearby spots close together
along the line. That is a **space-filling curve**. Visit the cities in the order
the line reaches them, and you have a tour.

```
A 4×4 Hilbert curve (starts bottom-left, ends bottom-right):

   5 ── 6    9 ── 10
   |    |    |    |
   4    7 ── 8    11
   |              |
   3 ── 2    13 ── 12
        |    |
   0 ── 1    14 ── 15
```

### How it works, step by step

1. **Put the map on a grid.** Every city is scaled onto a 65,536 × 65,536 grid
   (`ORDER 16`, so 2¹⁶ cells per side). Both axes use the same scale so the map
   isn't stretched.
2. **Find each city's position along the curve.** `hilbertIndex(x, y)` turns a
   cell into a number: 0 for the first cell the curve visits, 1 for the next,
   and so on.
3. **Sort the cities by that number.** The sorted order is the tour.

That's it. The method never computes a single distance.

### Deeper: how `hilbertIndex` finds the position

The Hilbert curve is built from quarters. The whole square splits into 4
quadrants, and the curve visits them in this order:

```
   1 │ 2          (y up)
  ───┼───
   0 │ 3
```

Each quadrant contains a smaller copy of the curve, rotated or flipped so the
pieces join up end to end. Inside each quadrant the same rule repeats: 4
sub-quadrants, and so on, 16 levels deep.

So to find a cell's position, the code reads the coordinates **one bit at a
time, from the biggest scale down**:

- The current bit of x and of y says which quadrant the cell is in:
  `(3·x_bit) XOR y_bit` gives 0, 1, 2 or 3 in the order shown above.
- That quadrant number × (cells per quadrant at this level) is added to the
  position.
- The coordinates are then rotated or flipped so the next, smaller level "sees"
  its quadrant the right way up.

Worked example on the 4×4 grid, cell (x=3, y=0), bottom-right:

```
level 1 (quadrants of size 2×2 = 4 cells):
  x bit = 1, y bit = 0  → quadrant 3 (bottom-right) → position += 3 × 4 = 12
  bottom-right quadrant is flipped: x,y → (0, 3), then swapped → (3, 0)
level 2 (single cells):
  x bit = 1, y bit = 0  → quadrant 3 → position += 3 × 1 = 3
total: 15   ✓ the last cell, matching the picture
```

### Why it's good, and why it isn't better

**Good:** the curve keeps neighbours together, so most steps in the tour are
short. Sorting N numbers takes O(N log N) time, so it handles 10 million cities in
about 20 seconds.

**Not better:** the curve has **seams**. Two cities can sit right next to each
other on the map but on opposite sides of a big-quadrant boundary, so the
curve visits one early and the other much later. The tour then has to jump
across the map and come back. That is why space-filling tours land about
**30–50% above the best possible length** (38% on E1M.0 and E10M.0).

**Restriction:** the method needs coordinates. Matrix instances (`M*`) only
list distances, so `spacefill` refuses them, and the benchmark marks those runs
`SKIPPED`.

---

## 5. Neighbour lists: only look nearby

> File: `d2/neighbors2d.c` (the 3D version is `d3/neighbors3d.c`)

### The idea

Good tours almost never connect a city to something far away. So instead of
considering all N−1 other cities every time, each city gets a short list of its
**16 nearest neighbours**, and the later steps only consider those. That turns
"check everything" (N work per city) into "check 16" (constant work per city).

### How the list is built: a grid

Comparing every city with every other city would take N² steps: 10¹⁴ at 10
million cities, which is far too slow. Instead:

1. **Cut the map into square cells**, sized so each holds about 2 cities.
2. **Drop each city into its cell.** This is a counting sort: count cities per
   cell, compute where each cell's group starts, then place them.
3. **For each city, search outward in rings:**

```
ring 0 = own cell      ring 1 = the 8 around it      ring 2 = the next 16 ...

          2 2 2 2 2
          2 1 1 1 2
          2 1 0 1 2
          2 1 1 1 2
          2 2 2 2 2
```

   Every city found goes into a sorted "best 16 so far" list (`offer()`, an
   insertion sort that keeps only the closest 16).

### Deeper: when is it safe to stop searching?

After finishing ring r, every city we *haven't* looked at is in ring r+1 or
further out. Our city sits somewhere inside its own cell, so anything in ring
r+1 is **at least r cell-widths away**. So:

> Once we have 16 candidates **and** the 16th is closer than r × cell width,
> nothing unseen can beat it, so we stop.

For evenly spread cities this stops after one or two rings, so building all
the lists costs roughly 16 × N work instead of N².

### Even deeper: the cache trick

Computers read memory fastest when consecutive reads are close together (the
**cache**). If we process cities in number order (1, 2, 3...), consecutive
cities are scattered randomly across the map, so each search reads a different
part of memory. Processing them **cell by cell** instead means consecutive
searches look at the same neighbouring cells, which are already in the cache.
That one change **halved** the build time on million-city instances.

### For matrix instances

There are no coordinates to grid, so each city simply scans its whole row of
the distance matrix. That's N² work, which is fine because the matrix instances
in the testbed have at most 10,000 cities.

---

## 6. Greedy start tour

> File: `common/greedystart.c`

### The idea

Look at all possible connections from **shortest to longest**, and keep each
one unless it breaks the rules:

1. **No city may get a third connection.** In a tour every city has exactly 2:
   one in, one out.
2. **No loop may close early.** Connecting two cities that are already linked
   through a chain would make a small circle that skips the other cities.

At first this creates many short chains (**fragments**). As longer
connections get added, the fragments merge, until one chain covers every
city. Connecting its two ends closes the tour.

```
early:     A─B    C─D─E    F    G─H        (fragments)
later:     A─B─C─D─E    F─G─H
end:       A─B─C─D─E─F─G─H   (+ H back to A)
```

### Deeper: "would this close a loop?" in almost no time

Rule 2 needs a fast answer to "are these two cities already in the same
fragment?". Walking along the fragment would be slow. The code uses a
**union-find** structure:

- Every city points to a "parent". Following parents leads to the fragment's
  **root**, which works like a name tag for the fragment.
- Two cities are in the same fragment exactly when they have the same root.
- Joining two fragments means pointing one root at the other.
- **Path halving:** while following parents, each city is re-pointed to its
  grandparent, so future lookups are shorter. In practice every lookup is
  almost instant.

### Deeper: which connections are considered?

Considering all N² possible connections is impossible at scale, so the work
happens in rounds:

| Round | Candidate connections |
|---|---|
| 1 | each city to its 10 nearest neighbours (from the neighbour lists) |
| 2, 3, … | each fragment **end** to its 8 nearest **other ends** (a fresh grid over just the ends) |
| final | whatever fragments remain get chained, always jumping to the nearest unused end |

After round 1 some cities are still fragment ends, because their nearby
cities were all "used up". Later rounds only look at those ends, so each round
is small. Each round merges most of the remaining fragments, so the final
chaining (which checks every remaining end) only handles a handful.

> **History:** the first version skipped the middle rounds and chained all
> leftover fragments directly. At 1M cities that took 25 s, because the chaining
> step's work grows with the *square* of the number of fragments. With the rounds
> it takes 4 s.

### Why start with greedy at all?

Section 7 explains how 2-opt improves a tour, but it can only make **small, local
changes**, so where it starts matters a lot. We measured it on E10k.0:

| Start tour | after 2-opt | after 2-opt + Or-opt |
|---|---|---|
| Space-filling curve | ~18% above best | ~8% above best |
| **Greedy** | **~5%** | **~4%** |

The space-filling tour's seams are long jumps across the map. Fixing a seam
would mean connecting cities that aren't in each other's 16-nearest lists, so
2-opt simply can't see the fix. Greedy tours are already built from short edges,
so the remaining problems are local, which is exactly what 2-opt is good at.

Our greedy comes within ±2% of the official DIMACS greedy program's tour
lengths, and is often slightly shorter.

---

## 7. 2-opt (`twoopt`)

> File: `common/localsearch.c` (the moves), `common/tour.c` (`exchange2`, `reversePath`)

### The idea

Look for two edges in the tour that would be better **swapped**:

```
before:   a ──── b              after:   a      b
                  ╲                       │      │
                   ╳          →           │      │
                  ╱                       │      │
          d ──── c                        c      d

removed: (a,b) and (c,d)          added: (a,c) and (b,d)
```

When two edges cross, swapping them **always** shortens the tour: the two
diagonals of a four-sided shape add up to more than either pair of opposite
sides.
2-opt is the systematic version of "untangle every crossing", and it also finds
improvements that aren't crossings.

### Why a swap is a reversal

The tour is a single loop. Remove (a,b) and (c,d) and you have two pieces. To
reconnect them with (a,c) and (b,d) **and keep one loop**, the piece between b
and c must be walked in the **opposite direction**:

```
before:   … a → b → x → y → c → d …
after:    … a → c → y → x → b → d …       (b..c reversed)
```

So every 2-opt move is: **reverse the path from b to c**. That's why
`reversePath` (section 2) is the core operation.

### Deeper: is a swap worth it? The gain formula

```
gain = (length of removed edges) − (length of added edges)
     = d(a,b) + d(c,d) − d(a,c) − d(b,d)
```

If the gain is positive, the tour gets shorter by exactly that amount. Distances
are whole numbers (DIMACS rounds them), so there are no rounding surprises.

### Deeper: finding swaps quickly, the "gain criterion"

For a city `a` and its tour neighbour `b`, we try each `c` from a's neighbour
list, **nearest first**, and compute

```
g = d(a,b) − d(a,c)        (what we gain by trading edge a–b for a–c)
```

If `g ≤ 0`, the new edge a–c is already at least as long as the old one. Every
later `c` in the list is even farther away, so `g` only gets worse, and we
**stop**. This is the standard *gain criterion* (Lin & Kernighan, 1973). It's
why a search usually checks only a few neighbours, not all 16.

Each city has two tour edges (to its successor and its predecessor), and both
get tried. Among all improving swaps found for that city, the **best** one is
applied.

### Deeper: `exchange2` and which way the array runs

After many reversals (and because we sometimes reverse the *other* side, see
section 2), the array might list the tour in either direction. `exchange2(x1,
x2, y1, y2)` means *"remove edges x1–x2 and y1–y2, add x1–y1 and x2–y2"*, and it
checks which way the array currently runs before choosing what to reverse. That
way the search code never has to care about direction.

### Deeper: the work queue (a.k.a. don't-look bits)

Naive 2-opt rescans every city after every change, which wastes almost all its
time: a swap only changes the tour around 4 cities.

So we keep a **queue of cities worth checking**:

1. At the start, every city is in the queue.
2. Take a city from the front. Search it.
3. **Improving swap found?** Apply it, and put its 4 endpoints (a, b, c, d) back
   in the queue, since their surroundings changed.
4. **Nothing found?** The city leaves the queue and is only re-added if a later
   swap changes its edges.
5. When the queue is empty, stop.

When it stops, the tour is a **local optimum**: no 2-opt swap among the
neighbour lists improves it. It isn't necessarily the best tour, just one that
this kind of move can't improve.

> **Honest caveat:** a city is re-queued only when *its own* edges change.
> Sometimes a neighbour's edge changes in a way that creates a new improving swap
> for a city that isn't re-queued. Checking every city's whole neighbourhood
> after every move would catch those, but it would be much slower. A brute-force
> check on pr1002 found 6 such leftover moves out of a 1,002-city tour, a tiny
> fraction. This trade-off is standard (Bentley 1992; Johnson & McGeoch 1997).

---

## 8. Or-opt (`twoopt_or`)

> File: `common/localsearch.c` (`orOptMove`)

### The idea

2-opt is bad at one particular fix: **a city (or a short run of cities) sitting
in the wrong place**. Picture a city that should be visited between two others
on the far side of its neighbourhood. Fixing that with 2-opt needs several
swaps in a row, and each one alone may make the tour *longer*, so 2-opt never
takes the first step.

Or-opt (Or, 1976) does it directly: **lift out a chain of 1–3 cities and put it
somewhere better**.

```
before:   a [s1 s2] b  ......  u  v
after:    a  b  ......  u [s1 s2] v          (or u [s2 s1] v, reversed)
```

### Deeper: is it worth it?

Two parts:

```
cut-out saving   g0 = d(a,s1) + d(s2,b) − d(a,b)
                 (remove the chain, close the gap with a–b)

insertion cost   d(u,s1) + d(s2,v) − d(u,v)        (same way round)
             or  d(u,s2) + d(s1,v) − d(u,v)        (reversed)

gain = g0 − insertion cost
```

The search tries chains of length 1, 2 and 3 starting at the current city,
extending in both tour directions. If `g0 ≤ 0` (cutting this chain out doesn't
even help), that chain is skipped. Insertion spots come from the neighbour
lists of the chain's two end cities, nearest first. The search stops early once
the new edge alone costs more than `g0` could ever save. Both orientations are
priced, and the best move overall is applied.

**Why only up to 3 cities?** Longer chains are rarely worth moving, and the
number of combinations to check grows quickly. Chains of 1–3 catch most
misplaced cities. This is the classic choice.

### Deeper: doing an Or-opt move with 2-opt swaps

Rather than writing separate "cut and paste" code (with its own bugs), an Or-opt
move is built from the **same `exchange2` swaps** that 2-opt uses. In the
general case:

```
start:     a [s1..s2] b ... u v
swap 1:    a  u ... b [s2..s1] v       remove (a,s1),(u,v)   add (a,u),(s1,v)
swap 2:    a  b ... u [s2..s1] v       remove (a,u),(b,s2)   add (a,b),(u,s2)
swap 3:    a  b ... u [s1..s2] v       remove (u,s2),(s1,v)  add (u,s1),(s2,v)
           (only if we want the chain the same way round)
```

Two special cases need only one swap before the optional third: when the
insertion spot is right next to the chain (`u` is `b`, or `v` is `a`). Since each
step is a valid 2-opt exchange, the tour stays a single valid loop throughout.

### How 2-opt and Or-opt share the work

For each city from the queue: **try 2-opt first. Only if that finds nothing,
try Or-opt.** 2-opt moves are cheaper to search, and they often fix things
before Or-opt is needed.

The result is consistently about **1 percentage point better** than 2-opt
alone (for example 4.2% vs 5.5% above the Held-Karp bound on E10k.0), for
roughly twice the search time (E1M.0: 5.9 s vs 3.4 s).

---

## 9. Making it fast on 10 million cities

Each fix below came from a real measurement.

| Problem | Fix | Effect |
|---|---|---|
| Checking all N cities for every move | **Neighbour lists** (16 per city) + **gain criterion** | each search checks a few cities, not N |
| Rescanning every city after each move | **Work queue** | only changed areas are re-examined |
| Reversals over most of the tour | **Reverse the shorter side** | each reversal ≤ N/2 |
| Fragment chaining took 25 s at 1M | **Join rounds** in greedy (section 6) | 25 s → 4 s |
| Neighbour search jumping around memory | **Process cities cell by cell** | build time halved |
| Search ~O(N²) on huge tours | **Reversal cap** (below) | E1M search: 18 s → 3 s |

### The reversal cap, in detail

With a greedy start, two cities that are close on the map are often far apart
in the tour list, because fragments got joined in arbitrary order. A swap
between them must reverse everything in between, which can be hundreds of
thousands of cities. On 1–10 million cities, these huge reversals dominate.

So `localsearch.c` **skips any move that would reverse more than
`MAX_REVERSAL = 100,000` cities**. Because we always reverse the shorter side
(at most half the tour), **instances up to 200,000 cities are completely
unaffected**. Their results are identical with or without the cap.

On bigger instances the cap trades a little quality for a lot of speed.
Measured on E1M.0 (1 million cities):

| Cap | 2-opt search | 2-opt gap | +Or-opt search | +Or-opt gap |
|---|---|---|---|---|
| 50,000 | 2.2 s | 6.16% | 3.9 s | 5.32% |
| **100,000 (used)** | **3.4 s** | **5.87%** | **5.9 s** | **4.86%** |
| 200,000 | 5.8 s | 5.62% | 12.2 s | 4.50% |
| none | 18.4 s | 5.01% | 33.5 s | 4.06% |

The cleaner long-term fix is a tour structure where reversal is cheap no matter
the length (a *2-level doubly-linked list*). That's a bigger project, and it's
listed as a speed trick we skipped in `docs/algorithms.md`.

---

## 10. Results, and why they look the way they do

"Gap" = how much longer than the best known tour (or above the Held-Karp lower
bound where no optimum is known).

| Instance | spacefill | twoopt | twoopt_or | twoopt_or time |
|---|---|---|---|---|
| pr1002 (real board) | 47.5% | 9.4% | 7.4% | < 0.1 s |
| E10k.0 (uniform random) | 39.4% | 5.5% | 4.2% | < 0.1 s |
| C10k.0 (clustered) | 52.5% | 11.3% | 10.4% | < 0.1 s |
| pla85900 (largest TSPLIB) | 32.5% | 5.5% | 4.6% | 0.7 s |
| E1M.0 | 38.2% | 5.9% | 4.9% | 14 s |
| E10M.0 | 38.2% | 6.2% | 5.4% | 130 s |

For reference, the DIMACS greedy program alone gets about 17%. Published
DIMACS results for neighbour-list 2-opt are about 5% on uniform instances
(Johnson & McGeoch), which matches ours.

**Why clustered instances (`C*`) are worse:** the 16 nearest neighbours of a city
inside a dense cluster are all in the same cluster. The edges *between*
clusters are what matter most, and they're never in anyone's 16-nearest list,
so 2-opt can't improve them.

**Why matrix instances (`M*`) are much worse (50–100%):** their distances are
random numbers with no geometry. "Close to A" and "close to B" tell you nothing
about A and B, so the triangle-inequality intuitions behind 2-opt break down.
Every heuristic struggles on them. Greedy alone is about 186% there.

**Why some long edges survive:** if fixing an edge needs a connection to a
city that isn't in anyone's 16-nearest list, neither 2-opt nor Or-opt can see
the move. Plot a tour with `bench/plot_tour.py` and you'll often spot one or
two.

---

## 11. Settings you can change

| Setting | Where | Default | What it controls |
|---|---|---|---|
| `NEIGHBORS` | `twoopt.c`, `twoopt_or.c` | 16 | neighbour list size. More = better tours, slower, more memory |
| `GREEDY_CAND` | `twoopt.c`, `twoopt_or.c` | 10 | candidate edges per city in greedy round 1 |
| `JOIN_K` | `greedystart.c` | 8 | candidate ends per end in greedy join rounds |
| `MAX_REVERSAL` | `localsearch.c` | 100000 | reversal cap (section 9). Can also be set with `make CFLAGS="-O2 -DMAX_REVERSAL=200000"` |
| `ORDER` | `hilbert.c` | 16 | Hilbert grid is 2^ORDER cells per side |

**Benchmark fairness:** the DIMACS rules require the same settings on every
instance. Changing a setting mid-benchmark makes earlier rows incomparable. If
you change one, either delete that solver's rows from `results/2d/runs.csv` or give
the new version a new name in `bench/config.json` (for example `twoopt_k20`).

---

## 12. Glossary

| Term | Meaning |
|---|---|
| **Tour** | a round trip visiting every city once |
| **Edge** | one leg of the tour, between two consecutive cities |
| **Heuristic** | a method that finds a good answer fast, without a guarantee it's the best |
| **Construction** | building a tour from scratch |
| **Local search / improvement** | repeatedly making small changes that shorten a tour |
| **Local optimum** | a tour that no single move of the chosen kind can improve |
| **Fragment** | a chain of connected cities during greedy construction |
| **Union-find** | a structure that answers "are these two in the same group?" almost instantly |
| **Neighbour list** | each city's k nearest cities, used to limit the search |
| **Gain criterion** | stop scanning neighbours once the new edge is no shorter than the one it replaces |
| **Don't-look bits / work queue** | only re-examine cities whose surroundings changed |
| **Held-Karp bound** | a mathematically proven lower limit: no tour can be shorter |
| **Gap** | (your length − reference) ÷ reference × 100 |
| **O(N²), O(N log N)** | how the work grows with N. O(N²): 10× the cities means 100× the work. O(N log N): about 10× |
| **Cache** | the CPU's small, very fast memory. Reading nearby data in a row is much faster than jumping around |

### Further reading

- J. L. Bentley, *Fast algorithms for geometric traveling salesman problems* (1992): neighbour lists, don't-look bits.
- D. S. Johnson & L. A. McGeoch, *The traveling salesman problem: a case study in local optimization* (1997): 2-opt and 3-opt in practice, reversal costs.
- J. J. Bartholdi & L. K. Platzman, *An O(N log N) planar travelling salesman heuristic based on spacefilling curves* (1982).
- I. Or, PhD thesis, Northwestern University (1976): Or-opt.
