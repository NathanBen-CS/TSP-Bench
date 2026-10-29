# Real-world maps: which algorithm is best for each type of problem

The 34 maps in `big/` are the TSPLIB maps used by the DIMACS challenge. Unlike
the generated E/C maps, they come from real applications. This page ranks every
algorithm in the bench on them, separately for each type of problem.

- **Chart:** `results/real_world.png` (top 5 per type)
- **Full rankings and per-map numbers:** `results/real_world.md`
- **Script:** `bench/real_world.py`. It only reads `results/runs.csv`; no new runs
  were needed (section 5).

Codes: a bare name = that order only; B = boxes (exact loops + merge); T = 2-opt;
O = Or-opt; R1 = one round of cleanup; Rn = rounds repeated until nothing
improves. So "greedy TO R1" is the bench's one-round `twoopt_or` solver and
"greedy TO Rn" is `boxtsp --start greedy --cleanup oropt`; every box setting is
Rn. The difference is explained in `docs/box_tsp_2d.md`, section 2.

## 1. What the maps are

Taken from the `COMMENT` line in each file's own header. rl5915's header only
says "5915-city TSP"; the University of Waterloo TSP site describes it as a
circuit-board drilling problem.

| Type | Maps | What the data is |
|---|---|---|
| **Cities and towns** (6) | usa13509, brd14051, d15112, d18512, fnl4461, nrw1379 | US towns with at least 500 people; West Germany (1989 borders); Germany; Germany including the former East Germany; the five former East German states; 1379 towns in North Rhine-Westphalia |
| **Circuit-board drilling** (14) | d1291, d1655, d2103, fl1400, fl1577, fl3795, pcb1173, pcb3038, u1060, u1432, u1817, u2152, u2319, rl5915 | holes a machine drills in a circuit board (Reinelt; Juenger/Reinelt) |
| **Chip layout** (3) | pla7397, pla33810, pla85900 | programmed logic array layouts (Johnson) |
| **Origin not stated** (9) | rl1304, rl1323, rl1889, rl5934, rl11849, pr1002, pr2392, vm1084, vm1748 | the headers only say "N-city problem". The rl maps may be drilling like rl5915, but no source says so, so they are not labelled. |
| **Not real data** (1) | dsj1000 | randomly generated clustered points (Johnson) |
| **Distance table only** (1) | si1032 | a 1032 × 1032 table of distances with no coordinates (Hofmeister) |

## 2. How the ranking works

- **Measure:** % above the optimal tour. All 34 TSPLIB maps have a proven optimal
  tour length (`refs/bounds`, third column), so 0% would be a perfect tour.
- **Rank:** the mean of that % over the maps of the type. Lower = better.
- **"Best on":** on how many of the type's maps the algorithm had the shortest
  tour of all 14.
- **All 14 algorithms** are ranked, nn included (these maps are at most 85,900
  cities, so nn finishes). On si1032 only greedy, nn, twoopt and twoopt_or can
  run, because the others need coordinates.

## 3. Top 5 for each type

**Cities and towns** (6 maps)

| Rank | Algorithm | Mean % above optimal | Best on |
|---|---|---|---|
| 1 | greedy TO R1 | 3.55 | 4 of 6 |
| 2 | greedy TO Rn | 3.57 | 2 of 6 |
| 3 | greedy T R1 | 4.67 | 0 |
| 4 | sierpinski BTO Rn M5 t=4 | 6.58 | 0 |
| 5 | sierpinski BTO Rn M1 t=10 | 7.05 | 0 |

**Circuit-board drilling** (14 maps)

| Rank | Algorithm | Mean % above optimal | Best on |
|---|---|---|---|
| 1 | greedy TO Rn | 5.54 | 8 of 14 |
| 2 | greedy TO R1 | 5.59 | 4 of 14 |
| 3 | greedy T R1 | 6.75 | 0 |
| 4 | sierpinski BTO Rn M5 t=4 | 9.01 | 0 |
| 5 | hilbert BTO Rn M5 t=10 | 9.04 | 1 of 14 |

**Chip layout** (3 maps)

| Rank | Algorithm | Mean % above optimal | Best on |
|---|---|---|---|
| 1 | greedy TO Rn | 4.79 | 3 of 3 |
| 2 | greedy TO R1 | 4.87 | 0 |
| 3 | greedy T R1 | 5.87 | 0 |
| 4 | hilbert BTO Rn M5 t=10 | 6.37 | 0 |
| 5 | sierpinski BTO Rn M5 t=4 | 6.75 | 0 |

**Origin not stated** (9 maps)

| Rank | Algorithm | Mean % above optimal | Best on |
|---|---|---|---|
| 1 | greedy TO Rn | 5.87 | 5 of 9 |
| 2 | greedy TO R1 | 6.65 | 2 of 9 |
| 3 | greedy T R1 | 8.02 | 0 |
| 4 | sierpinski BTO Rn M5 t=4 | 8.84 | 0 |
| 5 | hilbert BTO Rn M5 t=10 | 8.91 | 1 of 9 |

**Set apart** (one map each, so a rank can turn on one lucky or unlucky run):
- **dsj1000** (generated): 1. greedy TO R1 4.23, 2. sierpinski BTO Rn M1
  6.32, 3. knopp BTO Rn M1 6.67, 4. sierpinski BTO Rn M2 6.98, 5. sierpinski BTO Rn M5 7.52.
  greedy TO Rn is 7th (8.19): it starts from the same greedy tour as twoopt_or,
  but cleans it up in a different order and gets stuck in a worse place.
- **si1032** (distance table): 1. twoopt_or 0.41, 2. twoopt 0.43, 3. greedy 0.47,
  4. nn 2.07. Every solver is within 2.1% of optimal here, so this map says
  little about the differences between them.

## 4. What it shows

- **The same three solvers take places 1–3 on every real-world type:** greedy TO
  R1, greedy TO Rn and greedy T R1. greedy TO Rn is first on drilling, chip layout
  and the unstated maps; twoopt_or is first on cities, by 0.02 points.
- **The best box settings are 4th and 5th everywhere, and both use M5:** sierpinski
  M5 t=4 and hilbert M5 t=10 (on cities, sierpinski M1 takes 5th). They end
  3.0 (cities), 3.5 (drilling), 1.6 (chip) and 3.0 (unstated) points behind the
  winner.
- **M5 vs M1 matters most on drilling and chip maps**, as on the generated
  clustered maps: hilbert M5 against hilbert M1 is 9.0 vs 12.2% on drilling and
  6.4 vs 11.1% on chip layout, but 7.4 vs 7.4% on cities. Drilling and chip maps
  are made of tight rows and blocks of points, much like the clustered maps.
- **A box setting does beat greedy TO Rn on 4 of the 32 real-world maps** with
  coordinates: d1655 (hilbert M5, 6.6 vs 7.9%), fl1400 (sierpinski M1, 6.6 vs
  8.3%), rl1889 (hilbert M5, 8.3 vs 10.1%) and pr1002 (knopp M1, 7.2 vs 8.0%).
  On average, it never does.
- **Without cleanup the order is the same on every type:** greedy (15.7–17.8%),
  then nn (17.6–24.3%), then the plain curve orders (32.7–55.8%). Cleanup matters far more than
  which start is used.
- **Time doesn't decide anything here:** almost every solver takes well under a
  second on these maps. The exceptions are nn on the chip maps (14 s mean,
  because it is O(N²) and pla85900 has 85,900 cities) and sierpinski M2 t=12
  (1.7 s).

## 5. Caveats

- **One run per map.** Every solver is deterministic, so a repeat gives the
  same tour. But a type's mean is still a small sample (3–14 maps), and the
  rankings between greedy TO R1 and greedy TO Rn (0.02–0.8 points apart) and within the box
  settings (often under 0.5 points apart) would not survive other maps reliably.
- **Times are from the main bench run** (one run, core 0). Most are at the
  resolution limit of the timer.
- **Same files as before:** the bench's copies in `C:\Users\Natha\tsp-instances`
  are byte-identical to `big/*.tsp.gz` (checked for all 34).

## 6. Reproduce

```powershell
py bench\real_world.py        # writes results\real_world.md and results\real_world.png
```
