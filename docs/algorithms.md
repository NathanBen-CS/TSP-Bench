# Algorithm descriptions

DIMACS rules: at most two pages in total. One section per algorithm that appears
in `results/runs.csv`. For each one, say what it does, which data structures it
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

## nn: nearest neighbour (`solvers/nn.c`)

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

## <next algorithm>

- **Idea:**
- **Data structures:**
- **Complexity:**
- **Parameters / seeds:**
- **Starting tour** (local search methods):
- **Speed tricks used:**
- **Speed tricks skipped:**
