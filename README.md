# TSP benchmark: curves, boxes and local search, in 2D and 3D

| Folder | What is in it |
|---|---|
| `codes/` | the DIMACS programs (reader, generators, `greedy`, `length`), unchanged, plus the 3D generators `portgen3d.c` and `portcgen3d.c` |
| `solvers/common/` | everything that only uses `distance()` and the neighbour lists: tour, local search, greedy start, boxes, Held-Karp, merges, diagnostics. Compiled once, linked into both dimensions |
| `solvers/d2/` | 2D only: curves, Hilbert, grid neighbours, `boxtsp`, `spacefill` |
| `solvers/d3/` | 3D only: EUC_3D reader, kd-tree neighbours, curves, kd boxes, `boxtsp3d` |
| `bench/` | runner, plans, reports, charts |
| `tests/` | 2D regression, box tests (2D and 3D), 3D unit tests, Waterloo validation, link check |
| `docs/` | [algorithms](docs/algorithms.md), [2D boxes](docs/box_tsp.md), [2D curves and merges](docs/box_tsp_2d.md), [real-world maps](docs/real_world.md), [3D](docs/box_tsp_3d.md) |
| `results/2d/`, `results/3d/` | runs, summaries, charts |
| `refs/` | published bounds (`bounds`, `star_bounds.csv`) and our own Held-Karp bounds (`hk_own.csv`) |
| `data3d/` | 3D instances (not tracked; rebuilt by `bench/gen_instances3d.py` and `bench/get_star.py`) |
| `build/`, `bin/` | object files and programs, `2d` / `3d` suffix (not tracked) |
| `_to_delete/` | old build output moved out of `solvers/` |

```powershell
$env:PATH = "C:\msys64\ucrt64\bin;C:\msys64\usr\bin;" + $env:PATH
cd solvers; make; cd ..          # both dimensions; make 2d or make 3d for one
py tests\regress_2d.py           # every 2D tour identical to before the 3D work
bin\boxtest2d.exe; bin\unit3d.exe; bin\boxtest3d.exe; py tests\validate_star.py
```
