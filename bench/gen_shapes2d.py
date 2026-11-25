"""Extra 2D map families for the hull-ring study (docs/hull_tsp.md, section 6).

All in the 10^6 square with integer coordinates, like portgen; TSPLIB EUC_2D
files with the same header style.  10,000 cities, seeds .0 .1 .2 (seed .0 is
N itself, .k is N * 10 + k, as in the DIMACS set).

  ann10k.k    annulus: uniform by area between radius 250,000 and 500,000
              around the centre
  circ10k.k   concentric circles: 20 circles, radii 25,000 ... 500,000 in equal
              steps, 500 cities each at uniform random angles
  blob10k.k   one Gaussian blob at the centre, standard deviation 100,000
              (cities outside the square are drawn again)
  rect10k.k   thin rectangle 10:1: x uniform in [0, 10^6), y in [0, 10^5)

The coordinates are rounded to integers, so the circle cities are only
approximately on their circles.

    py bench/gen_shapes2d.py            # writes to C:/Users/Natha/tsp-instances-shapes
"""

import math
import random
from pathlib import Path

OUT = Path("C:/Users/Natha/tsp-instances-shapes")
N = 10000
SIDE = 1000000
C = SIDE / 2


def annulus(rng):
    r = math.sqrt(rng.random() * (500000 ** 2 - 250000 ** 2) + 250000 ** 2)
    a = rng.random() * 2 * math.pi
    return C + r * math.cos(a), C + r * math.sin(a)


def circles(rng, i):
    r = 25000 * (1 + i % 20)                       # city i on circle i mod 20
    a = rng.random() * 2 * math.pi
    return C + r * math.cos(a), C + r * math.sin(a)


def blob(rng):
    while True:
        x, y = rng.gauss(C, 100000), rng.gauss(C, 100000)
        if 0 <= x < SIDE and 0 <= y < SIDE:
            return x, y


def rect(rng):
    return rng.random() * SIDE, rng.random() * SIDE / 10


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for fam, make in (("ann", lambda r, i: annulus(r)), ("circ", circles),
                      ("blob", lambda r, i: blob(r)), ("rect", lambda r, i: rect(r))):
        for k in range(3):
            name = f"{fam}10k.{k}"
            path = OUT / name
            if path.exists():
                continue
            seed = N if k == 0 else N * 10 + k
            rng = random.Random(f"{fam}-{seed}")
            lines = [f"NAME : {fam}-{N}-{seed}", f"COMMENT : gen_shapes2d.py {fam} N={N}, seed={seed}",
                     "TYPE : TSP", f"DIMENSION : {N}", "EDGE_WEIGHT_TYPE : EUC_2D", "NODE_COORD_SECTION"]
            for i in range(N):
                x, y = make(rng, i)
                lines.append(f"{i + 1} {int(round(x))} {int(round(y))}")
            lines.append("EOF")
            tmp = path.with_name(name + ".part")
            tmp.write_text("\n".join(lines) + "\n")
            tmp.replace(path)
            print(f"  {name}")


if __name__ == "__main__":
    main()
