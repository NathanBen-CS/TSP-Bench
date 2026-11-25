"""Measure Waterloo's published star tours with our 3D reader and distance.

For every star instance: turn the published TSPLIB .tour file into a plain
city list, measure it with bin/length3d (the DIMACS tour checker linked with
the 3D reader), and compare with the published length (refs/star_bounds.csv,
in integer units of 1/10 parsec).  All must agree exactly: that proves the
reader and the rounding.

    py tests/validate_star.py       # exit 0 = every length reproduced
"""

import csv
import re
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
STAR = ROOT / "data3d" / "star"


def main():
    ok = True
    for r in csv.DictReader(open(ROOT / "refs" / "star_bounds.csv", newline="")):
        tsp, tour = STAR / r["file"], STAR / r["tour_file"]
        if not (tsp.exists() and tour.exists()):
            print(f"MISSING  {r['instance']}")
            ok = False
            continue
        text = tour.read_text()
        body = text.split("TOUR_SECTION", 1)[1].split("-1", 1)[0].split()
        claimed = re.search(r"Length\s*=?\s*(\d+)", text)
        plain = Path(tempfile.gettempdir()) / f"{r['instance']}.plain"
        plain.write_text("\n".join(body) + "\n")
        out = subprocess.run([str(ROOT / "bin" / "length3d.exe"), str(tsp), str(plain)],
                             capture_output=True, text=True)
        plain.unlink()
        m = re.search(r"Tour length is (\d+)", out.stdout)
        got = int(m.group(1)) if m else None
        good = got == int(r["best_tour"]) and len(body) == int(r["n"])
        ok &= good
        print(f"{'PASS' if good else 'FAIL'}  {r['instance']:10} measured {got}  published {r['best_tour']}"
              f"  (comment in the tour file: {claimed.group(1) if claimed else '-'}) {out.stderr.strip()}")
    print("Waterloo tours:", "every published length reproduced" if ok else "NOT all reproduced")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
