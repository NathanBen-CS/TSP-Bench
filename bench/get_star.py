"""Download the Waterloo star instances and their published tours.

    py bench/get_star.py            # into data3d/star, skips files already there

Source pages: https://math.uwaterloo.ca/tsp/star/index.html (the data pages
link the .tsp files, the tour pages the .tour files).  The .tsp.gz files are
unpacked next to the download.  star10m is not downloadable: its data page
asks for an email to the site's owner (433 MB), so it is not used here.
"""

import gzip
import shutil
import sys
import urllib.request

from common import ROOT

BASE = "https://math.uwaterloo.ca/tsp/star/data/"
FILES = ["star100.tsp", "star100.tour", "star1k.tsp", "star1k.tour", "star10k.tsp", "star10k.tour",
         "kj37859.tsp", "kj37859.tour", "hyg109399.tsp", "hyg109399.tour", "hyg119614.tsp",
         "hyg119614.tour", "star250k.tsp.gz", "star250k.9794191.tour.gz", "gaia2079471.tsp.gz",
         "gaia2079471.288843524.tour.gz"]


def main():
    out = ROOT / "data3d" / "star"
    out.mkdir(parents=True, exist_ok=True)
    failed = []
    for f in FILES:
        p = out / f
        if not p.exists():
            print("  downloading", f, flush=True)
            try:
                with urllib.request.urlopen(BASE + f, timeout=120) as r, open(str(p) + ".part", "wb") as w:
                    shutil.copyfileobj(r, w)
                (out / (f + ".part")).replace(p)
            except OSError as e:
                failed.append(f)
                print("    FAILED:", e)
                continue
        if f.endswith(".gz") and not p.with_suffix("").exists():
            with gzip.open(p, "rb") as r, open(p.with_suffix(""), "wb") as w:
                shutil.copyfileobj(r, w)
    if failed:
        sys.exit("Not downloaded: " + ", ".join(failed))
    print("All star files are in", out)


if __name__ == "__main__":
    main()
