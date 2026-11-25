"""Check that every path the docs and bench configs mention exists.

  * markdown links [text](path) in *.md (docs/, solvers/, results/, root)
  * back-quoted paths such as `solvers/common/tour.c` or `results/2d/runs.csv`
    (project-relative; a name with * is a glob; programs (.exe) may be unbuilt)
  * exe, length_exe, bounds_file and results_csv in bench/*.json

    py tests/check_links.py        # exit 0 = nothing broken
"""

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
bad = []
n = 0
mds = [p for p in ROOT.rglob("*.md") if "_to_delete" not in p.parts and ".git" not in p.parts]
for md in mds:
    text = md.read_text(encoding="utf-8")
    text_nocode = re.sub(r"```.*?```", "", text, flags=re.S)
    for link in re.findall(r"\]\(([^)\s]+)\)", text_nocode):
        if link.startswith(("http", "mailto", "#")):
            continue
        n += 1
        if not (md.parent / link.split("#")[0]).exists():
            bad.append(f"{md.relative_to(ROOT)}: link {link}")
    for p in re.findall(r"`((?:solvers|bench|results|docs|codes|tests|refs|bin|data3d)[/\][^`\s]+)`", text_nocode):
        q = p.replace("\\", "/").rstrip(".,")
        if "<" in q or "{" in q or q.startswith("data3d"):
            continue
        n += 1
        if "*" in q:
            ok = any(ROOT.glob(q))
        else:
            ok = (ROOT / q).exists() or (q.startswith("bin/") and not q.endswith(".exe")
                                         and (ROOT / (q + ".exe")).exists())
        if not ok:
            bad.append(f"{md.relative_to(ROOT)}: path {p}")
for cfg in (ROOT / "bench").glob("*.json"):
    c = json.loads(cfg.read_text(encoding="utf-8"))
    for p in [c["length_exe"], c["bounds_file"], c["results_csv"]] + [a["exe"] for a in c["algorithms"]]:
        n += 1
        if not (ROOT / p).exists() and not p.startswith("results/3d"):
            bad.append(f"{cfg.relative_to(ROOT)}: {p}")
for b in sorted(set(bad)):
    print("BROKEN ", b)
print(f"{n} paths checked, {len(set(bad))} broken")
sys.exit(1 if bad else 0)
