"""Copy finished runs from one results file into another, instead of re-running.

A run can be reused when it was made with the same program and options on
the same map file.  This copies rows for the given algorithms (optionally
renamed) from a source results CSV into a target one, for instances the
target's config actually uses, skipping any (algorithm, instance, seed)
the target already has.  bench/run.py then only runs what is missing.

    py bench/import_runs.py --from results/2d/box/runs.csv --to-config bench/config.json \
        --alg "greedy+2opt+oropt" --alg "boxes+2opt+oropt=box:hilbert:M1:t10"

Each --alg is NAME or SOURCE_NAME=TARGET_NAME.  Imported rows keep their
original timestamp, so they can be told apart.  Reused runs are only valid
when the solver code has not changed since: the existing solvers and the
boxtsp defaults were checked byte-identical after every refactor.
"""

import argparse
import csv

from common import list_instances, load_config, resolve
from run import COLUMNS, append_row, load_done


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--from", dest="source", required=True, help="source results CSV")
    ap.add_argument("--to-config", required=True, help="config of the target results file")
    ap.add_argument("--alg", action="append", required=True,
                    help="NAME or SOURCE_NAME=TARGET_NAME (repeatable)")
    args = ap.parse_args()

    cfg = load_config(args.to_config)
    target = resolve(cfg["results_csv"])
    wanted = cfg.get("instances", "all")
    names = {n for n, _, _ in list_instances(resolve(cfg["instance_dir"]))}
    if wanted != "all":
        names &= set(wanted)
    rename = dict((a.split("=", 1) + [a])[:2] if "=" in a else (a, a) for a in args.alg)
    done = load_done(target)

    copied = skipped = 0
    with open(resolve(args.source), newline="") as f:
        for r in csv.DictReader(f):
            if r["algorithm"] not in rename or r["instance"] not in names or r["status"] != "OK":
                continue
            r = {k: r.get(k, "") for k in COLUMNS}
            r["algorithm"] = rename[r["algorithm"]]
            if (r["algorithm"], r["instance"], r["seed"]) in done:
                skipped += 1
                continue
            append_row(target, r)
            done.add((r["algorithm"], r["instance"], r["seed"]))
            copied += 1
    print(f"{target}: copied {copied} run(s), {skipped} already there")


if __name__ == "__main__":
    main()
