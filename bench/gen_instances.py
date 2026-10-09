"""Build the DIMACS testbed in the instance directory from bench/config.json.

  * random instances: every "portgen/portcgen/portmgen N seed > NAME" line
    in codes/README, run with the compiled generators in codes/
  * TSPLIB instances: big/*.tsp.gz decompressed
  * self-check: E1k.1, C1k.1, M1k.1 must equal sample/ (ignoring CR)

Existing files are skipped, so the script can be re-run safely.

    python bench/gen_instances.py              # full testbed (~1 GB)
    python bench/gen_instances.py --max-n 10000
"""

import argparse
import gzip
import re
import shutil
import subprocess
import sys

from common import ROOT, load_config, resolve

README_CMD = re.compile(r"^(portgen|portcgen|portmgen) (\d+) (\d+) > (\S+)\s*$")


def readme_commands():
    cmds = []
    with open(ROOT / "codes" / "README") as f:
        for line in f:
            m = README_CMD.match(line)
            if m:
                cmds.append((m.group(1), int(m.group(2)), int(m.group(3)), m.group(4)))
    return cmds


def generate(gen, n, seed, out_path):
    tmp = out_path.with_name(out_path.name + ".part")
    with open(tmp, "wb") as out:
        subprocess.run([str(ROOT / "codes" / f"{gen}.exe"), str(n), str(seed)],
                       stdout=out, check=True)
    tmp.replace(out_path)        # only complete files get the real name


def same_ignoring_cr(a, b):
    with open(a, "rb") as fa, open(b, "rb") as fb:
        return fa.read().replace(b"\r", b"") == fb.read().replace(b"\r", b"")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--max-n", type=int, default=None,
                    help="skip instances with more cities than this")
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    inst_dir = resolve(load_config(args.config)["instance_dir"])
    inst_dir.mkdir(parents=True, exist_ok=True)
    print(f"Instance directory: {inst_dir}")

    for gen, n, seed, name in readme_commands():
        if args.max_n and n > args.max_n:
            continue
        out = inst_dir / name
        if out.exists():
            continue
        print(f"  {gen} {n} {seed} > {name}", flush=True)
        generate(gen, n, seed, out)

    for gz in sorted((ROOT / "big").glob("*.tsp.gz")):
        out = inst_dir / gz.name[:-3]
        if out.exists():
            continue
        n = int(re.search(r"(\d+)$", out.stem).group(1))   # TSPLIB names end in N
        if args.max_n and n > args.max_n:
            continue
        print(f"  gunzip {gz.name}", flush=True)
        with gzip.open(gz, "rb") as src, open(out, "wb") as dst:
            shutil.copyfileobj(src, dst)

    ok = True
    for name in ("E1k.1", "C1k.1", "M1k.1"):
        path = inst_dir / name
        if not path.exists():
            continue
        good = same_ignoring_cr(path, ROOT / "sample" / name)
        ok &= good
        print(f"Self-check {name}: {'matches sample/' if good else 'MISMATCH with sample/'}")
    if not ok:
        sys.exit("Generator self-check failed.")


if __name__ == "__main__":
    main()
