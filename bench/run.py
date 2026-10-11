"""Run every algorithm on every instance and append one row per run to the
results CSV (one row per seed for randomized algorithms, never averaged).

For each run:
  1. launch the solver pinned to CPU 0:  solver <instance> <tour_file> [seed]
  2. record user CPU time and peak memory of the whole program
  3. check the tour with codes/length.exe and record its length
  4. append the row, delete the tour file

What you see while it runs:
  * a plan first: per algorithm, how many runs are new, already done, or
    skipped (and why)
  * one line per run with a counter [12/245] and a live elapsed timer, then
    the result: status, tour length, gap to the reference, time, memory
  * a summary at the end
Everything printed for finished runs is also appended to run.log next to
the results CSV.

Runs already in the CSV are skipped, so an interrupted run can be resumed:
just start it again.  Ctrl+C stops cleanly (the current run is discarded).

    py bench/run.py                       # everything in bench/config.json
    py bench/run.py --alg nn --inst E1k.0,pr1002
    py bench/run.py --list                # show the plan, run nothing
"""

import argparse
import csv
import datetime
import os
import sys
import tempfile

from common import (instance_is_geometric, list_instances, load_bounds, load_config,
                    parse_extras, parse_tour_length, pin_to_core0, resolve, run_measured)

COLUMNS = ["algorithm", "instance", "n", "seed", "tour_length", "user_s", "kernel_s",
           "wall_s", "peak_mb", "status", "extras", "timestamp"]

TTY = sys.stdout.isatty()          # live timer only makes sense in a terminal


def load_done(csv_path):
    if not csv_path.exists():
        return set()
    with open(csv_path, newline="") as f:
        return {(r["algorithm"], r["instance"], r["seed"]) for r in csv.DictReader(f)}


def append_row(csv_path, row):
    new = not csv_path.exists()
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with open(csv_path, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS)
        if new:
            w.writeheader()
        w.writerow(row)


def fmt_time(seconds):
    seconds = int(seconds)
    if seconds < 60:
        return f"{seconds}s"
    if seconds < 3600:
        return f"{seconds // 60}m{seconds % 60:02d}s"
    return f"{seconds // 3600}h{seconds % 3600 // 60:02d}m"


def gap_text(length, name, bounds):
    """'+5.42% vs opt' using the optimal/best-known length, else the HK bound."""
    hk, opt = bounds.get(name, (None, None))
    ref, label = (opt, "opt") if opt else (hk, "HK")
    if not ref or length in ("", None):
        return ""
    return f"{(int(length) - ref) / ref * 100:+.2f}% vs {label}"


def skip_reason(alg, path, n):
    if n > alg.get("max_n", float("inf")):
        return f"n > max_n {alg['max_n']}"
    if alg.get("geometric_only") and not instance_is_geometric(path):
        return "matrix instance (geometric only)"
    return None


def run_one(alg, inst_path, seed, length_exe, heartbeat):
    fd, tour = tempfile.mkstemp(prefix="tour_", suffix=".txt")
    os.close(fd)
    try:
        cmd = [resolve(alg["exe"]), inst_path, tour]
        if seed != "":
            cmd.append(seed)
        r = run_measured(cmd, timeout=alg.get("timeout_s"), heartbeat=heartbeat)
        row = {"user_s": f"{r['user_s']:.3f}", "kernel_s": f"{r['kernel_s']:.3f}",
               "wall_s": f"{r['wall_s']:.3f}", "peak_mb": f"{r['peak_mb']:.1f}",
               "tour_length": "",
               "extras": ";".join(f"{k}={v}" for k, v in parse_extras(r["stdout"]).items())}
        problem = ""
        if r["timed_out"]:
            row["status"] = "DNF"
            problem = f"stopped after timeout_s={alg.get('timeout_s')}"
        elif r["exit_code"] != 0:
            row["status"] = "FAILED"
            problem = (r["stderr"] or r["stdout"]).strip()[-300:]
        else:
            check = run_measured([length_exe, inst_path, tour])
            length = parse_tour_length(check["stdout"])
            if check["exit_code"] == 0 and length is not None:
                row["status"], row["tour_length"] = "OK", length
            else:
                row["status"] = "INVALID"
                problem = "length.exe: " + (check["stderr"] or check["stdout"]).strip()[-300:]
        return row, problem
    finally:
        os.remove(tour)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    ap.add_argument("--alg", help="comma-separated algorithm names (default: all)")
    ap.add_argument("--inst", help="comma-separated instance names (default: config)")
    ap.add_argument("--list", action="store_true",
                    help="print every run that would happen, then exit")
    args = ap.parse_args()

    cfg = load_config(args.config)
    instances = list_instances(resolve(cfg["instance_dir"]))
    wanted = args.inst.split(",") if args.inst else cfg.get("instances", "all")
    if wanted != "all":
        instances = [t for t in instances if t[0] in wanted]
    algs = cfg["algorithms"]
    if args.alg:
        algs = [a for a in algs if a["name"] in args.alg.split(",")]
    if not algs or not instances:
        sys.exit("Nothing selected: check --alg / --inst names and instance_dir "
                 f"({resolve(cfg['instance_dir'])}).")

    csv_path = resolve(cfg["results_csv"])
    log_path = csv_path.parent / "run.log"
    length_exe = resolve(cfg["length_exe"])
    bounds = load_bounds(resolve(cfg["bounds_file"]))
    done = load_done(csv_path)

    # ---- 1. plan: sort every (algorithm, instance, seed) into done/skip/todo
    todo, skips = [], []
    print(f"Results file: {csv_path}")
    print(f"Plan ({len(algs)} algorithm(s) x {len(instances)} instance(s)):")
    for alg in algs:
        seeds = [str(s) for s in alg.get("seeds", [])] if alg.get("randomized") else [""]
        n_done = n_new = 0
        reasons = {}
        for name, path, n in instances:
            for seed in seeds:
                if (alg["name"], name, seed) in done:
                    n_done += 1
                    continue
                why = skip_reason(alg, path, n)
                if why:
                    skips.append((alg, name, n, seed, why))
                    reasons[why] = reasons.get(why, 0) + 1
                else:
                    todo.append((alg, name, path, n, seed))
                    n_new += 1
        skipped = ", ".join(f"{c} skipped: {w}" for w, c in reasons.items())
        print(f"  {alg['name']:10} {n_new:4} to run   {n_done:4} already done"
              + (f"   {skipped}" if skipped else ""))

    if args.list:
        for alg, name, path, n, seed in todo:
            print(f"  would run: {alg['name']:10} {name:12} n={n}"
                  + (f" seed={seed}" if seed else ""))
        return

    if not todo and not skips:
        print("Nothing to do: every selected run is already in the results file.")
        return

    # skipped runs are recorded right away, so the report shows them
    for alg, name, n, seed, why in skips:
        append_row(csv_path, {"algorithm": alg["name"], "instance": name, "n": n,
                              "seed": seed, "status": "SKIPPED",
                              "timestamp": datetime.datetime.now().isoformat(timespec="seconds")})

    # ---- 2. run them, smallest instances first within each algorithm
    pin_to_core0()
    print(f"\nRunning {len(todo)} run(s), pinned to CPU 0.  Ctrl+C stops; "
          "run again to resume.\n")
    counts = {}
    started = datetime.datetime.now()
    with open(log_path, "a") as log:
        log.write(f"\n=== {started:%Y-%m-%d %H:%M:%S} start: {len(todo)} run(s)\n")
        try:
            for k, (alg, name, path, n, seed) in enumerate(todo, 1):
                label = (f"[{k}/{len(todo)}] {alg['name']:10} {name:12} n={n:<9}"
                         + (f" seed={seed}" if seed else ""))

                minutes_shown = [0]

                def heartbeat(elapsed, label=label, minutes_shown=minutes_shown):
                    if TTY:                                  # rewrite the same line
                        print(f"\r{label} running {fmt_time(elapsed)} ",
                              end="", flush=True)
                    elif int(elapsed // 60) > minutes_shown[0]:   # not a terminal:
                        minutes_shown[0] = int(elapsed // 60)     #   a line per minute
                        print(f"\n{label} still running ({fmt_time(elapsed)})",
                              end="", flush=True)

                print(f"{label} running", end="", flush=True)
                row, problem = run_one(alg, path, seed, length_exe, heartbeat)
                row.update(algorithm=alg["name"], instance=name, n=n, seed=seed,
                           timestamp=datetime.datetime.now().isoformat(timespec="seconds"))
                append_row(csv_path, row)

                result = (f"{label} {row['status']:7} {str(row['tour_length']):>12}  "
                          f"{gap_text(row['tour_length'], name, bounds):18} "
                          f"user {float(row['user_s']):8.2f}s  peak {row['peak_mb']}MB")
                print(("\r" if TTY else "\n") + result + " " * 10)
                log.write(result + "\n")
                if problem:
                    print(f"      {problem}")
                    log.write(f"      {problem}\n")
                log.flush()
                counts[row["status"]] = counts.get(row["status"], 0) + 1
        except KeyboardInterrupt:
            print("\n\nStopped by Ctrl+C. The run in progress was discarded; "
                  "everything before it is saved. Run the same command to resume.")
            log.write("stopped by Ctrl+C\n")
            sys.exit(130)

        # ---- 3. summary
        took = datetime.datetime.now() - started
        summary = ", ".join(f"{c} {s}" for s, c in sorted(counts.items()))
        print(f"\nDone in {fmt_time(took.total_seconds())}: {summary}"
              + (f"; {len(skips)} SKIPPED" if skips else ""))
        print(f"Results: {csv_path}\nLog:     {log_path}\n"
              "Next:    py bench\\report.py")
        log.write(f"=== done in {fmt_time(took.total_seconds())}: {summary}\n")


if __name__ == "__main__":
    main()
