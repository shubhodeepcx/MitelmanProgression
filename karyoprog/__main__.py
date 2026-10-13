# Author: Shubhodeep Chatterjee
"""Command-line entry point: python -m karyoprog (or `karyoprog` after pip install).

    karyoprog run                 every step, 1 to 11
    karyoprog run --from 6        steps 6 to 11
    karyoprog run --only 8 AML    one step, with arguments passed to the script
    karyoprog steps               list the steps
    karyoprog report              rebuild report/main.tex and the DOCX
"""
import argparse
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STEPS = {
    1: ("01_build_clone_features.py", "per-clone features and case table"),
    2: ("02_power_threshold.py", "minimum event count N per disease"),
    3: ("03_freeze_vocabulary.py", "freeze the event vocabulary (run only to re-freeze)"),
    4: ("04_build_tables.py", "released event table, ordered pairs, matrices, case table"),
    5: ("05_descriptive.py", "selection-bias, era and data-quality tables"),
    6: ("06_cooccurrence.py", "co-occurrence baseline"),
    7: ("07_cbn.py", "conjunctive Bayesian networks"),
    8: ("08_mhn.py", "Mutual Hazard Networks"),
    9: ("09_clone_timing.py", "Bradley-Terry timing from clone orders"),
    10: ("10_validation.py", "validation tables and verdicts"),
    11: ("11_figures.py", "figures"),
}


def run_step(n, extra=()):
    script = ROOT / "scripts" / STEPS[n][0]
    print(f"[{n:>2}] {STEPS[n][1]}", flush=True)
    env = dict(os.environ, OMP_NUM_THREADS=os.environ.get("OMP_NUM_THREADS", "1"))
    subprocess.run([sys.executable, str(script), *extra], check=True, env=env)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="karyoprog", description="Per-disease cytogenetic progression models from Mitelman karyotypes.")
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="run pipeline steps")
    r.add_argument("--from", dest="start", type=int, default=1)
    r.add_argument("--to", dest="end", type=int, default=max(STEPS))
    r.add_argument("--only", type=int)
    r.add_argument("--refreeze", action="store_true", help="include step 3 (re-freeze the vocabulary)")
    r.add_argument("args", nargs="*", help="extra arguments for --only (e.g. disease names)")
    sub.add_parser("steps", help="list the steps")
    sub.add_parser("report", help="rebuild the report sources")
    a = ap.parse_args(argv)
    if a.cmd == "steps":
        for n, (f, what) in STEPS.items():
            print(f"{n:>2}  {f:<28} {what}")
    elif a.cmd == "report":
        subprocess.run([sys.executable, str(ROOT / "report" / "build_report.py")], check=True)
    elif a.only:
        run_step(a.only, a.args)
    else:
        for n in range(a.start, a.end + 1):
            if n == 3 and not a.refreeze:
                continue
            run_step(n)


if __name__ == "__main__":
    main()
