#!/usr/bin/env python3
"""Region TSV rerun, saved-value checks and Figures 1, 2 and 4 rendering."""
from pathlib import Path
import os
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
PROJECT = ROOT / "frozen_intermediates/region_composite_project"
ANALYSIS = ROOT / "nargab_region_composite/run_region_composite_v1.py"
RUN = ROOT / "reproduction_output"


def call(args, env=None):
    subprocess.run(args, check=True, env=env)


def main():
    av = RUN / "available_case"
    cc = RUN / "complete_case"
    av.mkdir(parents=True, exist_ok=True)
    cc.mkdir(parents=True, exist_ok=True)
    call([sys.executable, str(ANALYSIS), "--frozen-project", str(PROJECT), "--population-mode", "AVAILABLE_CASE", "--output", str(av)])
    call([sys.executable, str(ANALYSIS), "--frozen-project", str(PROJECT), "--population-mode", "COMPLETE_CASE", "--output", str(cc)])
    call([sys.executable, str(ROOT / "verify_nargab_outputs.py"), "--rerun-available", str(av), "--rerun-complete", str(cc)])
    env = os.environ.copy()
    env.setdefault("MPLCONFIGDIR", str(RUN / "matplotlib_cache"))
    call([sys.executable, str(ROOT / "figure_scripts/build_nargab_figures.py")], env=env)
    print(f"PASS: statistics and Figures 1, 2 and 4 rebuilt; Figure 3 is the unchanged frozen primary figure. Region TSV output: {RUN}; figure output: {ROOT.parent / 'Figures'}")


if __name__ == "__main__":
    main()
