import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config.paths import *

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path


# ---------------------------------------------------------------------
# ROOT PATH (fixed & stable)
# ---------------------------------------------------------------------

ROOT = Path(__file__).resolve().parents[2]
FIG_DIR = Path(__file__).resolve().parent
OUT_DIR = ROOT / "outputs" / "figures"

OUT_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------
# FIGURES
# ---------------------------------------------------------------------

FIGURES = [
    "make_figure1.py",
    "make_figure2.py",
    "make_figure3.py",
    "make_figure4.py",
    "make_figure5.py",
    "make_figure6.py",
]


# ---------------------------------------------------------------------
# RUN SCRIPT
# ---------------------------------------------------------------------

def run(script_path: Path, i: int, n: int):

    print("\n" + "=" * 70)
    print(f"[{i}/{n}] Running {script_path.name}")
    print("=" * 70)

    subprocess.run(
        [sys.executable, str(script_path)],
        cwd=FIG_DIR,
        check=True
    )


# ---------------------------------------------------------------------
# UNIFY OUTPUTS (STRICT FILTER, NO RANDOM FILES)
# ---------------------------------------------------------------------

def unify_outputs():
    """
    Copy only official figure outputs into a single folder.
    Prevents grabbing unrelated or duplicate files.
    """

    allowed_keywords = (
        "Figure1",
        "Figure2",
        "Figure3",
        "Figure4",
        "Figure5",
        "Figure6",
    )

    patterns = ["*.png", "*.tiff", "*.pdf"]

    for pattern in patterns:
        for f in ROOT.rglob(pattern):

            name = f.name

            # keep only real figure outputs
            if not name.startswith(allowed_keywords):
                continue

            try:
                target = OUT_DIR / name

                # overwrite safely (latest version wins)
                target.write_bytes(f.read_bytes())

            except Exception:
                pass


# ---------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------

def main():

    print("\n" + "=" * 70)
    print("FINAL FIGURE PIPELINE (CLEAN + UNIFIED OUTPUT)")
    print("=" * 70)

    t0 = time.perf_counter()

    for i, f in enumerate(FIGURES, 1):

        script_path = FIG_DIR / f

        if not script_path.exists():
            raise FileNotFoundError(script_path)

        run(script_path, i, len(FIGURES))

    print("\nCollecting outputs into single folder...")

    unify_outputs()

    print("\n" + "=" * 70)
    print("DONE - ALL FIGURES IN ONE FOLDER")
    print(f"OUTPUT: {OUT_DIR}")
    print(f"TIME: {time.perf_counter() - t0:.2f}s")
    print("=" * 70)


if __name__ == "__main__":
    main()