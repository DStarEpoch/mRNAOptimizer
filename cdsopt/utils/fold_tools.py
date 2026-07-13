# -*- coding: utf-8 -*-
"""
RNA folding tools with auto-detection of LinearFold and ViennaRNA fallback.
"""
from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

from ViennaRNA import RNA


def _find_linearfold_executable() -> str | None:
    """Locate a LinearFold executable.

    Search order:
      1. ``linearfold`` on PATH (common on Linux/macOS).
      2. ``linearfold_v`` / ``linearfold_v.exe`` on PATH.
      3. Bundled Windows binary next to this file:
         ``mRNAOptimizer/submodules/LinearFold/bin/linearfold_v.exe``.
      4. Bundled ``linearfold.exe`` in the same directory.
    """
    for name in ("linearfold", "linearfold_v", "linearfold_v.exe"):
        exe = shutil.which(name)
        if exe:
            return exe

    # Fallback to bundled binary relative to this module.
    # This file lives in mRNAOptimizer/cdsopt/utils/fold_tools.py.
    bundled = (
        Path(__file__).resolve().parents[2]
        / "submodules"
        / "LinearFold"
        / "bin"
        / "linearfold_v.exe"
    )
    if bundled.exists():
        return str(bundled)

    bundled_alt = bundled.with_name("linearfold.exe")
    if bundled_alt.exists():
        return str(bundled_alt)

    return None


def _fold_vienna(sequence: str, need_aup: bool = True) -> dict:
    """Fold with ViennaRNA; AUP from partition-function pairing probs."""
    fc = RNA.fold_compound(sequence)
    ss, mfe = fc.mfe()
    if need_aup:
        aup = _compute_aup_from_structure(ss)
    else:
        aup = 0.0
    return {"mfe": mfe, "structure": ss, "aup": aup}


def _compute_aup_from_structure(structure: str) -> float:
    """Compute AUP from a dot-bracket structure (MFE-based, consistent across engines)."""
    if not structure:
        return 0.0
    return structure.count(".") / len(structure)


def _fold_linearfold(sequence: str, need_aup: bool = True) -> dict:
    """Fold with LinearFold CLI; AUP from dot-counting."""
    executable = _find_linearfold_executable()
    if executable is None:
        raise FileNotFoundError("linearfold not found in PATH or bundled location")

    # The Windows build of linearfold_v.exe reads stdin and prints:
    #   STRUCTURE  (MFE)
    # e.g.  (((...)))  (-1.20)
    if sys.platform == "win32":
        result = subprocess.run(
            [executable],
            input=sequence + "\n",
            capture_output=True,
            text=True,
        )
    else:
        result = subprocess.run(
            [executable, "-V"],
            input=sequence + "\n",
            capture_output=True,
            text=True,
        )

    if result.returncode != 0:
        raise RuntimeError(
            f"LinearFold failed (exit {result.returncode}): {result.stderr}"
        )

    # Output format:  SEQUENCE  STRUCTURE  (MFE)
    # e.g.  GGGAAACCC  (((...)))  (-1.20)
    parts = result.stdout.strip().split()
    if len(parts) < 3:
        raise RuntimeError(f"Unexpected LinearFold output: {result.stdout!r}")

    # parts[0] is the input sequence, parts[1] is the structure,
    # parts[2] is the free energy like "(-1.20)"
    structure = parts[1]
    mfe_str = parts[2].strip("()")
    try:
        mfe = float(mfe_str)
    except ValueError:
        raise RuntimeError(f"Cannot parse MFE from {mfe_str!r}")

    aup = _compute_aup_from_structure(structure) if need_aup else 0.0
    return {"mfe": mfe, "structure": structure, "aup": aup}


def estimate_fold(sequence: str, engine: str = "auto", need_aup: bool = True) -> dict:
    """
    Estimate the secondary structure of an RNA sequence.

    :param sequence: RNA sequence to fold.
    :param engine: 'auto', 'vienna', or 'linearfold'.
                   Auto prefers linearfold if available, else ViennaRNA.
    :param need_aup: Whether to compute AUP (skip if not needed to save time).
    :return: dict with keys 'mfe', 'structure', 'aup'.
    """
    if engine == "linearfold":
        return _fold_linearfold(sequence, need_aup=need_aup)

    if engine == "vienna":
        return _fold_vienna(sequence, need_aup=need_aup)

    # auto-detect: prefer LinearFold if available
    if _find_linearfold_executable():
        return _fold_linearfold(sequence, need_aup=need_aup)
    return _fold_vienna(sequence, need_aup=need_aup)
