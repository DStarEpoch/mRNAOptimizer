# -*- coding: utf-8 -*-
"""
RNA folding tools with auto-detection of LinearFold and ViennaRNA fallback.
"""
from __future__ import annotations

import shutil
import subprocess

from ViennaRNA import RNA


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
    executable = shutil.which("linearfold")
    if executable is None:
        raise FileNotFoundError("linearfold not found in PATH")

    # Use input= instead of echo pipe for Windows compatibility.
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
        if not shutil.which("linearfold"):
            raise FileNotFoundError("linearfold not found in PATH")
        return _fold_linearfold(sequence, need_aup=need_aup)

    if engine == "vienna":
        return _fold_vienna(sequence, need_aup=need_aup)

    # auto-detect
    if shutil.which("linearfold"):
        return _fold_linearfold(sequence, need_aup=need_aup)
    return _fold_vienna(sequence, need_aup=need_aup)
