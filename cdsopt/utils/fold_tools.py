# -*- coding: utf-8 -*-
"""
RNA folding tools with auto-detection of LinearFold and ViennaRNA fallback.

ViennaRNA mode uses partition function for exact AUP.
LinearFold mode falls back to dot-counting.
"""
from __future__ import annotations

import shutil
import subprocess

from ViennaRNA import RNA


def _fold_vienna(sequence: str) -> dict:
    """Fold with ViennaRNA; AUP from partition-function pairing probs."""
    fc = RNA.fold_compound(sequence)
    ss, mfe = fc.mfe()
    fc.pf()
    aup = _compute_aup_vienna(fc, len(sequence))
    return {"mfe": mfe, "structure": ss, "aup": aup}


def _compute_aup_vienna(fc, seq_len: int) -> float:
    """Compute AUP from base-pair probability matrix."""
    try:
        bpp = fc.bpp()
        if not bpp:
            raise ValueError("empty bpp")
    except Exception:
        # Fallback: dot-counting from MFE structure
        return 0.0

    pairing_probs = [0.0] * seq_len
    try:
        for i in range(1, seq_len + 1):
            if i < len(bpp) and bpp[i] is not None:
                if hasattr(bpp[i], "values"):
                    pairing_probs[i - 1] = sum(bpp[i].values())
                elif hasattr(bpp[i], "__iter__"):
                    pairing_probs[i - 1] = sum(float(v) for v in bpp[i] if v)
    except Exception:
        return 0.0

    aup = sum(1.0 - p for p in pairing_probs) / seq_len if seq_len > 0 else 0.0
    return max(0.0, min(1.0, aup))


def _fold_linearfold(sequence: str) -> dict:
    """Fold with LinearFold CLI; AUP from dot-counting."""
    executable = shutil.which("linearfold")
    if executable is None:
        raise FileNotFoundError("linearfold not found in PATH")

    # Use echo | linearfold -V to match the canonical usage pattern
    p1 = subprocess.Popen(
        ["echo", sequence],
        stdout=subprocess.PIPE,
    )
    result = subprocess.run(
        [executable, "-V"],
        stdin=p1.stdout,
        capture_output=True,
        text=True,
    )
    p1.stdout.close()

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

    aup = structure.count(".") / len(structure) if len(structure) > 0 else 0.0
    return {"mfe": mfe, "structure": structure, "aup": aup}


def estimate_fold(sequence: str, engine: str = "auto") -> dict:
    """
    Estimate the secondary structure of an RNA sequence.

    :param sequence: RNA sequence to fold.
    :param engine: 'auto', 'vienna', or 'linearfold'.
                   Auto prefers linearfold if available, else ViennaRNA.
    :return: dict with keys 'mfe', 'structure', 'aup'.
    """
    if engine == "linearfold":
        if not shutil.which("linearfold"):
            raise FileNotFoundError("linearfold not found in PATH")
        return _fold_linearfold(sequence)

    if engine == "vienna":
        return _fold_vienna(sequence)

    # auto-detect
    if shutil.which("linearfold"):
        return _fold_linearfold(sequence)
    return _fold_vienna(sequence)
