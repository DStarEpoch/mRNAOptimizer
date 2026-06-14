# -*- coding: utf-8 -*-
"""
Thin wrapper around GEMORNA's ``src/generate.py`` for UTR generation.

Simply calls the official GEMORNA CLI ``n`` times and returns the generated
sequences.  Generation only works on Linux/macOS/WSL because GEMORNA ships
Linux-only ``.so`` files.
"""
from __future__ import annotations

import logging
import platform
import subprocess
from pathlib import Path
from typing import Literal

logger = logging.getLogger(__name__)

_GEMORNA_ROOT = Path(__file__).resolve().parent.parent / "submodules" / "GEMORNA"

UTRLength = Literal["short", "medium", "long"]


def _require_non_windows() -> None:
    if platform.system() == "Windows":
        raise RuntimeError(
            "GEMORNA UTR generation is not supported on Windows. "
            "Please run on Linux, macOS, or WSL2."
        )


def _generate(mode: str, length: UTRLength, n: int) -> list[str]:
    _require_non_windows()

    ckpt = _GEMORNA_ROOT / "checkpoints" / f"gemorna_{mode}.pt"
    script = _GEMORNA_ROOT / "src" / "generate.py"

    sequences: list[str] = []
    for i in range(n):
        result = subprocess.run(
            [
                "python",
                str(script),
                "--mode",
                mode,
                "--ckpt_path",
                str(ckpt),
                "--utr_length",
                length,
            ],
            cwd=str(_GEMORNA_ROOT),
            capture_output=True,
            text=True,
            check=False,
        )
        if result.returncode != 0:
            logger.error("GEMORNA %s generation failed: %s", mode, result.stderr)
            raise RuntimeError(f"GEMORNA {mode} generation failed: {result.stderr}")

        lines = [line.strip() for line in result.stdout.splitlines() if line.strip()]
        if not lines:
            raise RuntimeError(f"GEMORNA {mode} produced no output")
        sequences.append(lines[-1])
        logger.info("Generated %s %d/%d", mode, i + 1, n)

    return sequences


def generate_5utr_candidates(length: UTRLength = "short", n: int = 5) -> list[str]:
    """Generate ``n`` 5'UTR candidates using GEMORNA-5UTR."""
    return _generate("5utr", length, n)


def generate_3utr_candidates(length: UTRLength = "long", n: int = 5) -> list[str]:
    """Generate ``n`` 3'UTR candidates using GEMORNA-3UTR."""
    return _generate("3utr", length, n)
