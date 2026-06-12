# -*- coding: utf-8 -*-
"""
Make GEMORNA importable from utropt without subprocess.

GEMORNA lives in `submodules/GEMORNA`.  This adapter temporarily adds its
`src/` directory to `sys.path`, imports the prediction models and helper
utilities, and exposes a small functional API:

    from utropt.gemorna_adapter import predict_5utr, predict_3utr
    score = predict_5utr("GCCGCCACCATGGGAGAATAAACTAGT")

The first call loads the model checkpoint, which is cached for reuse.
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

# Path to the GEMORNA submodule source tree.
_GEMORNA_SRC = Path(__file__).resolve().parent.parent / "submodules" / "GEMORNA" / "src"
_GEMORNA_CKPT = Path(__file__).resolve().parent.parent / "submodules" / "GEMORNA" / "checkpoints"


def _ensure_importable() -> None:
    """Add GEMORNA src/ to sys.path if not already present."""
    src_str = str(_GEMORNA_SRC)
    if src_str not in sys.path:
        sys.path.insert(0, src_str)


# Lazily populated model singletons.
_5UTR_MODEL: Any | None = None
_3UTR_MODEL: Any | None = None


def _load_5utr_model() -> Any:
    global _5UTR_MODEL
    if _5UTR_MODEL is not None:
        return _5UTR_MODEL

    _ensure_importable()
    import torch
    import models.model_pred5UTR as model
    from shared.helper import kernel_sizes_5UTR, scale, tokenize, vocab

    class Args:
        embed_num = 10
        embed_dim = 64
        kernel_num = 128
        kernel_sizes = kernel_sizes_5UTR
        dropout = 0.1

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    predictor = model.Model(Args()).to(device)
    ckpt = _GEMORNA_CKPT / "5utr.pt"
    predictor.load_state_dict(torch.load(ckpt, map_location=device), strict=True)
    predictor.eval()

    _5UTR_MODEL = (predictor, device, tokenize, scale, vocab)
    logger.info("Loaded GEMORNA 5'UTR predictor from %s", ckpt)
    return _5UTR_MODEL


def _load_3utr_model() -> Any:
    global _3UTR_MODEL
    if _3UTR_MODEL is not None:
        return _3UTR_MODEL

    _ensure_importable()
    import torch
    import models.model_pred3UTR as model
    from shared.helper import kernel_sizes_3UTR, tokenize

    class Args:
        embed_num = 10
        embed_dim = 256
        kernel_num = 200
        kernel_sizes = kernel_sizes_3UTR
        dropout = 0.1

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    predictor = model.Model(Args()).to(device)
    ckpt = _GEMORNA_CKPT / "3utr.pt"
    predictor.load_state_dict(torch.load(ckpt, map_location=device), strict=True)
    predictor.eval()

    _3UTR_MODEL = (predictor, device, tokenize)
    logger.info("Loaded GEMORNA 3'UTR predictor from %s", ckpt)
    return _3UTR_MODEL


def predict_5utr(sequence: str) -> float:
    """Predict 5'UTR quality (TIE proxy) using GEMORNA.

    :param sequence: RNA sequence (A/C/G/U/T, length <= 100 for the model).
    :return: GEMORNA prediction score (scaled MRL-like value).
    """
    import torch

    predictor, device, tokenize, scale, vocab = _load_5utr_model()
    seq = sequence.upper().replace("T", "U")
    tokens = tokenize(seq)
    padded = tokens + [vocab["[PAD]"]] * (100 - len(tokens))
    with torch.no_grad():
        pred = predictor(torch.tensor([padded], device=device)).squeeze().cpu().numpy()
    return float(scale(pred))


def predict_3utr(sequence: str) -> float:
    """Predict 3'UTR quality (stability proxy) using GEMORNA.

    :param sequence: RNA sequence (A/C/G/U/T).
    :return: GEMORNA prediction score.
    """
    import torch

    predictor, device, tokenize = _load_3utr_model()
    seq = sequence.upper().replace("T", "U")
    tokens = tokenize(seq)
    with torch.no_grad():
        pred = predictor(torch.tensor([tokens], device=device)).squeeze().cpu().numpy()
    return float(pred)
