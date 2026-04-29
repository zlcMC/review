"""Shared helpers for projectmd/v2 refactored scripts.

Goals
-----
* Single source of truth for small statistical / IO utilities that were copy
  pasted across ``fig04_*`` scripts (``bh_adjust``, ``parse_tsv_floats`` etc.).
* A consistent way to write outputs into ``output/v2/...`` while still
  reading dependencies from either ``output/v2/...`` (preferred) or the
  original ``output/...`` legacy location.
* Keep the existing v1 scripts untouched so the user can compare before/after.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

# Make ``workspace_paths`` importable regardless of cwd.
_REPO_ROOT = Path(__file__).resolve().parents[3]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from workspace_paths import OUTPUT_DIR, output_path, raw_data_path  # noqa: E402

V2_TAG = "v2"


def v2_dir(*parts: str) -> Path:
    """Return (and create) the directory ``output/v2/<parts>/``."""
    path = OUTPUT_DIR / V2_TAG / Path(*parts)
    path.mkdir(parents=True, exist_ok=True)
    return path


def v2_file(*parts: str) -> Path:
    """Return ``output/v2/<parts>`` and ensure its parent directory exists."""
    path = OUTPUT_DIR / V2_TAG / Path(*parts)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def dep_file(*parts: str) -> Path:
    """Locate a dependency file produced by another script.

    Prefers ``output/v2/<parts>`` (so v2 scripts chain naturally); falls back
    to the legacy ``output/<parts>`` location so each v2 script can run on its
    own without first re-running its v2 upstream sibling.
    """

    v2_candidate = OUTPUT_DIR / V2_TAG / Path(*parts)
    if v2_candidate.exists():
        return v2_candidate
    return OUTPUT_DIR / Path(*parts)


def parse_tsv_floats(values: str, dtype=np.float32) -> np.ndarray:
    """Drop-in replacement for the deprecated ``np.fromstring(s, sep='\\t')``.

    NumPy 2.0 has removed ``np.fromstring``; using ``str.split`` + ``np.array``
    is portable and roughly as fast for the small per-line arrays we use.
    """

    return np.array(values.split("\t"), dtype=dtype)


def bh_adjust(p_values: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR correction; preserves NaN entries."""

    p_values = np.asarray(p_values, dtype=np.float64)
    adjusted = np.full_like(p_values, np.nan, dtype=np.float64)
    valid = ~np.isnan(p_values)
    if valid.sum() == 0:
        return adjusted
    valid_p = p_values[valid]
    order = np.argsort(valid_p)
    ranked = valid_p[order]
    q_values = ranked * len(ranked) / (np.arange(len(ranked)) + 1)
    q_values = np.minimum.accumulate(q_values[::-1])[::-1]
    q_values = np.clip(q_values, 0, 1)
    valid_indices = np.where(valid)[0]
    adjusted[valid_indices[order]] = q_values
    return adjusted


# Cell-type / subtype normalization ------------------------------------------------

NON_TUMOR_CELL_TYPES = {"Myeloid", "Lymphocytes", "Oligodendrocytes", "Doublets"}


def subtype_group(value: object) -> str:
    """Normalize tumor_subtype string to canonical labels (PF-A/PF-B/ST-RELA/...)."""

    text = str(value).upper().replace("_", "-")
    if "PFA" in text or "PF-A" in text:
        return "PF-A"
    if "PFB" in text or "PF-B" in text:
        return "PF-B"
    if "RELA" in text or "ZFTA" in text:
        return "ST-RELA"
    if "YAP" in text:
        return "ST-YAP"
    if text in {"NAN", "", "UNKNOWN", "NONE"}:
        return "Unknown"
    return "Other"


# Key ligand/receptor genes (union of Gillen + Gojo cohort lists) -----------------

KEY_LR_GENES_CORE = [
    "MDK", "NCL", "PTN", "PTPRZ1", "LRP1", "JAG1", "NOTCH3", "NOTCH1",
    "APP", "CD74", "CXCR4", "PPIA", "BSG", "COL1A1", "COL1A2", "COL6A1",
    "SPP1", "TREM2", "TYROBP", "HBEGF", "EGFR", "CD44", "SDC2", "SDC3",
]

KEY_LR_GENES_GOJO_EXTRA = ["L1CAM", "COL9A2", "COL9A3", "MIF", "APOE"]

KEY_LR_GENES_UNION = KEY_LR_GENES_CORE + KEY_LR_GENES_GOJO_EXTRA


__all__ = [
    "OUTPUT_DIR",
    "output_path",
    "raw_data_path",
    "V2_TAG",
    "v2_dir",
    "v2_file",
    "dep_file",
    "parse_tsv_floats",
    "bh_adjust",
    "subtype_group",
    "NON_TUMOR_CELL_TYPES",
    "KEY_LR_GENES_CORE",
    "KEY_LR_GENES_GOJO_EXTRA",
    "KEY_LR_GENES_UNION",
]
