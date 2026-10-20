"""Independent target/background statistical sensitivity analysis."""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Any

import pandas as pd

from .differential import _edgeR_version, _r_script, _write_design
from .errors import ConfigurationError, ValidationFailed


def _safe_spearman(left: pd.Series, right: pd.Series) -> float | None:
    value = left.corr(right, method="spearman")
    return None if pd.isna(value) else float(value)


def _safe_pearson(left: pd.Series, right: pd.Series) -> float | None:
    value = left.corr(right, method="pearson")
    return None if pd.isna(value) else float(value)


def run_target_background_sensitivity(config: dict[str, Any]) -> dict[str, Any]:
    """Fit differential feature usage and compare it with the offset model."""

    _edgeR_version()
    output = Path(config["output"]["directory"])
    profiles = output / "tmd_profiles.tsv.gz"
    primary_path = output / "differential" / "primary_offset_model.tsv"
    if not profiles.is_file() or not primary_path.is_file():
        raise ValidationFailed(
            "The profile and primary differential stages must complete before sensitivity analysis."
        )
    differential_dir = output / "differential"
    design_path = differential_dir / "resolved_design.tsv"
    factors = _write_design(design_path, config)
    references = config["design"].get("reference_levels", {})
    if any(name not in references for name in factors):
        raise ConfigurationError("design.reference_levels must define every factor.")
    result_path = differential_dir / "target_background_sensitivity.tsv"
    script = _r_script("fit_target_background_sensitivity.R")
    command = [
        "Rscript",
        str(script),
        str(profiles),
        str(design_path),
        str(result_path),
        factors[0],
        factors[1],
        str(references[factors[0]]),
        str(references[factors[1]]),
        str(int(config["profiles"]["primary_center_nt"])),
    ]
    completed = subprocess.run(command, check=False, capture_output=True, text=True)
    if completed.returncode != 0:
        raise ValidationFailed(
            "Target/background sensitivity model failed.\n"
            f"STDOUT:\n{completed.stdout}\nSTDERR:\n{completed.stderr}"
        )

    primary = pd.read_csv(primary_path, sep="\t")
    sensitivity = pd.read_csv(result_path, sep="\t")
    merged = primary.merge(sensitivity, on=["tmd_id", "gene_id", "transcript_id"])
    primary_effect = merged["interaction_log2_effect"]
    sensitivity_effect = merged[
        "interaction_target_minus_background_log2_effect"
    ]
    fdr_threshold = float(config["statistics"]["fdr_threshold"])
    primary_sig = merged["interaction_fdr_x"] < fdr_threshold
    sensitivity_sig = merged["interaction_fdr_y"] < fdr_threshold
    primary_tail = merged.loc[primary_sig].copy()
    tail_primary_effect = primary_tail["interaction_log2_effect"]
    tail_sensitivity_effect = primary_tail[
        "interaction_target_minus_background_log2_effect"
    ]
    tail_same_direction = (
        (tail_primary_effect * tail_sensitivity_effect) > 0
    )
    comparison = {
        "fdr_threshold": fdr_threshold,
        "tested_in_both_models": int(len(merged)),
        "effect_spearman": _safe_spearman(primary_effect, sensitivity_effect),
        "effect_direction_agreement_fraction": float(
            ((primary_effect * sensitivity_effect) > 0).mean()
        ),
        "primary_significant": int(primary_sig.sum()),
        "sensitivity_significant": int(sensitivity_sig.sum()),
        "significant_in_both": int((primary_sig & sensitivity_sig).sum()),
        "primary_significant_same_direction": int(
            (primary_sig & ((primary_effect * sensitivity_effect) > 0)).sum()
        ),
        "primary_significant_subset_n": int(len(primary_tail)),
        "primary_significant_subset_effect_spearman": _safe_spearman(
            tail_primary_effect, tail_sensitivity_effect
        ),
        "primary_significant_subset_effect_pearson": _safe_pearson(
            tail_primary_effect, tail_sensitivity_effect
        ),
        "primary_significant_subset_same_direction": int(
            tail_same_direction.sum()
        ),
        "primary_significant_subset_significant_in_both": int(
            (primary_tail["interaction_fdr_y"] < fdr_threshold).sum()
        ),
    }
    comparison_path = differential_dir / "model_sensitivity_comparison.tsv"
    pd.DataFrame([comparison]).to_csv(comparison_path, sep="\t", index=False)
    merged.to_csv(
        differential_dir / "primary_and_target_background_results.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    tail_path = differential_dir / "primary_significant_model_comparison.tsv"
    primary_tail.to_csv(tail_path, sep="\t", index=False)
    return {
        "results": str(result_path),
        "comparison": str(comparison_path),
        "primary_significant_comparison": str(tail_path),
        **comparison,
    }
