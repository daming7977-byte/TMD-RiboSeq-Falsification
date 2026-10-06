"""Overlap analysis across basal and differential TMD signal sets."""

from __future__ import annotations

from itertools import combinations
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

from .errors import ValidationFailed


def _pairwise_overlap(
    membership: pd.DataFrame,
    identifier: str,
    set_columns: list[str],
    level: str,
    eligibility_columns: Mapping[str, str] | None = None,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    eligibility_columns = eligibility_columns or {}
    for first, second in combinations(set_columns, 2):
        required_eligibility = sorted(
            {
                eligibility_columns[set_name]
                for set_name in (first, second)
                if set_name in eligibility_columns
            }
        )
        pair_membership = membership
        for eligibility_column in required_eligibility:
            pair_membership = pair_membership[
                pair_membership[eligibility_column].astype(bool)
            ]
        universe_size = len(pair_membership)
        in_first = pair_membership[first].astype(bool)
        in_second = pair_membership[second].astype(bool)
        first_n = int(in_first.sum())
        second_n = int(in_second.sum())
        overlap_n = int((in_first & in_second).sum())
        union_n = int((in_first | in_second).sum())
        neither_n = int((~in_first & ~in_second).sum())
        first_only_n = int((in_first & ~in_second).sum())
        second_only_n = int((~in_first & in_second).sum())
        odds_ratio, p_value = fisher_exact(
            [[overlap_n, first_only_n], [second_only_n, neither_n]]
        )
        rows.append(
            {
                "level": level,
                "universe_size": universe_size,
                "universe_filter": (
                    ",".join(required_eligibility)
                    if required_eligibility
                    else "all_primary_tested"
                ),
                "set_a": first,
                "set_b": second,
                "set_a_n": first_n,
                "set_b_n": second_n,
                "overlap_n": overlap_n,
                "overlap_fraction_of_a": overlap_n / first_n if first_n else np.nan,
                "overlap_fraction_of_b": overlap_n / second_n if second_n else np.nan,
                "jaccard_index": overlap_n / union_n if union_n else np.nan,
                "fisher_odds_ratio": float(odds_ratio),
                "fisher_p": float(p_value),
            }
        )
    result = pd.DataFrame(rows)
    if not result.empty:
        result["fisher_fdr"] = multipletests(
            result["fisher_p"], method="fdr_bh"
        )[1]
    return result


def run_overlap_stage(config: dict[str, Any]) -> dict[str, Any]:
    """Compare basal-high, condition-responsive, and interaction TMD sets."""

    output = Path(config["output"]["directory"])
    differential_path = output / "differential" / "primary_offset_model.tsv"
    baseline_path = output / "downstream" / "topology" / "baseline_tmd_scores.tsv"
    if not differential_path.is_file():
        raise ValidationFailed("Run the primary differential stage before overlap analysis.")
    if not baseline_path.is_file():
        raise ValidationFailed("Run the topology stage before overlap analysis.")

    threshold = float(config["statistics"]["fdr_threshold"])
    differential = pd.read_csv(differential_path, sep="\t")
    baseline = pd.read_csv(baseline_path, sep="\t")
    required = {
        "tmd_id",
        "gene_id",
        "reference_b_fdr",
        "alternate_b_fdr",
        "interaction_fdr",
    }
    missing = required - set(differential.columns)
    if missing:
        raise ValidationFailed(
            f"Primary differential table lacks overlap columns: {sorted(missing)}"
        )
    if "high_signal_tmd" not in baseline.columns:
        raise ValidationFailed("Baseline topology table lacks high_signal_tmd.")

    membership = differential[
        [
            "tmd_id",
            "gene_id",
            "reference_b_fdr",
            "alternate_b_fdr",
            "interaction_fdr",
        ]
    ].copy()
    basal = baseline[["tmd_id", "high_signal_tmd"]].copy()
    basal["high_signal_tmd"] = basal["high_signal_tmd"].astype(str).str.lower().isin(
        {"true", "1", "yes"}
    )
    membership = membership.merge(basal, on="tmd_id", how="left", validate="one_to_one")
    membership["basal_eligible"] = membership["high_signal_tmd"].notna()
    membership["basal_high_signal"] = membership.pop("high_signal_tmd").eq(True)
    membership["reference_condition"] = (
        pd.to_numeric(membership.pop("reference_b_fdr"), errors="coerce") < threshold
    )
    membership["alternate_condition"] = (
        pd.to_numeric(membership.pop("alternate_b_fdr"), errors="coerce") < threshold
    )
    membership["genotype_condition_interaction"] = (
        pd.to_numeric(membership.pop("interaction_fdr"), errors="coerce") < threshold
    )
    set_columns = [
        "basal_high_signal",
        "reference_condition",
        "alternate_condition",
        "genotype_condition_interaction",
    ]

    membership_columns = ["basal_eligible", *set_columns]
    gene_membership = membership.groupby("gene_id", as_index=False)[
        membership_columns
    ].any()
    eligibility_columns = {"basal_high_signal": "basal_eligible"}
    tmd_overlap = _pairwise_overlap(
        membership,
        "tmd_id",
        set_columns,
        "TMD",
        eligibility_columns,
    )
    gene_overlap = _pairwise_overlap(
        gene_membership,
        "gene_id",
        set_columns,
        "gene",
        eligibility_columns,
    )
    pairwise = pd.concat([tmd_overlap, gene_overlap], ignore_index=True)

    result_dir = output / "downstream" / "overlap"
    result_dir.mkdir(parents=True, exist_ok=True)
    membership_path = result_dir / "tmd_set_membership.tsv"
    gene_membership_path = result_dir / "gene_set_membership.tsv"
    pairwise_path = result_dir / "pairwise_overlap.tsv"
    membership.to_csv(membership_path, sep="\t", index=False)
    gene_membership.to_csv(gene_membership_path, sep="\t", index=False)
    pairwise.to_csv(pairwise_path, sep="\t", index=False)
    return {
        "primary_tested_tmd_universe": int(len(membership)),
        "primary_tested_gene_universe": int(len(gene_membership)),
        "basal_eligible_tmd_universe": int(membership["basal_eligible"].sum()),
        "basal_eligible_gene_universe": int(gene_membership["basal_eligible"].sum()),
        "fdr_threshold": threshold,
        "tmd_membership": str(membership_path),
        "gene_membership": str(gene_membership_path),
        "pairwise": str(pairwise_path),
    }
