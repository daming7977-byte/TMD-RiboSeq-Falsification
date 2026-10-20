"""Protein-aware topology summaries for baseline TMD signals."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import fisher_exact, mannwhitneyu, norm

from .errors import ConfigurationError, ValidationFailed


def run_topology_stage(config: dict[str, Any]) -> dict[str, Any]:
    """Summarize primary-window signals with equal protein weighting."""

    output = Path(config["output"]["directory"])
    profiles_path = output / "tmd_profiles.tsv.gz"
    if not profiles_path.is_file():
        raise ValidationFailed(f"Profile table not found: {profiles_path}")
    options = config.get("optional_modules", {}).get("topology", {})
    role = str(options.get("role", config["samples"].get("baseline_role", "baseline")))
    baseline_filters = {
        str(key): str(value)
        for key, value in config["samples"].get("baseline_filters", {}).items()
    }
    protein_metric = str(options.get("protein_metric", "mean"))
    if protein_metric not in {"mean", "median", "maximum"}:
        raise ConfigurationError(
            "optional_modules.topology.protein_metric must be mean, median, or maximum."
        )
    high_quantile = float(options.get("high_signal_quantile", 0.90))
    if not 0 < high_quantile < 1:
        raise ConfigurationError("high_signal_quantile must be in (0, 1).")

    profiles = pd.read_csv(profiles_path, sep="\t", low_memory=False)
    center = int(config["profiles"]["primary_center_nt"])
    baseline = profiles[
        (profiles["center_nt"] == center) & (profiles["role"].astype(str) == role)
    ].copy()
    for factor, value in baseline_filters.items():
        if factor not in baseline.columns:
            raise ConfigurationError(
                f"samples.baseline_filters refers to absent profile column {factor!r}."
            )
        baseline = baseline[baseline[factor].astype(str) == value].copy()
    if baseline.empty:
        raise ValidationFailed(
            f"No profile rows have role={role!r} at primary center {center}."
        )
    sample_ids = sorted(baseline["sample_id"].astype(str).unique())
    pass_counts = baseline.groupby("tmd_id")["coverage_pass"].agg(
        rows="size", passed="sum"
    )
    eligible_ids = pass_counts.index[
        (pass_counts["rows"] == len(sample_ids))
        & (pass_counts["passed"] == len(sample_ids))
    ]
    baseline = baseline[baseline["tmd_id"].isin(eligible_ids)].copy()
    if baseline.empty:
        raise ValidationFailed("No TMD passes coverage in every selected baseline sample.")

    tmd_scores = (
        baseline.groupby(["tmd_id", "gene_id", "transcript_id"], as_index=False)
        .agg(
            mean_signal=("local_log2_occupancy", "mean"),
            mean_background_count=("background_count", "mean"),
            n_samples=("sample_id", "nunique"),
        )
    )
    all_tmds = profiles[["gene_id", "tmd_id"]].drop_duplicates()
    tmd_counts = all_tmds.groupby("gene_id")["tmd_id"].nunique()
    tmd_scores["annotated_tmd_count"] = tmd_scores["gene_id"].map(tmd_counts)
    tmd_scores["topology"] = np.where(
        tmd_scores["annotated_tmd_count"] == 1, "single-pass", "multi-pass"
    )
    tmd_high_cutoff = float(tmd_scores["mean_signal"].quantile(high_quantile))
    tmd_scores["high_signal_tmd"] = tmd_scores["mean_signal"] >= tmd_high_cutoff

    protein = (
        tmd_scores.groupby(["gene_id", "topology", "annotated_tmd_count"], as_index=False)
        .agg(
            protein_mean_signal=("mean_signal", "mean"),
            protein_median_signal=("mean_signal", "median"),
            protein_maximum_signal=("mean_signal", "max"),
            fraction_high_signal_tmds=("high_signal_tmd", "mean"),
            observed_tmd_count=("tmd_id", "nunique"),
        )
    )
    metric_column = f"protein_{protein_metric}_signal"
    protein_cutoff = float(protein[metric_column].quantile(high_quantile))
    protein["high_signal_protein"] = protein[metric_column] >= protein_cutoff

    single = protein.loc[protein["topology"] == "single-pass", metric_column]
    multi = protein.loc[protein["topology"] == "multi-pass", metric_column]
    mannwhitney_p = np.nan
    if len(single) > 0 and len(multi) > 0:
        mannwhitney_p = float(mannwhitneyu(single, multi, alternative="two-sided").pvalue)
    table = pd.crosstab(protein["topology"], protein["high_signal_protein"]).reindex(
        index=["single-pass", "multi-pass"], columns=[True, False], fill_value=0
    )
    odds_ratio, fisher_p = fisher_exact(
        [
            [table.loc["single-pass", True], table.loc["single-pass", False]],
            [table.loc["multi-pass", True], table.loc["multi-pass", False]],
        ]
    )

    sensitivity_quantiles = options.get(
        "tail_sensitivity_quantiles", [0.95, 0.90, 0.85, 0.80]
    )
    tail_rows: list[dict[str, float | int]] = []
    for quantile in [float(value) for value in sensitivity_quantiles]:
        if not 0 < quantile < 1:
            raise ConfigurationError(
                "tail_sensitivity_quantiles values must be in (0, 1)."
            )
        cutoff = float(protein[metric_column].quantile(quantile))
        high = protein[metric_column] >= cutoff
        sensitivity_table = pd.crosstab(protein["topology"], high).reindex(
            index=["single-pass", "multi-pass"],
            columns=[True, False],
            fill_value=0,
        )
        cells = np.array(
            [
                sensitivity_table.loc["single-pass", True],
                sensitivity_table.loc["single-pass", False],
                sensitivity_table.loc["multi-pass", True],
                sensitivity_table.loc["multi-pass", False],
            ],
            dtype=float,
        )
        sensitivity_or, sensitivity_p = fisher_exact(
            [[cells[0], cells[1]], [cells[2], cells[3]]]
        )
        ci_cells = cells + 0.5 if (cells == 0).any() else cells
        log_or = np.log((ci_cells[0] * ci_cells[3]) / (ci_cells[1] * ci_cells[2]))
        log_or_se = np.sqrt((1.0 / ci_cells).sum())
        critical = float(norm.ppf(0.975))
        tail_rows.append(
            {
                "upper_tail_percent": round(100.0 * (1.0 - quantile), 6),
                "quantile": quantile,
                "signal_cutoff": cutoff,
                "single_high_n": int(cells[0]),
                "single_total_n": int(cells[0] + cells[1]),
                "single_high_fraction": float(cells[0] / (cells[0] + cells[1])),
                "multi_high_n": int(cells[2]),
                "multi_total_n": int(cells[2] + cells[3]),
                "multi_high_fraction": float(cells[2] / (cells[2] + cells[3])),
                "odds_ratio": float(sensitivity_or),
                "odds_ratio_ci95_low": float(np.exp(log_or - critical * log_or_se)),
                "odds_ratio_ci95_high": float(np.exp(log_or + critical * log_or_se)),
                "fisher_p": float(sensitivity_p),
            }
        )

    tmd_scores["single_pass"] = (tmd_scores["topology"] == "single-pass").astype(int)
    tmd_scores["log2_background"] = np.log2(tmd_scores["mean_background_count"] + 1.0)
    clustered_beta = clustered_se = clustered_p = np.nan
    if tmd_scores["single_pass"].nunique() == 2 and tmd_scores["gene_id"].nunique() >= 3:
        design = sm.add_constant(tmd_scores[["single_pass", "log2_background"]])
        fit = sm.OLS(tmd_scores["mean_signal"], design).fit(
            cov_type="cluster", cov_kwds={"groups": tmd_scores["gene_id"]}
        )
        clustered_beta = float(fit.params["single_pass"])
        clustered_se = float(fit.bse["single_pass"])
        clustered_p = float(fit.pvalues["single_pass"])

    summary = pd.DataFrame(
        [
            {
                "center_nt": center,
                "baseline_role": role,
                "baseline_samples": len(sample_ids),
                "strict_baseline_tmds": len(tmd_scores),
                "proteins": len(protein),
                "single_pass_proteins": int((protein["topology"] == "single-pass").sum()),
                "multi_pass_proteins": int((protein["topology"] == "multi-pass").sum()),
                "protein_metric": protein_metric,
                "high_signal_quantile": high_quantile,
                "protein_high_signal_cutoff": protein_cutoff,
                "single_high_n": int(table.loc["single-pass", True]),
                "single_total_n": int(table.loc["single-pass"].sum()),
                "single_high_fraction": float(
                    table.loc["single-pass", True] / max(1, table.loc["single-pass"].sum())
                ),
                "multi_high_n": int(table.loc["multi-pass", True]),
                "multi_total_n": int(table.loc["multi-pass"].sum()),
                "multi_high_fraction": float(
                    table.loc["multi-pass", True] / max(1, table.loc["multi-pass"].sum())
                ),
                "high_tail_odds_ratio": float(odds_ratio),
                "high_tail_fisher_p": float(fisher_p),
                "continuous_mannwhitney_p": mannwhitney_p,
                "clustered_single_pass_beta": clustered_beta,
                "clustered_single_pass_se": clustered_se,
                "clustered_single_pass_p": clustered_p,
            }
        ]
    )
    result_dir = output / "downstream" / "topology"
    result_dir.mkdir(parents=True, exist_ok=True)
    tmd_path = result_dir / "baseline_tmd_scores.tsv"
    protein_path = result_dir / "protein_summary.tsv"
    summary_path = result_dir / "topology_summary.tsv"
    tail_path = result_dir / "tail_cutoff_sensitivity.tsv"
    tmd_scores.to_csv(tmd_path, sep="\t", index=False)
    protein.to_csv(protein_path, sep="\t", index=False)
    summary.to_csv(summary_path, sep="\t", index=False)
    pd.DataFrame(tail_rows).to_csv(tail_path, sep="\t", index=False)
    return {
        "baseline_tmds": int(len(tmd_scores)),
        "proteins": int(len(protein)),
        "tmd_scores": str(tmd_path),
        "protein_summary": str(protein_path),
        "summary": str(summary_path),
        "tail_cutoff_sensitivity": str(tail_path),
    }
