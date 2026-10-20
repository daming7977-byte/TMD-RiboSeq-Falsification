#!/usr/bin/env python3
"""Create compact diagnostics from frozen first-reveal machine outputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", required=True, type=Path)
    args = parser.parse_args()
    result = args.project_root / "phase_b_first_reveal"
    scans = pd.read_csv(result / "spatial_scans.tsv", sep="\t")
    bands_path = result / "spatial_bands.tsv"
    bands = pd.read_csv(bands_path, sep="\t") if bands_path.stat().st_size else pd.DataFrame()
    primary_mean = scans.loc[scans.analysis == "primary", "effect"].mean()
    summaries = []
    for name, rows in scans.groupby("analysis", sort=False):
        peak = rows.loc[rows.effect.abs().idxmax()]
        selected_bands = bands.loc[bands.analysis == name] if len(bands) else pd.DataFrame()
        summaries.append(
            {
                "analysis": name,
                "domain_mean_effect": rows.effect.mean(),
                "minimum_effect": rows.effect.min(),
                "maximum_effect": rows.effect.max(),
                "absolute_peak_center_nt": int(peak.center_nt),
                "absolute_peak_effect": peak.effect,
                "primary_mean_attenuation_fraction": 1.0 - abs(rows.effect.mean()) / abs(primary_mean),
                "minimum_maxT_FWER_p": rows.maxT_FWER_p.min(),
                "simultaneously_significant_centers": int(((rows.simultaneous_95_low > 0) | (rows.simultaneous_95_high < 0)).sum()),
                "band_count": len(selected_bands),
                "bands": ";".join(
                    f"{int(row.start_nt)}:{int(row.end_nt)}:{int(row.direction)}"
                    for row in selected_bands.itertuples()
                ),
                "minimum_TMDs": int(rows.n_tmds.min()),
                "maximum_TMDs": int(rows.n_tmds.max()),
                "minimum_genes": int(rows.n_genes.min()),
                "maximum_genes": int(rows.n_genes.max()),
            }
        )
    pd.DataFrame(summaries).to_csv(
        result / "sensitivity_summary.tsv", sep="\t", index=False
    )

    coverage = pd.read_csv(result / "primary_coverage_by_position.tsv", sep="\t")
    fixed = pd.read_csv(result / "fixed_all_positions_population_curve.tsv", sep="\t")
    loo = pd.read_csv(result / "leave_one_gene_out.tsv", sep="\t")
    membership = pd.read_csv(result / "primary_coverage_membership.tsv.gz", sep="\t")
    eligible = membership.loc[membership.eligible & membership.effect.notna()].copy()
    eligible["local_pooled"] = eligible.local_0 + eligible.local_30
    per_tmd = eligible.groupby("tmd_id", as_index=False).agg(
        band_effect=("effect", "mean"),
        mean_local_pooled=("local_pooled", "mean"),
        background_0=("background_0", "first"),
        background_30=("background_30", "first"),
    )
    influence = {
        "adjacent_membership_jaccard_min": float(coverage.adjacent_membership_jaccard.min()),
        "adjacent_membership_jaccard_max": float(coverage.adjacent_membership_jaccard.max()),
        "primary_TMD_count_min": int(coverage.n_tmds.min()),
        "primary_TMD_count_max": int(coverage.n_tmds.max()),
        "primary_gene_count_min": int(coverage.n_genes.min()),
        "primary_gene_count_max": int(coverage.n_genes.max()),
        "fixed_population_TMDs": int(fixed.n_tmds.iloc[0]),
        "fixed_population_genes": int(fixed.n_genes.iloc[0]),
        "fixed_population_domain_mean_effect": float(fixed.effect.mean()),
        "primary_vs_fixed_curve_Pearson": float(np.corrcoef(scans.loc[scans.analysis == "primary", "effect"], fixed.effect)[0, 1]),
        "leave_one_gene_out_max_attenuation": float(loo.proportional_attenuation.max()),
        "leave_one_gene_out_direction_reversals": int(loo.direction_reversal.sum()),
        "effect_vs_local_coverage_Spearman": float(spearmanr(per_tmd.band_effect, per_tmd.mean_local_pooled).statistic),
        "effect_vs_background0_Spearman": float(spearmanr(per_tmd.band_effect, per_tmd.background_0).statistic),
        "effect_vs_background30_Spearman": float(spearmanr(per_tmd.band_effect, per_tmd.background_30).statistic),
    }
    (result / "coverage_influence_summary.json").write_text(
        json.dumps(influence, indent=2, sort_keys=True)
    )

    pseudo = pd.read_csv(result / "pseudo_anchor_scan.tsv", sep="\t")
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    first = axes[0]
    for name, label, color in [
        ("primary", "Primary non-deduplicated", "#b2182b"),
        ("umi_primary", "UMI sensitivity", "#2166ac"),
        ("Sec63_noCHX_30_minus_0", "Sec63 no-CHX", "#4d9221"),
    ]:
        rows = scans.loc[scans.analysis == name]
        first.plot(rows.center_nt, rows.effect, label=label, color=color, linewidth=1.8)
    first.plot(pseudo.center_nt, pseudo.effect, label="Matched pseudo-anchors", color="#762a83", linewidth=1.8)
    first.axhline(0, color="black", linewidth=0.7)
    first.set_ylabel("Δlog2(local/background)")
    first.legend(frameon=False, ncol=2)
    first.set_title("Primary signal versus prespecified falsification controls")

    second = axes[1]
    for name, label in [
        ("sample_specific_psite", "Sample-specific P-site"),
        ("author_L_minus_15", "Author-like L−15"),
        ("adjacent_L_minus_16", "L−16"),
        ("window_15nt", "15-nt window"),
        ("window_31nt", "31-nt window"),
        ("neighboring_TMD_excluded", "Neighboring TMD excluded"),
    ]:
        rows = scans.loc[scans.analysis == name]
        second.plot(rows.center_nt, rows.effect, label=label, linewidth=1.2)
    second.axhline(0, color="black", linewidth=0.7)
    second.set_xlabel("Position relative to TMD start (nt)")
    second.set_ylabel("Δlog2(local/background)")
    second.legend(frameon=False, ncol=3, fontsize=8)
    second.set_title("Prespecified coordinate, window, and coverage sensitivities")
    fig.suptitle("Phase B first-reveal spatial diagnostics", y=0.995)
    fig.tight_layout()
    fig.savefig(result / "first_reveal_spatial_diagnostics.png", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
