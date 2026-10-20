"""Gene-set enrichment using an explicit tested gene universe."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import fisher_exact, hypergeom
from statsmodels.stats.multitest import multipletests

from .errors import ConfigurationError, ValidationFailed


def _opportunity_stratified_exact_p(
    gene_table: pd.DataFrame, term_genes: set[str]
) -> float:
    """Exact upper-tail test conditional on positives in each opportunity stratum.

    Within every distinct tested-TMD-count stratum, the number of selected
    genes is fixed.  Term hits therefore follow a hypergeometric distribution.
    Convolution across strata gives the exact null distribution for the total
    number of selected term genes without relying on an asymptotic small-sample
    logistic approximation.
    """

    observed = int(
        gene_table.loc[gene_table["gene_id"].isin(term_genes), "selected"].sum()
    )
    total_pmf = np.array([1.0])
    for _, stratum in gene_table.groupby("tested_tmd_count"):
        population = int(len(stratum))
        selected = int(stratum["selected"].sum())
        term_size = int(stratum["gene_id"].isin(term_genes).sum())
        lower = max(0, selected - (population - term_size))
        upper = min(selected, term_size)
        support = np.arange(lower, upper + 1, dtype=int)
        stratum_pmf = np.zeros(upper + 1, dtype=float)
        stratum_pmf[support] = hypergeom.pmf(
            support, population, term_size, selected
        )
        total_pmf = np.convolve(total_pmf, stratum_pmf)
    return float(total_pmf[observed:].sum())


def _read_mapping(options: dict[str, Any]) -> pd.DataFrame:
    mapping = options.get("mapping")
    if not isinstance(mapping, dict):
        raise ConfigurationError("enrichment requires an explicit mapping block.")
    path = Path(str(mapping.get("path", "")))
    if not path.is_file():
        raise ValidationFailed(f"Enrichment mapping table not found: {path}")
    separator = "," if path.suffix.lower() == ".csv" else "\t"
    table = pd.read_csv(path, sep=separator, dtype=str)
    gene_column = str(mapping.get("gene_column", "gene_id"))
    term_column = str(mapping.get("term_column", "term_id"))
    ontology_column = mapping.get("ontology_column")
    term_name_column = mapping.get("term_name_column")
    required = [gene_column, term_column]
    for optional in [ontology_column, term_name_column]:
        if optional:
            required.append(str(optional))
    missing = [column for column in required if column not in table.columns]
    if missing:
        raise ValidationFailed(f"Enrichment mapping lacks columns: {missing}")
    renamed = table[required].rename(
        columns={
            gene_column: "gene_id",
            term_column: "term_id",
            **({str(ontology_column): "ontology"} if ontology_column else {}),
            **({str(term_name_column): "term_name"} if term_name_column else {}),
        }
    )
    renamed = renamed.dropna(subset=["gene_id", "term_id"]).drop_duplicates()
    return renamed


def run_enrichment_stage(config: dict[str, Any]) -> dict[str, Any]:
    """Test interaction genes against terms within the analyzable universe."""

    output = Path(config["output"]["directory"])
    differential_path = output / "differential" / "primary_offset_model.tsv"
    if not differential_path.is_file():
        raise ValidationFailed("Run the primary differential stage before enrichment.")
    options = config.get("optional_modules", {}).get("enrichment", {})
    annotations = _read_mapping(options)
    differential = pd.read_csv(differential_path, sep="\t")
    if not {"gene_id", "tmd_id", "interaction_fdr"}.issubset(differential.columns):
        raise ValidationFailed(
            "Primary differential output lacks gene_id, tmd_id, or interaction_fdr."
        )
    threshold = float(config["statistics"]["fdr_threshold"])
    gene_table = (
        differential.groupby("gene_id", as_index=False)
        .agg(
            tested_tmd_count=("tmd_id", "nunique"),
            minimum_interaction_fdr=("interaction_fdr", "min"),
        )
    )
    gene_table["selected"] = gene_table["minimum_interaction_fdr"] < threshold
    gene_table["log1p_tested_tmd_count"] = np.log1p(gene_table["tested_tmd_count"])
    universe = set(gene_table["gene_id"].astype(str))
    annotations = annotations[annotations["gene_id"].isin(universe)].copy()
    minimum_term_genes = int(options.get("minimum_term_genes", 3))
    if minimum_term_genes < 2:
        raise ConfigurationError("enrichment.minimum_term_genes must be at least 2.")

    selected = set(gene_table.loc[gene_table["selected"], "gene_id"])
    rows: list[dict[str, Any]] = []
    group_columns = ["term_id"]
    if "ontology" in annotations.columns:
        group_columns.insert(0, "ontology")
    for key, group in annotations.groupby(group_columns, dropna=False):
        if not isinstance(key, tuple):
            key = (key,)
        key_values = dict(zip(group_columns, key))
        term_genes = set(group["gene_id"])
        if len(term_genes) < minimum_term_genes:
            continue
        selected_in_term = len(selected & term_genes)
        selected_not_term = len(selected - term_genes)
        unselected_in_term = len((universe - selected) & term_genes)
        unselected_not_term = len((universe - selected) - term_genes)
        odds_ratio, fisher_p = fisher_exact(
            [[selected_in_term, selected_not_term],
             [unselected_in_term, unselected_not_term]]
        )

        model = gene_table.copy()
        model["in_term"] = model["gene_id"].isin(term_genes).astype(int)
        adjusted_beta = adjusted_se = adjusted_p = np.nan
        if model["selected"].nunique() == 2 and model["in_term"].nunique() == 2:
            try:
                design = sm.add_constant(
                    model[["in_term", "log1p_tested_tmd_count"]], has_constant="add"
                )
                fit = sm.GLM(
                    model["selected"].astype(int),
                    design,
                    family=sm.families.Binomial(),
                ).fit()
                adjusted_beta = float(fit.params["in_term"])
                adjusted_se = float(fit.bse["in_term"])
                adjusted_p = float(fit.pvalues["in_term"])
            except (ValueError, np.linalg.LinAlgError, RuntimeError):
                pass
        row = {
            **key_values,
            "term_name": (
                str(group["term_name"].dropna().iloc[0])
                if "term_name" in group and not group["term_name"].dropna().empty
                else ""
            ),
            "universe_size": len(universe),
            "selected_size": len(selected),
            "term_size": len(term_genes),
            "selected_in_term": selected_in_term,
            "fisher_odds_ratio": float(odds_ratio),
            "fisher_p": float(fisher_p),
            "opportunity_adjusted_log_odds": adjusted_beta,
            "opportunity_adjusted_se": adjusted_se,
            "opportunity_adjusted_p": adjusted_p,
            "opportunity_stratified_exact_p": _opportunity_stratified_exact_p(
                gene_table, term_genes
            ),
        }
        rows.append(row)
    results = pd.DataFrame(rows)
    if results.empty:
        raise ValidationFailed("No enrichment term met the configured minimum size.")
    results["fisher_fdr"] = multipletests(results["fisher_p"], method="fdr_bh")[1]
    valid_adjusted = results["opportunity_adjusted_p"].notna()
    results["opportunity_adjusted_fdr"] = np.nan
    if valid_adjusted.any():
        results.loc[valid_adjusted, "opportunity_adjusted_fdr"] = multipletests(
            results.loc[valid_adjusted, "opportunity_adjusted_p"], method="fdr_bh"
        )[1]
    results["opportunity_stratified_exact_fdr"] = multipletests(
        results["opportunity_stratified_exact_p"], method="fdr_bh"
    )[1]
    results = results.sort_values(["fisher_fdr", "fisher_p"])

    result_dir = output / "downstream" / "enrichment"
    result_dir.mkdir(parents=True, exist_ok=True)
    universe_path = result_dir / "tested_gene_universe.tsv"
    results_path = result_dir / "gene_set_enrichment.tsv"
    gene_table.to_csv(universe_path, sep="\t", index=False)
    results.to_csv(results_path, sep="\t", index=False)
    return {
        "tested_gene_universe": int(len(universe)),
        "selected_interaction_genes": int(len(selected)),
        "tested_terms": int(len(results)),
        "universe": str(universe_path),
        "results": str(results_path),
    }
