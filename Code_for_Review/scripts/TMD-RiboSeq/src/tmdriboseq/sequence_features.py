"""Optional sequence and physicochemical associations with baseline signals."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.outliers_influence import variance_inflation_factor

from .errors import ConfigurationError, ValidationFailed
from .tmds import load_canonical_tmds, load_paper_legacy_tmds


HYDROPATHY = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5,
    "M": 1.9, "A": 1.8, "G": -0.4, "T": -0.7, "S": -0.8,
    "W": -0.9, "Y": -1.3, "P": -1.6, "H": -3.2, "E": -3.5,
    "Q": -3.5, "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.5,
}


def _read_fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, list[str]] = {}
    identifier: str | None = None
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                identifier = line[1:].split()[0]
                if identifier in sequences:
                    raise ValidationFailed(f"Duplicate protein FASTA ID: {identifier}")
                sequences[identifier] = []
            elif identifier is not None:
                sequences[identifier].append(line.upper())
    return {key: "".join(value) for key, value in sequences.items()}


def _fraction(sequence: str, residues: str) -> float:
    if not sequence:
        return float("nan")
    return sum(residue in residues for residue in sequence) / len(sequence)


def _net_charge(sequence: str) -> int:
    return sum(residue in "KR" for residue in sequence) - sum(
        residue in "DE" for residue in sequence
    )


def _protein_mapping(options: dict[str, Any]) -> tuple[str, dict[str, str], str]:
    mapping_options = options.get("protein_mapping")
    if not isinstance(mapping_options, dict):
        raise ConfigurationError(
            "sequence_features requires an explicit protein_mapping block."
        )
    path = Path(str(mapping_options.get("path", "")))
    key_type = str(mapping_options.get("key", "transcript_id"))
    key_column = str(mapping_options.get("key_column", key_type))
    protein_column = str(mapping_options.get("protein_id_column", "protein_id"))
    key_suffix_to_strip = str(mapping_options.get("key_suffix_to_strip", ""))
    if key_type not in {"gene_id", "transcript_id"}:
        raise ConfigurationError("protein_mapping.key must be gene_id or transcript_id.")
    if not path.is_file():
        raise ValidationFailed(f"Protein mapping table not found: {path}")
    delimiter = "," if path.suffix.lower() == ".csv" else "\t"
    mapping: dict[str, str] = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if key_column not in (reader.fieldnames or []) or protein_column not in (
            reader.fieldnames or []
        ):
            raise ValidationFailed(
                f"Protein mapping requires columns {key_column!r} and {protein_column!r}."
            )
        for row in reader:
            key = row[key_column].strip()
            value = row[protein_column].strip()
            if key in mapping and mapping[key] != value:
                raise ValidationFailed(f"Ambiguous protein mapping for {key!r}.")
            mapping[key] = value
    return key_type, mapping, key_suffix_to_strip


def _load_tmd_records(config: dict[str, Any]):
    system = config["tmd_annotation"]["coordinate_system"]
    if system == "cds_nt_1based_inclusive":
        return load_canonical_tmds(
            config["tmd_annotation"]["path"],
            config["tmd_annotation"].get("columns"),
        )
    if system == "paper_legacy_genomic":
        return load_paper_legacy_tmds(
            config["tmd_annotation"]["path"], config["reference"]["annotation"]
        )
    raise ConfigurationError(f"Unsupported TMD coordinate system: {system}")


def run_sequence_feature_stage(config: dict[str, Any]) -> dict[str, Any]:
    output = Path(config["output"]["directory"])
    signal_path = output / "downstream" / "topology" / "baseline_tmd_scores.tsv"
    if not signal_path.is_file():
        raise ValidationFailed("Run the topology stage before sequence features.")
    options = config.get("optional_modules", {}).get("sequence_features", {})
    peptide = config.get("reference", {}).get("peptide_fasta")
    if not peptide or not Path(peptide).is_file():
        raise ValidationFailed("reference.peptide_fasta is required for sequence features.")
    flank = int(options.get("flank_length_aa", 10))
    if flank < 0:
        raise ConfigurationError("flank_length_aa cannot be negative.")
    key_type, protein_mapping, key_suffix_to_strip = _protein_mapping(options)
    sequences = _read_fasta(Path(peptide))
    records = {record.tmd_id: record for record in _load_tmd_records(config)}
    signals = pd.read_csv(signal_path, sep="\t")

    rows: list[dict[str, Any]] = []
    missing_mapping = missing_sequence = outside_sequence = 0
    for row in signals.itertuples(index=False):
        record = records.get(str(row.tmd_id))
        if record is None:
            continue
        mapping_key = record.gene_id if key_type == "gene_id" else record.transcript_id
        if key_suffix_to_strip and mapping_key.endswith(key_suffix_to_strip):
            mapping_key = mapping_key[: -len(key_suffix_to_strip)]
        protein_id = protein_mapping.get(mapping_key)
        if protein_id is None:
            missing_mapping += 1
            continue
        sequence = sequences.get(protein_id)
        if sequence is None:
            missing_sequence += 1
            continue
        start = (record.start_cds_nt - 1) // 3
        end = (record.end_cds_nt + 2) // 3
        if start < 0 or end > len(sequence) or end <= start:
            outside_sequence += 1
            continue
        tmd = sequence[start:end]
        n_flank = sequence[max(0, start - flank):start]
        c_flank = sequence[end:min(len(sequence), end + flank)]
        rows.append(
            {
                "tmd_id": record.tmd_id,
                "gene_id": record.gene_id,
                "transcript_id": record.transcript_id,
                "protein_id": protein_id,
                "mean_signal": float(row.mean_signal),
                "mean_background_count": float(row.mean_background_count),
                "annotated_tmd_count": int(row.annotated_tmd_count),
                "tmd_sequence": tmd,
                "n_flank_sequence": n_flank,
                "c_flank_sequence": c_flank,
                "tmd_length_aa": len(tmd),
                "tmd_gravy": float(np.mean([HYDROPATHY.get(x, np.nan) for x in tmd])),
                "tmd_hydrophobic_fraction": _fraction(tmd, "AILMFWVY"),
                "tmd_aromatic_fraction": _fraction(tmd, "FWY"),
                "tmd_bulky_fraction": _fraction(tmd, "FILMWVY"),
                "tmd_pro_gly_fraction": _fraction(tmd, "PG"),
                "tmd_charged_fraction": _fraction(tmd, "KRDE"),
                "n_flank_positive_fraction": _fraction(n_flank, "KR") if n_flank else np.nan,
                "n_flank_negative_fraction": _fraction(n_flank, "DE") if n_flank else np.nan,
                "n_flank_net_charge": _net_charge(n_flank),
                "c_flank_positive_fraction": _fraction(c_flank, "KR") if c_flank else np.nan,
                "c_flank_negative_fraction": _fraction(c_flank, "DE") if c_flank else np.nan,
                "c_flank_net_charge": _net_charge(c_flank),
            }
        )
    features = pd.DataFrame(rows)
    if features.empty:
        raise ValidationFailed("No baseline TMD could be matched to a protein sequence.")
    features["single_pass"] = (features["annotated_tmd_count"] == 1).astype(int)
    features["log2_background"] = np.log2(features["mean_background_count"] + 1.0)

    feature_names = [
        "tmd_length_aa", "tmd_gravy", "tmd_hydrophobic_fraction",
        "tmd_aromatic_fraction", "tmd_bulky_fraction", "tmd_pro_gly_fraction",
        "tmd_charged_fraction", "n_flank_positive_fraction",
        "n_flank_negative_fraction", "n_flank_net_charge",
        "c_flank_positive_fraction", "c_flank_negative_fraction",
        "c_flank_net_charge",
    ]
    association_rows: list[dict[str, Any]] = []
    for feature in feature_names:
        data = features.dropna(subset=[feature, "mean_signal"]).copy()
        standard_deviation = float(data[feature].std(ddof=0))
        if len(data) < 3 or standard_deviation == 0:
            continue
        data["feature_z"] = (data[feature] - data[feature].mean()) / standard_deviation
        design = sm.add_constant(data[["feature_z", "single_pass", "log2_background"]])
        fit = sm.OLS(data["mean_signal"], design).fit(
            cov_type="cluster", cov_kwds={"groups": data["gene_id"]}
        )
        association_rows.append(
            {
                "feature": feature,
                "n_tmds": len(data),
                "n_genes": data["gene_id"].nunique(),
                "coefficient_per_sd": float(fit.params["feature_z"]),
                "clustered_se": float(fit.bse["feature_z"]),
                "p_value": float(fit.pvalues["feature_z"]),
                "covariates": "single_pass;log2_background",
            }
        )
    associations = pd.DataFrame(association_rows)
    if associations.empty:
        raise ValidationFailed(
            "No sequence feature had sufficient non-missing variation for association testing."
        )
    associations["fdr"] = multipletests(
        associations["p_value"], method="fdr_bh"
    )[1]
    associations = associations.sort_values("fdr")

    joint_features = list(
        options.get(
            "joint_features",
            [
                "tmd_length_aa",
                "tmd_aromatic_fraction",
                "tmd_pro_gly_fraction",
                "tmd_charged_fraction",
                "c_flank_net_charge",
            ],
        )
    )
    unknown = sorted(set(joint_features) - set(feature_names))
    if unknown:
        raise ConfigurationError(f"Unknown joint sequence features: {unknown}")
    joint_data = features.dropna(subset=[*joint_features, "mean_signal"]).copy()
    if len(joint_data) < 3 or joint_data["gene_id"].nunique() < 2:
        raise ValidationFailed(
            "The joint sequence model requires at least three TMDs from two genes."
        )
    standardized: list[str] = []
    for feature in joint_features:
        sd = float(joint_data[feature].std(ddof=0))
        if sd > 0:
            column = f"{feature}_z"
            joint_data[column] = (joint_data[feature] - joint_data[feature].mean()) / sd
            standardized.append(column)
    joint_design = sm.add_constant(
        joint_data[[*standardized, "single_pass", "log2_background"]]
    )
    if not standardized:
        raise ValidationFailed(
            "None of the configured joint sequence features varies in the matched data."
        )
    joint_fit = sm.OLS(joint_data["mean_signal"], joint_design).fit(
        cov_type="cluster", cov_kwds={"groups": joint_data["gene_id"]}
    )
    joint_vif = {
        term: (
            np.nan
            if term == "const"
            else float(
                variance_inflation_factor(
                    joint_design.to_numpy(), joint_design.columns.get_loc(term)
                )
            )
        )
        for term in joint_design.columns
    }
    joint_rows = []
    for term in joint_fit.params.index:
        joint_rows.append(
            {
                "term": term,
                "coefficient": float(joint_fit.params[term]),
                "clustered_se": float(joint_fit.bse[term]),
                "p_value": float(joint_fit.pvalues[term]),
                "n_tmds": len(joint_data),
                "n_genes": joint_data["gene_id"].nunique(),
                "variance_inflation_factor": joint_vif[term],
            }
        )
    joint = pd.DataFrame(joint_rows)
    feature_mask = joint["term"].str.endswith("_z")
    joint["fdr_among_sequence_terms"] = np.nan
    if feature_mask.any():
        joint.loc[feature_mask, "fdr_among_sequence_terms"] = multipletests(
            joint.loc[feature_mask, "p_value"], method="fdr_bh"
        )[1]

    result_dir = output / "downstream" / "sequence_features"
    result_dir.mkdir(parents=True, exist_ok=True)
    feature_path = result_dir / "tmd_sequence_features.tsv"
    association_path = result_dir / "sequence_feature_associations.tsv"
    joint_path = result_dir / "joint_sequence_model.tsv"
    features.to_csv(feature_path, sep="\t", index=False)
    associations.to_csv(association_path, sep="\t", index=False)
    joint.to_csv(joint_path, sep="\t", index=False)
    return {
        "matched_tmds": int(len(features)),
        "matched_genes": int(features["gene_id"].nunique()),
        "missing_mapping": missing_mapping,
        "missing_sequence": missing_sequence,
        "outside_sequence": outside_sequence,
        "features": str(feature_path),
        "associations": str(association_path),
        "joint_model": str(joint_path),
    }
