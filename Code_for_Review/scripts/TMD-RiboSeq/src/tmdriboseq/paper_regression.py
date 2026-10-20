"""Regression checks against the frozen paper reference manifest."""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml
import pandas as pd


@dataclass(frozen=True)
class Mismatch:
    path: str
    expected: Any
    observed: Any


def _compare(
    expected: Any,
    observed: Any,
    path: str,
    mismatches: list[Mismatch],
    relative_tolerance: float,
    absolute_tolerance: float,
) -> None:
    if isinstance(expected, dict):
        if not isinstance(observed, dict):
            mismatches.append(Mismatch(path, expected, observed))
            return
        for key, value in expected.items():
            if key not in observed:
                mismatches.append(Mismatch(f"{path}.{key}", value, "<missing>"))
            else:
                _compare(
                    value,
                    observed[key],
                    f"{path}.{key}",
                    mismatches,
                    relative_tolerance,
                    absolute_tolerance,
                )
        return
    if isinstance(expected, list):
        if expected != observed:
            mismatches.append(Mismatch(path, expected, observed))
        return
    if isinstance(expected, (int, float)) and not isinstance(expected, bool):
        try:
            is_close = math.isclose(
                float(expected),
                float(observed),
                rel_tol=relative_tolerance,
                abs_tol=absolute_tolerance,
            )
        except (TypeError, ValueError):
            is_close = False
        if not is_close:
            mismatches.append(Mismatch(path, expected, observed))
        return
    if expected != observed:
        mismatches.append(Mismatch(path, expected, observed))


def check_manifest(
    expected_path: str | Path,
    manifest_path: str | Path,
    relative_tolerance: float = 1e-8,
    absolute_tolerance: float = 1e-10,
) -> list[Mismatch]:
    expected_document = yaml.safe_load(Path(expected_path).read_text(encoding="utf-8"))
    manifest = json.loads(Path(manifest_path).read_text(encoding="utf-8"))
    expected = dict(expected_document["expected"])
    expected["psite_consensus_offsets"] = expected_document[
        "psite_consensus_offsets"
    ]
    observed = manifest["observed"]
    mismatches: list[Mismatch] = []
    _compare(
        expected,
        observed,
        "expected",
        mismatches,
        relative_tolerance,
        absolute_tolerance,
    )
    return mismatches


def check_corrected_outputs(
    expected_path: str | Path,
    output_directory: str | Path,
    relative_tolerance: float = 1e-8,
    absolute_tolerance: float = 1e-10,
) -> list[Mismatch]:
    """Compare corrected workflow outputs with the audited reference state."""

    expected_document = yaml.safe_load(Path(expected_path).read_text(encoding="utf-8"))
    expected = expected_document["corrected"]
    output = Path(output_directory)
    primary_summary_table = pd.read_csv(
        output / "differential" / "primary_offset_model_summary.tsv", sep="\t"
    )
    primary_summary = dict(
        zip(primary_summary_table["metric"], primary_summary_table["value"])
    )
    primary = pd.read_csv(
        output / "differential" / "primary_offset_model.tsv", sep="\t"
    )
    run_manifest = json.loads((output / "run_manifest.json").read_text(encoding="utf-8"))
    run_block = run_manifest["run"]
    profile_summary = run_block.get("stage_summaries", {}).get("profile", run_block)
    threshold = float(primary_summary["fdr_threshold"])
    significant = primary[primary["interaction_fdr"] < threshold]
    observed: dict[str, Any] = {
        "annotated_tmds": int(profile_summary["tmds"]),
        "strict_coverage_eligible_tmds": int(primary_summary["coverage_eligible_tmds"]),
        "primary_tested_tmds": int(primary_summary["tested_tmds"]),
        "wt_stress_fdr_lt_0_05": int(primary_summary["reference_b_significant"]),
        "mutant_stress_fdr_lt_0_05": int(primary_summary["alternate_b_significant"]),
        "interaction_fdr_lt_0_05": int(primary_summary["interaction_significant"]),
        "interaction_gene_count": int(significant["gene_id"].nunique()),
        "interaction_genes": sorted(significant["gene_id"].astype(str).unique()),
        "interaction_tmd_ids": sorted(significant["tmd_id"].astype(str).unique()),
    }
    mismatches: list[Mismatch] = []
    _compare(
        expected,
        observed,
        "corrected",
        mismatches,
        relative_tolerance,
        absolute_tolerance,
    )
    return mismatches
