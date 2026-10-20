"""Transparent WARN-first QC gating for supported workflow runs."""

from __future__ import annotations

from itertools import combinations
from pathlib import Path
from typing import Any

import pandas as pd

from .errors import ValidationFailed
from .profiles import load_samples


def run_qc_gate_stage(config: dict[str, Any]) -> dict[str, Any]:
    """Write QC decisions without silently excluding biological samples."""

    output = Path(config["output"]["directory"])
    offset_path = Path(config["psite"]["offset_table"])
    correlation_path = output / "qc" / "sample_primary_signal_correlations.tsv"
    assignment_path = output / "qc" / "psite_assignment_summary.tsv"
    if not offset_path.is_file() or not correlation_path.is_file():
        raise ValidationFailed("P-site offsets and sample correlations are required for QC gating.")

    options = config.get("qc", {})
    frame_warn = float(
        options.get(
            "minimum_frame0_fraction_warn",
            config.get("psite", {}).get("warning_frame0_fraction", 0.5),
        )
    )
    correlation_warn = float(options.get("minimum_replicate_correlation_warn", 0.3))
    common_warn = int(options.get("minimum_common_tmds_warn", 100))

    offsets = pd.read_csv(offset_path, sep="\t")
    sample_column = "sample_id" if "sample_id" in offsets.columns else "accession"
    rows: list[dict[str, Any]] = []
    for row in offsets.itertuples(index=False):
        sample = str(getattr(row, sample_column))
        read_length = int(getattr(row, "read_length"))
        fraction = getattr(row, "frame0_fraction", None)
        if fraction is None or pd.isna(fraction):
            status = "NOT_EVALUATED"
            reason = "frame0_fraction was not supplied in the offset table"
        elif float(fraction) < frame_warn:
            status = "WARN"
            reason = f"frame0_fraction {float(fraction):.3f} < {frame_warn:.3f}"
        else:
            status = "PASS"
            reason = ""
        rows.append(
            {
                "scope": "sample_read_length",
                "sample_id": sample,
                "comparison": f"read_length={read_length}",
                "metric": "frame0_fraction",
                "value": fraction,
                "threshold": frame_warn,
                "status": status,
                "reason": reason,
            }
        )

    samples = load_samples(config)
    model_role = str(config["samples"].get("model_role", "model"))
    factor_names = [str(value) for value in config["design"]["factors"]]
    by_cell: dict[tuple[str, ...], list[str]] = {}
    for sample in samples:
        if sample.role != model_role:
            continue
        cell = tuple(sample.factors[name] for name in factor_names)
        by_cell.setdefault(cell, []).append(sample.sample_id)
    correlations = pd.read_csv(correlation_path, sep="\t")
    correlation_lookup = {
        (str(row.sample_a), str(row.sample_b)): (int(row.n_common_tmds), row.pearson_r)
        for row in correlations.itertuples(index=False)
    }
    for cell, sample_ids in sorted(by_cell.items()):
        for left, right in combinations(sorted(sample_ids), 2):
            common, correlation = correlation_lookup.get((left, right), (0, None))
            status = "PASS"
            reasons: list[str] = []
            if common < common_warn:
                status = "WARN"
                reasons.append(f"common TMDs {common} < {common_warn}")
            if correlation is None or pd.isna(correlation):
                status = "WARN"
                reasons.append("correlation unavailable")
            elif float(correlation) < correlation_warn:
                status = "WARN"
                reasons.append(
                    f"Pearson r {float(correlation):.3f} < {correlation_warn:.3f}"
                )
            rows.append(
                {
                    "scope": "replicate_pair",
                    "sample_id": f"{left};{right}",
                    "comparison": " | ".join(cell),
                    "metric": "primary_signal_pearson_r",
                    "value": correlation,
                    "threshold": correlation_warn,
                    "status": status,
                    "reason": "; ".join(reasons),
                }
            )

    if assignment_path.is_file():
        assignment = pd.read_csv(assignment_path, sep="\t")
        for row in assignment.itertuples(index=False):
            eligible = int(row.eligible_alignments)
            assigned = int(row.assigned_psites)
            rows.append(
                {
                    "scope": "sample",
                    "sample_id": str(row.sample_id),
                    "comparison": "assignment",
                    "metric": "assigned_fraction_of_eligible",
                    "value": assigned / eligible if eligible else 0.0,
                    "threshold": "not_predeclared",
                    "status": "NOT_EVALUATED",
                    "reason": "reported for audit; no universal exclusion threshold is imposed",
                }
            )

    gate = pd.DataFrame(rows)
    result_path = output / "qc" / "qc_gate.tsv"
    result_path.parent.mkdir(parents=True, exist_ok=True)
    gate.to_csv(result_path, sep="\t", index=False)
    warn_count = int((gate["status"] == "WARN").sum())
    fail_count = int((gate["status"] == "FAIL").sum())
    not_evaluated_count = int((gate["status"] == "NOT_EVALUATED").sum())
    overall = "FAIL" if fail_count else "WARN" if warn_count else "PASS"
    if fail_count or (warn_count and bool(options.get("fail_on_warning", False))):
        raise ValidationFailed(
            f"QC gate did not permit formal inference: {warn_count} WARN and "
            f"{fail_count} FAIL result(s). Review {result_path}."
        )
    return {
        "overall_status": overall,
        "warnings": warn_count,
        "failures": fail_count,
        "not_evaluated": not_evaluated_count,
        "fail_on_warning": bool(options.get("fail_on_warning", False)),
        "table": str(result_path),
    }
