"""Coverage eligibility summaries and strict comparison intersections."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable


def summarize_coverage(
    profile_rows: Iterable[dict[str, Any]],
    primary_center_nt: int,
    model_role: str,
    factor_names: list[str],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    primary = [
        row
        for row in profile_rows
        if int(row["center_nt"]) == primary_center_nt and row["role"] == model_role
    ]
    sample_tmd_pass: dict[str, dict[str, bool]] = defaultdict(dict)
    sample_factors: dict[str, tuple[str, ...]] = {}
    all_tmds: set[str] = set()
    for row in primary:
        sample = str(row["sample_id"])
        tmd_id = str(row["tmd_id"])
        sample_tmd_pass[sample][tmd_id] = bool(row["coverage_pass"])
        sample_factors[sample] = tuple(str(row[name]) for name in factor_names)
        all_tmds.add(tmd_id)

    samples = sorted(sample_tmd_pass)
    strict = {
        tmd_id
        for tmd_id in all_tmds
        if samples and all(sample_tmd_pass[sample].get(tmd_id, False) for sample in samples)
    }

    flow: list[dict[str, Any]] = []
    for sample in samples:
        flow.append(
            {
                "scope": "sample",
                "group": sample,
                "n_required_samples": 1,
                "eligible_tmds": sum(sample_tmd_pass[sample].values()),
                "total_tmds": len(all_tmds),
            }
        )

    cells: dict[tuple[str, ...], list[str]] = defaultdict(list)
    for sample, factors in sample_factors.items():
        cells[factors].append(sample)
    for cell, cell_samples in sorted(cells.items()):
        eligible = sum(
            all(sample_tmd_pass[sample].get(tmd_id, False) for sample in cell_samples)
            for tmd_id in all_tmds
        )
        flow.append(
            {
                "scope": "design_cell",
                "group": " | ".join(cell),
                "n_required_samples": len(cell_samples),
                "eligible_tmds": eligible,
                "total_tmds": len(all_tmds),
            }
        )
    flow.append(
        {
            "scope": "strict_all_model_samples",
            "group": "all model samples",
            "n_required_samples": len(samples),
            "eligible_tmds": len(strict),
            "total_tmds": len(all_tmds),
        }
    )

    membership: list[dict[str, Any]] = []
    for tmd_id in sorted(all_tmds):
        row: dict[str, Any] = {
            "tmd_id": tmd_id,
            "strict_all_model_samples": tmd_id in strict,
        }
        for sample in samples:
            row[f"pass__{sample}"] = sample_tmd_pass[sample].get(tmd_id, False)
        membership.append(row)
    return flow, membership


def write_tsv(path: str | Path, rows: list[dict[str, Any]]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"Cannot write empty coverage table: {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

