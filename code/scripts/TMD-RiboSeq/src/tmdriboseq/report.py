"""Minimal dependency-free HTML run report."""

from __future__ import annotations

import csv
import html
import json
from pathlib import Path
from typing import Any


def _read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _table(rows: list[dict[str, Any]], maximum_rows: int = 30) -> str:
    if not rows:
        return '<p class="not-run">Not run or no rows available.</p>'
    columns = list(rows[0])
    header = "".join(f"<th>{html.escape(str(column))}</th>" for column in columns)
    body_rows = []
    for row in rows[:maximum_rows]:
        cells = "".join(
            f"<td>{html.escape(str(row.get(column, '')))}</td>" for column in columns
        )
        body_rows.append(f"<tr>{cells}</tr>")
    note = ""
    if len(rows) > maximum_rows:
        note = f"<p>Showing {maximum_rows} of {len(rows)} rows.</p>"
    return f"{note}<div class='table-wrap'><table><thead><tr>{header}</tr></thead><tbody>{''.join(body_rows)}</tbody></table></div>"


def _link(output: Path, target: Path, label: str) -> str:
    if not target.is_file():
        return f"<li>{html.escape(label)}: not generated</li>"
    relative = target.relative_to(output)
    return f'<li><a href="{html.escape(str(relative))}">{html.escape(label)}</a></li>'


def generate_html_report(output_directory: str | Path) -> Path:
    output = Path(output_directory)
    manifest_path = output / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.is_file() else {}
    coverage = _read_tsv(output / "coverage" / "coverage_flow.tsv")
    mapping = _read_tsv(output / "qc" / "mapping_summary.tsv")
    if not mapping:
        mapping = _read_tsv(output / "qc" / "psite_assignment_summary.tsv")
    offsets = _read_tsv(output / "qc" / "psite_offsets.tsv")
    cds_start = _read_tsv(output / "qc" / "cds_start_metagene.tsv")
    correlations = _read_tsv(
        output / "qc" / "sample_primary_signal_correlations.tsv"
    )
    qc_gate = _read_tsv(output / "qc" / "qc_gate.tsv")
    baseline_metagene = _read_tsv(
        output / "qc" / "baseline_equal_weight_metagene_summary.tsv"
    )
    design = _read_tsv(output / "differential" / "resolved_design.tsv")
    primary_summary = _read_tsv(
        output / "differential" / "primary_offset_model_summary.tsv"
    )
    positional_summary = _read_tsv(
        output / "positional_model" / "positional_model_summary.tsv"
    )
    sensitivity = _read_tsv(
        output / "differential" / "model_sensitivity_comparison.tsv"
    )
    topology = _read_tsv(
        output / "downstream" / "topology" / "topology_summary.tsv"
    )
    sequence = _read_tsv(
        output / "downstream" / "sequence_features" / "sequence_feature_associations.tsv"
    )
    overlap = _read_tsv(
        output / "downstream" / "overlap" / "pairwise_overlap.tsv"
    )
    enrichment = _read_tsv(
        output / "downstream" / "enrichment" / "gene_set_enrichment.tsv"
    )
    warnings = [
        row for row in offsets if str(row.get("qc_status", "")).upper() == "WARN"
    ]

    files = [
        _link(output, output / "tmd_profiles.tsv.gz", "TMD positional profiles"),
        _link(output, output / "coverage" / "coverage_flow.tsv", "Coverage flow"),
        _link(output, output / "qc" / "psite_offsets.tsv", "P-site offsets"),
        _link(
            output,
            output / "qc" / "cds_start_metagene.tsv",
            "CDS-start metagene table",
        ),
        _link(
            output,
            output / "qc" / "sample_primary_signal_correlations.tsv",
            "Sample signal correlations",
        ),
        _link(output, output / "qc" / "qc_gate.tsv", "QC gate decisions"),
        _link(
            output,
            output / "qc" / "baseline_equal_weight_metagene_summary.tsv",
            "Equal-weight baseline metagene",
        ),
        _link(
            output,
            output / "differential" / "primary_offset_model.tsv",
            "Primary edgeR results",
        ),
        _link(
            output,
            output / "differential" / "target_background_sensitivity.tsv",
            "Target/background sensitivity results",
        ),
        _link(
            output,
            output / "positional_model" / "positional_model_results.tsv.gz",
            "All-position edgeR results",
        ),
        _link(
            output,
            output / "positional_model" / "positional_model_summary.tsv",
            "All-position summary",
        ),
        _link(
            output,
            output / "downstream" / "topology" / "protein_summary.tsv",
            "Protein-level topology summary",
        ),
        _link(
            output,
            output / "downstream" / "sequence_features" / "sequence_feature_associations.tsv",
            "Sequence-feature associations",
        ),
        _link(
            output,
            output / "downstream" / "overlap" / "pairwise_overlap.tsv",
            "Pairwise set overlaps",
        ),
        _link(
            output,
            output / "downstream" / "enrichment" / "gene_set_enrichment.tsv",
            "Gene-set enrichment",
        ),
    ]
    software = manifest.get("software", {})
    run = manifest.get("run", {})
    document = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>TMD-RiboSeq run report</title>
<style>
body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; margin: 2rem auto; max-width: 1200px; padding: 0 1rem; color: #1f2933; }}
h1, h2 {{ color: #17324d; }}
.summary {{ background: #eef5fb; border-left: 5px solid #2b6f9f; padding: 1rem; }}
.warning {{ background: #fff6df; border-left: 5px solid #c88400; padding: 1rem; }}
.not-run {{ color: #6b7280; font-style: italic; }}
.table-wrap {{ overflow-x: auto; }}
table {{ border-collapse: collapse; width: 100%; font-size: 0.88rem; }}
th, td {{ border: 1px solid #d4dbe2; padding: 0.4rem; text-align: left; }}
th {{ background: #edf1f5; }}
code {{ background: #f2f4f6; padding: 0.1rem 0.25rem; }}
</style>
</head>
<body>
<h1>TMD-centered Ribo-seq analysis</h1>
<div class="summary">
<p><strong>Software:</strong> {html.escape(str(software.get('name', 'tmdriboseq')))} {html.escape(str(software.get('version', 'development')))}</p>
<p><strong>Completed stage:</strong> {html.escape(str(run.get('completed_stage', 'unknown')))}</p>
<p>This report summarizes local ribosome occupancy. It does not by itself prove mechanistic translational pausing.</p>
</div>
<h2>Warnings and limitations</h2>
{('<div class="warning">' + _table(warnings) + '</div>') if warnings else '<p>No offset-level QC warnings were recorded.</p>'}
<h2>QC and mapping</h2>
{_table(mapping)}
<h2>Resolved statistical design</h2>
{_table(design)}
<h2>P-site calibration</h2>
{_table(offsets)}
<h2>QC gate decisions</h2>
{_table(qc_gate)}
<h2>CDS-start metagene</h2>
{_table(cds_start)}
<h2>Sample correlations</h2>
{_table(correlations)}
<h2>Equal-weight baseline TMD metagene</h2>
{_table(baseline_metagene)}
<h2>Coverage flow</h2>
{_table(coverage)}
<h2>Primary differential model</h2>
{_table(primary_summary)}
<h2>All-position model</h2>
{_table(positional_summary)}
<h2>Target/background sensitivity</h2>
{_table(sensitivity)}
<h2>Protein/topology analysis</h2>
{_table(topology)}
<h2>Sequence-feature associations</h2>
{_table(sequence)}
<h2>Dataset overlap analysis</h2>
{_table(overlap)}
<h2>Gene-set enrichment</h2>
{_table(enrichment)}
<h2>Machine-readable outputs</h2>
<ul>{''.join(files)}</ul>
</body>
</html>
"""
    report_path = output / "report.html"
    report_path.write_text(document, encoding="utf-8")
    return report_path
