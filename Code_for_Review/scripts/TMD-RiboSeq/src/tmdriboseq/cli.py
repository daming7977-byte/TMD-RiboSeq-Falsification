"""Command-line interface for the staged framework implementation."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__
from .config import load_config
from .differential import run_offset_differential_stage
from .enrichment import run_enrichment_stage
from .errors import TMDRiboSeqError
from .paper_regression import check_corrected_outputs, check_manifest
from .overlap import run_overlap_stage
from .profiles import run_profile_stage
from .positional import run_positional_model_stage
from .provenance import create_run_manifest
from .qc import run_calibration_stage
from .qc_gate import run_qc_gate_stage
from .sensitivity import run_target_background_sensitivity
from .sequence_features import run_sequence_feature_stage
from .topology import run_topology_stage
from .report import generate_html_report
from .workflow import run_workflow
from .validation import format_report, validate_config


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="tmdriboseq",
        description="TMD-centered Ribo-seq analysis framework",
    )
    parser.add_argument("--version", action="version", version=__version__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    validate = subparsers.add_parser("validate", help="validate config and inputs")
    validate.add_argument("--config", required=True, type=Path)
    validate.add_argument("--json", action="store_true", help="print JSON report")

    paper = subparsers.add_parser(
        "paper-check", help="compare frozen paper metrics with expectations"
    )
    paper.add_argument(
        "--expected", type=Path, default=Path("paper/expected_results.yaml")
    )
    paper.add_argument(
        "--manifest", type=Path, default=Path("paper/reference_manifest.json")
    )

    corrected = subparsers.add_parser(
        "paper-check-corrected",
        help="compare corrected workflow outputs with the audited reference",
    )
    corrected.add_argument(
        "--expected", type=Path, default=Path("paper/corrected_results.yaml")
    )
    corrected.add_argument(
        "--output", type=Path, default=Path("paper/reproduction_results")
    )

    profile = subparsers.add_parser(
        "profile", help="count canonical TMD-centered target/background profiles"
    )
    profile.add_argument("--config", required=True, type=Path)

    calibrate = subparsers.add_parser(
        "calibrate", help="run CDS-periodicity-based P-site calibration and QC"
    )
    calibrate.add_argument("--config", required=True, type=Path)

    qc = subparsers.add_parser(
        "qc", help="evaluate P-site and replicate QC using declared warning gates"
    )
    qc.add_argument("--config", required=True, type=Path)

    differential = subparsers.add_parser(
        "differential", help="fit the validated edgeR 2 × 2 primary model"
    )
    differential.add_argument("--config", required=True, type=Path)

    sensitivity = subparsers.add_parser(
        "sensitivity", help="fit the target/background differential-usage model"
    )
    sensitivity.add_argument("--config", required=True, type=Path)

    positional = subparsers.add_parser(
        "positional", help="fit strict models at all scan centers with global FDR"
    )
    positional.add_argument("--config", required=True, type=Path)

    topology = subparsers.add_parser(
        "topology", help="run equal-weight protein/topology summaries"
    )
    topology.add_argument("--config", required=True, type=Path)

    sequence = subparsers.add_parser(
        "sequence-features", help="associate TMD sequence features with baseline signal"
    )
    sequence.add_argument("--config", required=True, type=Path)

    overlap = subparsers.add_parser(
        "overlap", help="compare basal, condition, and interaction signal sets"
    )
    overlap.add_argument("--config", required=True, type=Path)

    enrichment = subparsers.add_parser(
        "enrichment", help="test interaction genes in the analyzable gene universe"
    )
    enrichment.add_argument("--config", required=True, type=Path)

    report_parser = subparsers.add_parser(
        "report", help="generate a minimal HTML report from completed outputs"
    )
    report_parser.add_argument("--config", required=True, type=Path)

    run = subparsers.add_parser("run", help="run the supported end-to-end workflow")
    run.add_argument("--config", required=True, type=Path)
    run_mode = run.add_mutually_exclusive_group()
    run_mode.add_argument("--resume", action="store_true")
    run_mode.add_argument("--force", action="store_true")
    return parser


def _run_validate(args: argparse.Namespace) -> int:
    config = load_config(args.config)
    report = validate_config(config)
    print(report.to_json() if args.json else format_report(report))
    return 0 if report.passed else 2


def _run_paper_check(args: argparse.Namespace) -> int:
    mismatches = check_manifest(args.expected, args.manifest)
    if not mismatches:
        print("Paper regression reference: PASS")
        return 0
    print(f"Paper regression reference: FAIL ({len(mismatches)} mismatch(es))")
    for item in mismatches:
        print(f"ERROR: {item.path}: expected={item.expected!r}, observed={item.observed!r}")
    return 2


def _run_corrected_paper_check(args: argparse.Namespace) -> int:
    mismatches = check_corrected_outputs(args.expected, args.output)
    if not mismatches:
        print("Corrected paper workflow reference: PASS")
        return 0
    print(f"Corrected paper workflow reference: FAIL ({len(mismatches)} mismatch(es))")
    for item in mismatches:
        print(f"ERROR: {item.path}: expected={item.expected!r}, observed={item.observed!r}")
    return 2


def _run_profile(args: argparse.Namespace) -> int:
    config = load_config(args.config)
    report = validate_config(config)
    if not report.passed:
        print(format_report(report))
        return 2
    summary = run_profile_stage(config)
    output_directory = Path(config["output"]["directory"])
    input_paths = [
        config["reference"]["annotation"],
        config["tmd_annotation"]["path"],
        config["samples"]["metadata"],
        config["psite"]["offset_table"],
    ]
    create_run_manifest(
        output_directory / "run_manifest.json",
        args.config,
        input_paths,
        command=["tmdriboseq", "profile", "--config", str(args.config)],
        extra={"completed_stage": "profile", **summary},
    )
    print("Profile stage: PASS")
    for key, value in summary.items():
        print(f"  {key}: {value}")
    return 0


def _run_calibrate(args: argparse.Namespace) -> int:
    config = load_config(args.config)
    report = validate_config(config)
    if not report.passed:
        print(format_report(report))
        return 2
    summary = run_calibration_stage(config)
    print("Calibration stage: PASS")
    for key, value in summary.items():
        print(f"  {key}: {value}")
    return 0


def _run_qc(args: argparse.Namespace) -> int:
    config = load_config(args.config)
    report = validate_config(config)
    if not report.passed:
        print(format_report(report))
        return 2
    summary = run_qc_gate_stage(config)
    print("QC gate stage: PASS")
    for key, value in summary.items():
        print(f"  {key}: {value}")
    return 0


def _run_differential(args: argparse.Namespace) -> int:
    config = load_config(args.config)
    report = validate_config(config)
    if not report.passed:
        print(format_report(report))
        return 2
    summary = run_offset_differential_stage(config)
    print("Differential stage: PASS")
    for key, value in summary.items():
        print(f"  {key}: {value}")
    return 0


def _run_sensitivity(args: argparse.Namespace) -> int:
    config = load_config(args.config)
    report = validate_config(config)
    if not report.passed:
        print(format_report(report))
        return 2
    summary = run_target_background_sensitivity(config)
    print("Sensitivity stage: PASS")
    for key, value in summary.items():
        print(f"  {key}: {value}")
    return 0


def _run_positional(args: argparse.Namespace) -> int:
    config = load_config(args.config)
    report = validate_config(config)
    if not report.passed:
        print(format_report(report))
        return 2
    summary = run_positional_model_stage(config)
    print("Positional model stage: PASS")
    for key, value in summary.items():
        print(f"  {key}: {value}")
    return 0


def _run_topology(args: argparse.Namespace) -> int:
    config = load_config(args.config)
    report = validate_config(config)
    if not report.passed:
        print(format_report(report))
        return 2
    summary = run_topology_stage(config)
    print("Topology stage: PASS")
    for key, value in summary.items():
        print(f"  {key}: {value}")
    return 0


def _run_sequence_features(args: argparse.Namespace) -> int:
    config = load_config(args.config)
    report = validate_config(config)
    if not report.passed:
        print(format_report(report))
        return 2
    summary = run_sequence_feature_stage(config)
    print("Sequence-feature stage: PASS")
    for key, value in summary.items():
        print(f"  {key}: {value}")
    return 0


def _run_overlap(args: argparse.Namespace) -> int:
    config = load_config(args.config)
    report = validate_config(config)
    if not report.passed:
        print(format_report(report))
        return 2
    summary = run_overlap_stage(config)
    print("Overlap stage: PASS")
    for key, value in summary.items():
        print(f"  {key}: {value}")
    return 0


def _run_enrichment(args: argparse.Namespace) -> int:
    config = load_config(args.config)
    report = validate_config(config)
    if not report.passed:
        print(format_report(report))
        return 2
    summary = run_enrichment_stage(config)
    print("Enrichment stage: PASS")
    for key, value in summary.items():
        print(f"  {key}: {value}")
    return 0


def _run_report(args: argparse.Namespace) -> int:
    config = load_config(args.config)
    path = generate_html_report(config["output"]["directory"])
    print(f"Report: {path}")
    return 0


def _run_all(args: argparse.Namespace) -> int:
    config = load_config(args.config)
    report = validate_config(config)
    if not report.passed:
        print(format_report(report))
        return 2
    summaries = run_workflow(
        config, args.config, resume=bool(args.resume), force=bool(args.force)
    )
    print("Workflow: PASS")
    for stage, summary in summaries.items():
        print(f"  {stage}: {summary}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "validate":
            return _run_validate(args)
        if args.command == "paper-check":
            return _run_paper_check(args)
        if args.command == "paper-check-corrected":
            return _run_corrected_paper_check(args)
        if args.command == "profile":
            return _run_profile(args)
        if args.command == "calibrate":
            return _run_calibrate(args)
        if args.command == "qc":
            return _run_qc(args)
        if args.command == "differential":
            return _run_differential(args)
        if args.command == "sensitivity":
            return _run_sensitivity(args)
        if args.command == "positional":
            return _run_positional(args)
        if args.command == "topology":
            return _run_topology(args)
        if args.command == "sequence-features":
            return _run_sequence_features(args)
        if args.command == "overlap":
            return _run_overlap(args)
        if args.command == "enrichment":
            return _run_enrichment(args)
        if args.command == "report":
            return _run_report(args)
        if args.command == "run":
            return _run_all(args)
    except TMDRiboSeqError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    parser.error(f"Unsupported command: {args.command}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
