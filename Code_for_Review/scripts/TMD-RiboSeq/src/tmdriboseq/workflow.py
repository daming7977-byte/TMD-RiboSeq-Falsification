"""Conservative stage orchestration for the supported development pipeline."""

from __future__ import annotations

import hashlib
import json
from copy import deepcopy
from pathlib import Path
from typing import Any

import yaml

from . import __version__
from .differential import run_offset_differential_stage
from .enrichment import run_enrichment_stage
from .errors import ConfigurationError
from .overlap import run_overlap_stage
from .profiles import load_samples, run_profile_stage
from .positional import run_positional_model_stage
from .provenance import create_run_manifest, file_sha256
from .qc import run_calibration_stage
from .qc_gate import run_qc_gate_stage
from .report import generate_html_report
from .sensitivity import run_target_background_sensitivity
from .sequence_features import run_sequence_feature_stage
from .topology import run_topology_stage


def _file_identity(
    path: str | Path, cache: dict[str, dict[str, Any]]
) -> dict[str, Any]:
    path = Path(path).resolve()
    key = str(path)
    if key in cache:
        return cache[key]
    stat = path.stat()
    identity = {
        "path": key,
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "sha256": file_sha256(path),
    }
    cache[key] = identity
    return identity


def _stage_signature(
    stage: str,
    config_path: str | Path,
    inputs: list[str | Path],
    upstream: str | None,
    identity_cache: dict[str, dict[str, Any]],
) -> str:
    config_path = Path(config_path)
    payload = {
        "stage": stage,
        "software_version": __version__,
        "config_sha256": hashlib.sha256(config_path.read_bytes()).hexdigest(),
        "inputs": [_file_identity(path, identity_cache) for path in inputs],
        "upstream": upstream,
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def _load_state(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {"schema_version": 1, "stages": {}}
    return json.loads(path.read_text(encoding="utf-8"))


def _save_state(path: Path, state: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(state, indent=2) + "\n", encoding="utf-8")


def _stage_is_reusable(
    state: dict[str, Any], stage: str, signature: str, outputs: list[Path]
) -> bool:
    record = state.get("stages", {}).get(stage, {})
    return record.get("signature") == signature and all(path.is_file() for path in outputs)


def _record_stage(
    state_path: Path,
    state: dict[str, Any],
    stage: str,
    signature: str,
    outputs: list[Path],
    summary: dict[str, Any],
) -> None:
    state.setdefault("stages", {})[stage] = {
        "signature": signature,
        "outputs": [str(path) for path in outputs],
        "summary": summary,
    }
    _save_state(state_path, state)


def run_workflow(
    config: dict[str, Any],
    config_path: str | Path,
    *,
    resume: bool = False,
    force: bool = False,
) -> dict[str, Any]:
    """Run validated stages with conservative signature-based reuse."""

    output = Path(config["output"]["directory"])
    state_path = output / "stage_state.json"
    if state_path.exists() and not resume and not force:
        raise ConfigurationError(
            f"A prior run exists at {output}. Use --resume to reuse matching stages "
            "or --force to intentionally recompute and overwrite outputs."
        )
    state = _load_state(state_path) if resume else {"schema_version": 1, "stages": {}}
    resolved = deepcopy(config)
    summaries: dict[str, Any] = {}
    upstream: str | None = None
    completed_stage = "validation"
    identity_cache: dict[str, dict[str, Any]] = {}
    common_inputs = [
        resolved["reference"]["annotation"],
        resolved["tmd_annotation"]["path"],
        resolved["samples"]["metadata"],
    ]
    bam_inputs = [sample.bam for sample in load_samples(resolved)]
    bai_inputs: list[Path] = []
    for bam in bam_inputs:
        bam = Path(bam)
        candidates = [Path(f"{bam}.bai"), bam.with_suffix(".bai")]
        bai_inputs.append(next(path for path in candidates if path.is_file()))
    raw_read_inputs = [*bam_inputs, *bai_inputs]

    if resolved["psite"].get("mode") == "calibrate":
        signature = _stage_signature(
            "calibrate",
            config_path,
            [*common_inputs, *raw_read_inputs],
            upstream,
            identity_cache,
        )
        offset_path = output / "qc" / "psite_offsets.tsv"
        if resume and _stage_is_reusable(state, "calibrate", signature, [offset_path]):
            summary = state["stages"]["calibrate"]["summary"]
            summary["resumed"] = True
        else:
            summary = run_calibration_stage(resolved)
            _record_stage(state_path, state, "calibrate", signature, [offset_path], summary)
        summaries["calibrate"] = summary
        completed_stage = "calibrate"
        resolved["psite"]["mode"] = "fixed_table"
        resolved["psite"]["offset_table"] = str(offset_path)
        upstream = signature
    elif resolved["psite"].get("mode") == "fixed_table":
        offset_path = Path(resolved["psite"]["offset_table"])
    else:
        raise ConfigurationError("psite.mode must be calibrate or fixed_table.")

    profile_inputs = [*common_inputs, *raw_read_inputs, offset_path]
    signature = _stage_signature(
        "profile", config_path, profile_inputs, upstream, identity_cache
    )
    profile_outputs = [
        output / "tmd_profiles.tsv.gz",
        output / "coverage" / "coverage_flow.tsv",
        output / "qc" / "baseline_equal_weight_metagene_summary.tsv",
    ]
    if resume and _stage_is_reusable(state, "profile", signature, profile_outputs):
        summary = state["stages"]["profile"]["summary"]
        summary["resumed"] = True
    else:
        summary = run_profile_stage(resolved)
        _record_stage(state_path, state, "profile", signature, profile_outputs, summary)
    summaries["profile"] = summary
    completed_stage = "profile"
    upstream = signature

    signature = _stage_signature(
        "qc_gate",
        config_path,
        [
            offset_path,
            output / "qc" / "sample_primary_signal_correlations.tsv",
            output / "qc" / "psite_assignment_summary.tsv",
        ],
        upstream,
        identity_cache,
    )
    qc_gate_outputs = [output / "qc" / "qc_gate.tsv"]
    if resume and _stage_is_reusable(
        state, "qc_gate", signature, qc_gate_outputs
    ):
        summary = state["stages"]["qc_gate"]["summary"]
        summary["resumed"] = True
    else:
        summary = run_qc_gate_stage(resolved)
        _record_stage(
            state_path, state, "qc_gate", signature, qc_gate_outputs, summary
        )
    summaries["qc_gate"] = summary
    completed_stage = "qc_gate"
    upstream = signature

    signature = _stage_signature(
        "differential", config_path, profile_outputs, upstream, identity_cache
    )
    differential_outputs = [output / "differential" / "primary_offset_model.tsv"]
    if resume and _stage_is_reusable(
        state, "differential", signature, differential_outputs
    ):
        summary = state["stages"]["differential"]["summary"]
        summary["resumed"] = True
    else:
        summary = run_offset_differential_stage(resolved)
        _record_stage(
            state_path, state, "differential", signature, differential_outputs, summary
        )
    summaries["differential"] = summary
    completed_stage = "differential"
    upstream = signature

    if bool(resolved["statistics"].get("positional_global_fdr", True)):
        signature = _stage_signature(
            "positional", config_path, profile_outputs, upstream, identity_cache
        )
        positional_outputs = [
            output / "positional_model" / "positional_model_results.tsv.gz",
            output / "positional_model" / "positional_model_summary.tsv",
        ]
        if resume and _stage_is_reusable(
            state, "positional", signature, positional_outputs
        ):
            summary = state["stages"]["positional"]["summary"]
            summary["resumed"] = True
        else:
            summary = run_positional_model_stage(resolved)
            _record_stage(
                state_path, state, "positional", signature, positional_outputs, summary
            )
        summaries["positional"] = summary
        completed_stage = "positional"
        upstream = signature

    if bool(resolved["statistics"].get("run_target_background_sensitivity", True)):
        signature = _stage_signature(
            "sensitivity",
            config_path,
            [*profile_outputs, *differential_outputs],
            upstream,
            identity_cache,
        )
        sensitivity_outputs = [
            output / "differential" / "target_background_sensitivity.tsv",
            output / "differential" / "model_sensitivity_comparison.tsv",
        ]
        if resume and _stage_is_reusable(
            state, "sensitivity", signature, sensitivity_outputs
        ):
            summary = state["stages"]["sensitivity"]["summary"]
            summary["resumed"] = True
        else:
            summary = run_target_background_sensitivity(resolved)
            _record_stage(
                state_path, state, "sensitivity", signature, sensitivity_outputs, summary
            )
        summaries["sensitivity"] = summary
        completed_stage = "sensitivity"
        upstream = signature

    topology_options = resolved.get("optional_modules", {}).get("topology", {})
    if bool(topology_options.get("enabled", False)):
        signature = _stage_signature(
            "topology", config_path, profile_outputs, upstream, identity_cache
        )
        topology_outputs = [
            output / "downstream" / "topology" / "baseline_tmd_scores.tsv",
            output / "downstream" / "topology" / "protein_summary.tsv",
            output / "downstream" / "topology" / "topology_summary.tsv",
        ]
        if resume and _stage_is_reusable(
            state, "topology", signature, topology_outputs
        ):
            summary = state["stages"]["topology"]["summary"]
            summary["resumed"] = True
        else:
            summary = run_topology_stage(resolved)
            _record_stage(
                state_path, state, "topology", signature, topology_outputs, summary
            )
        summaries["topology"] = summary
        completed_stage = "topology"
        upstream = signature

    sequence_options = resolved.get("optional_modules", {}).get(
        "sequence_features", {}
    )
    if bool(sequence_options.get("enabled", False)):
        sequence_inputs = [
            output / "downstream" / "topology" / "baseline_tmd_scores.tsv",
            resolved["reference"]["peptide_fasta"],
            sequence_options["protein_mapping"]["path"],
            resolved["tmd_annotation"]["path"],
        ]
        signature = _stage_signature(
            "sequence_features", config_path, sequence_inputs, upstream, identity_cache
        )
        sequence_outputs = [
            output / "downstream" / "sequence_features" / "tmd_sequence_features.tsv",
            output / "downstream" / "sequence_features" / "sequence_feature_associations.tsv",
            output / "downstream" / "sequence_features" / "joint_sequence_model.tsv",
        ]
        if resume and _stage_is_reusable(
            state, "sequence_features", signature, sequence_outputs
        ):
            summary = state["stages"]["sequence_features"]["summary"]
            summary["resumed"] = True
        else:
            summary = run_sequence_feature_stage(resolved)
            _record_stage(
                state_path,
                state,
                "sequence_features",
                signature,
                sequence_outputs,
                summary,
            )
        summaries["sequence_features"] = summary
        completed_stage = "sequence_features"
        upstream = signature

    overlap_options = resolved.get("optional_modules", {}).get("overlap", {})
    if bool(overlap_options.get("enabled", False)):
        overlap_inputs = [
            output / "differential" / "primary_offset_model.tsv",
            output / "downstream" / "topology" / "baseline_tmd_scores.tsv",
        ]
        signature = _stage_signature(
            "overlap", config_path, overlap_inputs, upstream, identity_cache
        )
        overlap_outputs = [
            output / "downstream" / "overlap" / "tmd_set_membership.tsv",
            output / "downstream" / "overlap" / "gene_set_membership.tsv",
            output / "downstream" / "overlap" / "pairwise_overlap.tsv",
        ]
        if resume and _stage_is_reusable(
            state, "overlap", signature, overlap_outputs
        ):
            summary = state["stages"]["overlap"]["summary"]
            summary["resumed"] = True
        else:
            summary = run_overlap_stage(resolved)
            _record_stage(
                state_path, state, "overlap", signature, overlap_outputs, summary
            )
        summaries["overlap"] = summary
        completed_stage = "overlap"
        upstream = signature

    enrichment_options = resolved.get("optional_modules", {}).get("enrichment", {})
    if bool(enrichment_options.get("enabled", False)):
        enrichment_inputs = [
            output / "differential" / "primary_offset_model.tsv",
            enrichment_options["mapping"]["path"],
        ]
        signature = _stage_signature(
            "enrichment", config_path, enrichment_inputs, upstream, identity_cache
        )
        enrichment_outputs = [
            output / "downstream" / "enrichment" / "tested_gene_universe.tsv",
            output / "downstream" / "enrichment" / "gene_set_enrichment.tsv",
        ]
        if resume and _stage_is_reusable(
            state, "enrichment", signature, enrichment_outputs
        ):
            summary = state["stages"]["enrichment"]["summary"]
            summary["resumed"] = True
        else:
            summary = run_enrichment_stage(resolved)
            _record_stage(
                state_path,
                state,
                "enrichment",
                signature,
                enrichment_outputs,
                summary,
            )
        summaries["enrichment"] = summary
        completed_stage = "enrichment"
        upstream = signature

    resolved_copy = output / "resolved_config.yaml"
    serializable = {key: value for key, value in resolved.items() if key != "_meta"}
    resolved_copy.parent.mkdir(parents=True, exist_ok=True)
    resolved_copy.write_text(yaml.safe_dump(serializable, sort_keys=False), encoding="utf-8")

    manifest_inputs = [*common_inputs, *raw_read_inputs, offset_path]
    if bool(sequence_options.get("enabled", False)):
        manifest_inputs.extend(
            [
                resolved["reference"]["peptide_fasta"],
                sequence_options["protein_mapping"]["path"],
            ]
        )
    if bool(enrichment_options.get("enabled", False)):
        manifest_inputs.append(enrichment_options["mapping"]["path"])
    manifest_hashes = {
        key: str(value["sha256"])
        for key, value in identity_cache.items()
        if "sha256" in value
    }
    manifest = create_run_manifest(
        output / "run_manifest.json",
        config_path,
        manifest_inputs,
        command=["tmdriboseq", "run", "--config", str(config_path)],
        extra={"completed_stage": completed_stage, "stage_summaries": summaries},
        precomputed_sha256=manifest_hashes,
    )
    report_path = generate_html_report(output)
    summaries["report"] = {"path": str(report_path)}
    return summaries
