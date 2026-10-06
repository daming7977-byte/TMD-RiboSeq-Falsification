"""Python/R interface for the validated edgeR 2 × 2 primary model."""

from __future__ import annotations

import csv
import shutil
import subprocess
from pathlib import Path
from typing import Any

from .errors import ConfigurationError, ValidationFailed
from .profiles import load_samples


def _r_script(name: str) -> Path:
    """Return an R backend shipped inside the installed Python package."""

    path = Path(__file__).with_name("r_scripts") / name
    if not path.is_file():
        raise ConfigurationError(f"Packaged R backend is missing: {path}")
    return path


def _write_design(path: Path, config: dict[str, Any]) -> list[str]:
    factors = [str(value) for value in config["design"]["factors"]]
    if config["design"]["type"] != "factorial_2x2" or len(factors) != 2:
        raise ConfigurationError(
            "The current edgeR backend supports validated factorial_2x2 designs only."
        )
    model_role = str(config["samples"].get("model_role", "model"))
    samples = [sample for sample in load_samples(config) if sample.role == model_role]
    if not samples:
        raise ValidationFailed("No model samples were found.")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample_id", *factors], delimiter="\t")
        writer.writeheader()
        for sample in samples:
            writer.writerow(
                {"sample_id": sample.sample_id, **{name: sample.factors[name] for name in factors}}
            )
    return factors


def _edgeR_version() -> str:
    if shutil.which("Rscript") is None:
        raise ConfigurationError("Rscript is required but was not found on PATH.")
    result = subprocess.run(
        [
            "Rscript",
            "-e",
            'if (!requireNamespace("edgeR", quietly=TRUE)) quit(status=2); '
            'cat(as.character(packageVersion("edgeR")))',
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise ConfigurationError("The required Bioconductor package edgeR is not available.")
    return result.stdout.strip()


def run_offset_differential_stage(config: dict[str, Any]) -> dict[str, Any]:
    """Fit the primary target-count/background-offset interaction model."""

    edgeR_version = _edgeR_version()
    output = Path(config["output"]["directory"])
    profiles = output / "tmd_profiles.tsv.gz"
    if not profiles.is_file():
        raise ValidationFailed(
            f"Profile table not found: {profiles}. Run the profile stage first."
        )
    differential_dir = output / "differential"
    design_path = differential_dir / "resolved_design.tsv"
    factors = _write_design(design_path, config)
    references = config["design"].get("reference_levels", {})
    if any(name not in references for name in factors):
        raise ConfigurationError("design.reference_levels must define every factor.")
    results = differential_dir / "primary_offset_model.tsv"
    summary = differential_dir / "primary_offset_model_summary.tsv"
    script = _r_script("fit_edger_offset.R")
    command = [
        "Rscript",
        str(script),
        str(profiles),
        str(design_path),
        str(results),
        str(summary),
        factors[0],
        factors[1],
        str(references[factors[0]]),
        str(references[factors[1]]),
        str(int(config["profiles"]["primary_center_nt"])),
        str(float(config["statistics"]["fdr_threshold"])),
    ]
    completed = subprocess.run(command, check=False, capture_output=True, text=True)
    if completed.returncode != 0:
        raise ValidationFailed(
            "edgeR primary model failed.\n"
            f"STDOUT:\n{completed.stdout}\nSTDERR:\n{completed.stderr}"
        )
    return {
        "backend": "edgeR_ql_offset",
        "edgeR_version": edgeR_version,
        "results": str(results),
        "summary": str(summary),
        "stdout": completed.stdout.strip(),
    }
