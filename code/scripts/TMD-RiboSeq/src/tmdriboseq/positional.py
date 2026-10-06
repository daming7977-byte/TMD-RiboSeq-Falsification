"""Matched strict positional edgeR scan and global multiple testing."""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Any

from .differential import _edgeR_version, _r_script, _write_design
from .errors import ConfigurationError, ValidationFailed


def run_positional_model_stage(config: dict[str, Any]) -> dict[str, Any]:
    _edgeR_version()
    output = Path(config["output"]["directory"])
    profiles = output / "tmd_profiles.tsv.gz"
    if not profiles.is_file():
        raise ValidationFailed(f"Profile table not found: {profiles}")
    result_dir = output / "positional_model"
    design_path = result_dir / "resolved_design.tsv"
    factors = _write_design(design_path, config)
    references = config["design"].get("reference_levels", {})
    if any(name not in references for name in factors):
        raise ConfigurationError("design.reference_levels must define every factor.")
    results = result_dir / "positional_model_results.tsv.gz"
    summary = result_dir / "positional_model_summary.tsv"
    script = _r_script("fit_edger_positional_scan.R")
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
        str(float(config["statistics"]["fdr_threshold"])),
    ]
    completed = subprocess.run(command, check=False, capture_output=True, text=True)
    if completed.returncode != 0:
        raise ValidationFailed(
            "Positional edgeR scan failed.\n"
            f"STDOUT:\n{completed.stdout}\nSTDERR:\n{completed.stderr}"
        )
    return {"results": str(results), "summary": str(summary)}
