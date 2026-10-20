"""Run provenance and file hashing."""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import platform
import shlex
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping

from . import __version__


def file_sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _package_version(name: str) -> str | None:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return None


def _command_output(command: list[str], cwd: Path) -> str | None:
    try:
        result = subprocess.run(
            command, cwd=cwd, check=True, capture_output=True, text=True
        )
    except (FileNotFoundError, subprocess.CalledProcessError):
        return None
    return result.stdout.strip()


def create_run_manifest(
    output_path: str | Path,
    config_path: str | Path,
    input_paths: Iterable[str | Path],
    command: list[str] | None = None,
    extra: dict[str, Any] | None = None,
    precomputed_sha256: Mapping[str, str] | None = None,
) -> dict[str, Any]:
    """Write a machine-readable provenance manifest for a workflow run."""

    output_path = Path(output_path)
    config_path = Path(config_path).resolve()
    root = config_path.parent
    inputs = []
    for value in input_paths:
        path = Path(value).resolve()
        checksum = None
        if precomputed_sha256 is not None:
            checksum = precomputed_sha256.get(str(path))
        inputs.append(
            {
                "path": str(path),
                "size_bytes": path.stat().st_size,
                "sha256": checksum or file_sha256(path),
            }
        )
    manifest: dict[str, Any] = {
        "schema_version": 1,
        "software": {"name": "tmdriboseq", "version": __version__},
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "command": shlex.join(command or sys.argv),
        "config": {
            "path": str(config_path),
            "sha256": file_sha256(config_path),
        },
        "inputs": inputs,
        "runtime": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "numpy": _package_version("numpy"),
            "pandas": _package_version("pandas"),
            "pysam": _package_version("pysam"),
            "pyyaml": _package_version("PyYAML"),
            "scipy": _package_version("scipy"),
            "statsmodels": _package_version("statsmodels"),
            "R": _command_output(
                ["Rscript", "-e", "cat(as.character(getRversion()))"], root
            ),
            "edgeR": _command_output(
                [
                    "Rscript",
                    "-e",
                    'cat(as.character(packageVersion("edgeR")))',
                ],
                root,
            ),
        },
        "git": {
            "commit": _command_output(["git", "rev-parse", "HEAD"], root),
            "dirty": bool(_command_output(["git", "status", "--porcelain"], root)),
        },
    }
    if extra:
        manifest["run"] = extra
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest
