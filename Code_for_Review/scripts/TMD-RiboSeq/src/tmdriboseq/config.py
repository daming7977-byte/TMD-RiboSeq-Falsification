"""Configuration loading and path resolution."""

from __future__ import annotations

import os
from copy import deepcopy
from pathlib import Path
from typing import Any

import yaml

from .errors import ConfigurationError


SUPPORTED_SCHEMA_VERSION = 1


def _expand_environment(value: Any) -> Any:
    if isinstance(value, str):
        return os.path.expanduser(os.path.expandvars(value))
    if isinstance(value, list):
        return [_expand_environment(item) for item in value]
    if isinstance(value, dict):
        return {key: _expand_environment(item) for key, item in value.items()}
    return value


def _resolve_path(base: Path, value: str | None) -> str | None:
    if value in (None, ""):
        return value
    path = Path(value)
    if not path.is_absolute():
        path = base / path
    return str(path.resolve(strict=False))


def _resolve_known_paths(config: dict[str, Any], config_path: Path) -> None:
    base = config_path.parent
    path_fields = (
        ("reference", "annotation"),
        ("reference", "genome_fasta"),
        ("reference", "peptide_fasta"),
        ("tmd_annotation", "path"),
        ("samples", "metadata"),
        ("psite", "offset_table"),
        ("output", "directory"),
    )
    for section, field in path_fields:
        block = config.get(section)
        if isinstance(block, dict) and field in block:
            block[field] = _resolve_path(base, block.get(field))
    optional = config.get("optional_modules")
    if isinstance(optional, dict):
        sequence = optional.get("sequence_features")
        if isinstance(sequence, dict):
            mapping = sequence.get("protein_mapping")
            if isinstance(mapping, dict) and "path" in mapping:
                mapping["path"] = _resolve_path(base, mapping.get("path"))
        enrichment = optional.get("enrichment")
        if isinstance(enrichment, dict):
            mapping = enrichment.get("mapping")
            if isinstance(mapping, dict) and "path" in mapping:
                mapping["path"] = _resolve_path(base, mapping.get("path"))


def load_config(path: str | Path) -> dict[str, Any]:
    """Load YAML, expand environment variables, and resolve known paths."""

    config_path = Path(path).expanduser().resolve(strict=False)
    if not config_path.is_file():
        raise ConfigurationError(f"Configuration file does not exist: {config_path}")

    try:
        document = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        raise ConfigurationError(f"Invalid YAML in {config_path}: {exc}") from exc

    if not isinstance(document, dict):
        raise ConfigurationError("The YAML root must be a mapping/object.")

    config = _expand_environment(deepcopy(document))
    schema_version = config.get("schema_version")
    if schema_version != SUPPORTED_SCHEMA_VERSION:
        raise ConfigurationError(
            f"Unsupported schema_version={schema_version!r}; "
            f"expected {SUPPORTED_SCHEMA_VERSION}."
        )

    _resolve_known_paths(config, config_path)
    config["_meta"] = {
        "config_path": str(config_path),
        "config_directory": str(config_path.parent),
    }
    return config


def unresolved_environment_tokens(config: Any, location: str = "config") -> list[str]:
    """Return paths to strings that still contain an unexpanded ${VAR}."""

    found: list[str] = []
    if isinstance(config, str) and "${" in config:
        found.append(location)
    elif isinstance(config, list):
        for index, value in enumerate(config):
            found.extend(unresolved_environment_tokens(value, f"{location}[{index}]"))
    elif isinstance(config, dict):
        for key, value in config.items():
            found.extend(unresolved_environment_tokens(value, f"{location}.{key}"))
    return found
