# Administrative portability helpers

The current helper is `prepare_workspace_v3.py`; its input inventory is `../provenance/ORIGIN_MANIFEST_v2.tsv`, relative to this README. From the `Code_for_Review/` directory, invoke `python portability_v2/prepare_workspace_v3.py --output /path/to/new/isolated_workspace` (the output directory must not exist). It verifies and stages 145 preserved source/configuration/contract entries without running scientific code. The older prepare_workspace.py targets the earlier inventory and is retained as historical administrative tooling.

The new output directory must not exist. Optional --runtime-root and --external-inputs arguments bind externally supplied runtimes and frozen inputs. The external-input columns are logical_path, source_path and sha256. Original-source digest checks remain active for replacement records. Do not change embedded scientific expected hashes.

The vendor directory supplies PyYAML 6.0.3 for scoped imports. The supplementary exporter reads frozen tables for presentation only. No helper constitutes a new scientific pipeline or a guarantee of full raw-data reproduction. Current scope, input limitations and legacy path bindings are described in `../../README.md`, relative to this README.
