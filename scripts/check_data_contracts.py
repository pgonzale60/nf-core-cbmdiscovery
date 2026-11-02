#!/usr/bin/env python3
"""Authoritative data contracts and schema validation for nf-core/cbmdiscovery.

Validates:
  1. `assets/schemas/inter_module_tsv_contracts.json` and `assets/schemas/cbm_census.schema.json`
     exist, are valid JSON, and contain expected top-level properties and required keys.
  2. Cross-checks live Python module code to guarantee zero schema drift:
     - `modules/local/baitgermlinejunctions/resources/usr/bin/cbm_bait_germline_junctions.py` (lines 260-280)
       against `junction_catalog.tsv` columns.
     - `modules/local/memecoregrid/resources/usr/lib/memegrid/grid_runner.py` (lines 340-355)
       against `core_grid_summary.tsv` columns.
     - `modules/local/bipartiteprofiler/resources/usr/lib/bipartite/spacing.py` (lines 200-225)
       against `bipartite_spacing_summary.tsv` columns.
  3. If `scripts/extract_code_contracts.py` exists, executes `python3 scripts/extract_code_contracts.py --check`
     and asserts exit code 0.
"""

from __future__ import annotations

import json
from pathlib import Path
import re
import subprocess
import sys


REPO_ROOT = Path(__file__).resolve().parent.parent
TSV_CONTRACTS_PATH = REPO_ROOT / "assets" / "schemas" / "inter_module_tsv_contracts.json"
CENSUS_SCHEMA_PATH = REPO_ROOT / "assets" / "schemas" / "cbm_census.schema.json"


def extract_fieldnames_from_lines(
    file_path: Path, start_line: int, end_line: int
) -> list[str]:
    """Extract list of strings from fieldnames = [...] in the specified line range (1-indexed)."""
    assert file_path.is_file(), f"Source file does not exist: {file_path}"
    lines = file_path.read_text(encoding="utf-8").splitlines()
    # 1-indexed line range with some padding
    sub_lines = lines[max(0, start_line - 10) : min(len(lines), end_line + 10)]
    sub_text = "\n".join(sub_lines)

    match = re.search(r"fieldnames\s*=\s*\[(.*?)\]", sub_text, re.DOTALL)
    assert match, f"Could not find 'fieldnames = [...]' in {file_path} (lines {start_line}-{end_line})"

    items = re.findall(r"['\"]([^'\"]+)['\"]", match.group(1))
    assert items, f"No fieldname strings extracted from {file_path}"
    return items


def verify_schemas_structure() -> tuple[dict, dict]:
    """Validate existence, JSON syntax, and key requirements for schema files."""
    assert TSV_CONTRACTS_PATH.is_file(), f"Missing schema: {TSV_CONTRACTS_PATH}"
    assert CENSUS_SCHEMA_PATH.is_file(), f"Missing schema: {CENSUS_SCHEMA_PATH}"

    try:
        tsv_contracts = json.loads(TSV_CONTRACTS_PATH.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise AssertionError(f"Malformed JSON in {TSV_CONTRACTS_PATH}: {exc}") from exc

    try:
        census_schema = json.loads(CENSUS_SCHEMA_PATH.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise AssertionError(f"Malformed JSON in {CENSUS_SCHEMA_PATH}: {exc}") from exc

    # 1. Validate inter_module_tsv_contracts.json
    assert tsv_contracts.get("schema_version") == "1.0.0", (
        f"Unexpected schema_version in {TSV_CONTRACTS_PATH}: {tsv_contracts.get('schema_version')}"
    )
    assert "description" in tsv_contracts, f"Missing 'description' in {TSV_CONTRACTS_PATH}"
    contracts = tsv_contracts.get("contracts", {})
    assert isinstance(contracts, dict), f"'contracts' must be an object in {TSV_CONTRACTS_PATH}"

    expected_contracts = {
        "junction_catalog.tsv",
        "core_grid_summary.tsv",
        "bipartite_spacing_summary.tsv",
    }
    missing_contracts = expected_contracts - set(contracts.keys())
    assert not missing_contracts, f"Missing contracts in {TSV_CONTRACTS_PATH}: {missing_contracts}"

    # Check junction_catalog.tsv contract specification
    jc = contracts["junction_catalog.tsv"]
    assert jc.get("producer_module") == "BAITGERMLINEJUNCTIONS"
    assert jc.get("status_column") == "classification"
    assert set(jc.get("allowed_statuses", [])) == {
        "PDE_BREAK_SITE",
        "GERMLINE_NATIVE_TELOMERE",
        "SHARED_NATIVE_TELOMERE",
        "TANDEM_SATELLITE_ARRAY",
        "UNBROKEN_GERMLINE_LOCUS",
        "LOW_COHERENCE_JUNCTION",
        "CANDIDATE_TELOMERE",
    }
    assert isinstance(jc.get("columns"), list)

    # Check core_grid_summary.tsv contract specification
    cg = contracts["core_grid_summary.tsv"]
    assert cg.get("producer_module") == "MEMECOREGRID"
    assert cg.get("allowed_models") == ["zoops", "anr", "oops"]
    assert isinstance(cg.get("columns"), list)

    # Check bipartite_spacing_summary.tsv contract specification
    bp = contracts["bipartite_spacing_summary.tsv"]
    assert bp.get("producer_module") == "BIPARTITEPROFILER"
    assert isinstance(bp.get("columns"), list)

    # 2. Validate cbm_census.schema.json
    assert census_schema.get("$schema") == "https://json-schema.org/draft/2020-12/schema", (
        f"Unexpected $schema in {CENSUS_SCHEMA_PATH}"
    )
    assert census_schema.get("title") == "CBM Coverage and Chromosome Census Schema"
    assert census_schema.get("type") == "object"

    required_census_keys = [
        "median_c_soma",
        "median_c_germ",
        "germline_fraction",
        "n_breaks",
        "n_germ_telo",
        "estimated_germline_chromosomes",
        "n_shared_telo",
        "n_satellite_excluded",
    ]
    actual_required = census_schema.get("required", [])
    assert actual_required == required_census_keys, (
        f"Mismatch in census schema 'required' fields: expected {required_census_keys}, got {actual_required}"
    )

    props = census_schema.get("properties", {})
    for key in required_census_keys:
        assert key in props, f"Missing property definition '{key}' in {CENSUS_SCHEMA_PATH}"
        assert "type" in props[key], f"Missing 'type' in property '{key}'"
        assert "description" in props[key], f"Missing 'description' in property '{key}'"

    print("✓ Schema files exist, are valid JSON, and meet structural specifications.")
    return tsv_contracts, census_schema


def verify_live_code_contracts(tsv_contracts: dict) -> None:
    """Verify live Python module code against contract column specifications."""
    contracts = tsv_contracts["contracts"]

    # 1. BAITGERMLINEJUNCTIONS (junction_catalog.tsv)
    jc_path = (
        REPO_ROOT
        / "modules"
        / "local"
        / "baitgermlinejunctions"
        / "resources"
        / "usr"
        / "bin"
        / "cbm_bait_germline_junctions.py"
    )
    jc_cols = extract_fieldnames_from_lines(jc_path, start_line=260, end_line=280)
    expected_jc = contracts["junction_catalog.tsv"]["columns"]
    assert jc_cols == expected_jc, (
        f"Schema drift in {jc_path}: code fieldnames {jc_cols} != contract columns {expected_jc}"
    )
    print("✓ BAITGERMLINEJUNCTIONS (junction_catalog.tsv) columns match live code (lines 260-280).")

    # 2. MEMECOREGRID (core_grid_summary.tsv)
    cg_path = (
        REPO_ROOT
        / "modules"
        / "local"
        / "memecoregrid"
        / "resources"
        / "usr"
        / "lib"
        / "memegrid"
        / "grid_runner.py"
    )
    cg_cols = extract_fieldnames_from_lines(cg_path, start_line=340, end_line=355)
    expected_cg = contracts["core_grid_summary.tsv"]["columns"]
    assert cg_cols == expected_cg, (
        f"Schema drift in {cg_path}: code fieldnames {cg_cols} != contract columns {expected_cg}"
    )
    print("✓ MEMECOREGRID (core_grid_summary.tsv) columns match live code (lines 340-355).")

    # 3. BIPARTITEPROFILER (bipartite_spacing_summary.tsv)
    bp_path = (
        REPO_ROOT
        / "modules"
        / "local"
        / "bipartiteprofiler"
        / "resources"
        / "usr"
        / "lib"
        / "bipartite"
        / "spacing.py"
    )
    bp_cols = extract_fieldnames_from_lines(bp_path, start_line=200, end_line=225)
    expected_bp = contracts["bipartite_spacing_summary.tsv"]["columns"]
    assert bp_cols == expected_bp, (
        f"Schema drift in {bp_path}: code fieldnames {bp_cols} != contract columns {expected_bp}"
    )
    print("✓ BIPARTITEPROFILER (bipartite_spacing_summary.tsv) columns match live code (lines 200-225).")


def verify_code_contracts_extractor() -> None:
    """Run extract_code_contracts.py --check if present."""
    extractor_script = REPO_ROOT / "scripts" / "extract_code_contracts.py"
    if extractor_script.is_file():
        res = subprocess.run(
            [sys.executable, str(extractor_script), "--check"],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
        )
        assert res.returncode == 0, (
            f"extract_code_contracts.py --check failed with returncode {res.returncode}:\n"
            f"STDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}"
        )
        print("✓ scripts/extract_code_contracts.py --check passed with exit code 0.")
    else:
        print("ℹ scripts/extract_code_contracts.py not present; skipping check.")


def main() -> int:
    tsv_contracts, _ = verify_schemas_structure()
    verify_live_code_contracts(tsv_contracts)
    verify_code_contracts_extractor()
    print("check_data_contracts: PASS (all TSV/JSON schemas and live code contracts verified)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
