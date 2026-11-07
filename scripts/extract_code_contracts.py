#!/usr/bin/env python3
"""Deterministically extract Nextflow topology, inter-module data contracts, and Python AST symbols into docs/CODE_CONTRACTS.md.

Provides a zero-drift architectural map for developers and AI agents:
  1. Nextflow Workflow, Subworkflow & Module Topology (`include`, `take:`, `emit:`, `input:`, `output:`)
  2. Inter-Module Data Contracts (rendered from `assets/schemas/cbm_census.schema.json` and
     `assets/schemas/inter_module_tsv_contracts.json`)
  3. Python Module AST Symbol Table (exact signatures, `#Lstart-Lend` line ranges, and 1-line docstrings
     for all `.py` files under `modules/local/*/resources/usr/{bin,lib}/**/*.py`)

Uses Python stdlib only (`ast`, `json`, `pathlib`, `re`, `argparse`, `subprocess`).

Usage:
    python3 scripts/extract_code_contracts.py --all          # Regenerate docs/CODE_CONTRACTS.md
    python3 scripts/extract_code_contracts.py --pre-commit   # Regenerate & git-add if staged files affect contracts
    python3 scripts/extract_code_contracts.py --check        # Verify docs/CODE_CONTRACTS.md is up to date
"""

import argparse
import ast
import json
import pathlib
import re
import subprocess
import sys


REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
DOCS_OUT = REPO_ROOT / "docs" / "CODE_CONTRACTS.md"
CENSUS_SCHEMA_PATH = REPO_ROOT / "assets" / "schemas" / "cbm_census.schema.json"
TSV_CONTRACTS_PATH = REPO_ROOT / "assets" / "schemas" / "inter_module_tsv_contracts.json"

INCLUDE_RE = re.compile(
    r"^[ \t]*include\s*\{\s*([^}]+?)\s*\}\s*from\s*['\"]([^'\"]+)['\"]",
    re.MULTILINE,
)
WORKFLOW_HEADER_RE = re.compile(r"^[ \t]*workflow(?:\s+([A-Za-z0-9_]+))?\s*\{", re.MULTILINE)
PROCESS_HEADER_RE = re.compile(r"^[ \t]*process\s+([A-Za-z0-9_]+)\s*\{", re.MULTILINE)


def _strip_nf_comments(line: str) -> str:
    """Strip trailing // comments outside quotes and collapse whitespace."""
    in_single = False
    in_double = False
    for idx, ch in enumerate(line):
        if ch == "'" and not in_double:
            in_single = not in_single
        elif ch == '"' and not in_single:
            in_double = not in_double
        elif ch == "/" and not in_single and not in_double:
            if idx + 1 < len(line) and line[idx + 1] == "/":
                return line[:idx].strip()
    return line.strip()


def _extract_nf_block_statements(lines: list[str], start_idx: int, stop_directives: set[str]) -> list[str]:
    """Extract logical declarations from a Nextflow directive section (`take:`, `emit:`, `input:`, `output:`)."""
    items: list[str] = []
    cur_parts: list[str] = []
    paren_depth = 0

    for raw in lines[start_idx:]:
        cleaned = _strip_nf_comments(raw)
        if not cleaned:
            continue
        if paren_depth == 0 and (cleaned in stop_directives or cleaned == "}"):
            break

        open_parens = cleaned.count("(") + cleaned.count("[")
        close_parens = cleaned.count(")") + cleaned.count("]")
        paren_depth += open_parens - close_parens
        cur_parts.append(cleaned)

        if paren_depth <= 0:
            paren_depth = 0
            stmt = " ".join(cur_parts)
            stmt = re.sub(r"\s+", " ", stmt).strip()
            if stmt:
                items.append(stmt)
            cur_parts = []

    if cur_parts:
        stmt = re.sub(r"\s+", " ", " ".join(cur_parts)).strip()
        if stmt:
            items.append(stmt)
    return items


def parse_nextflow_file(nf_path: pathlib.Path) -> dict:
    """Parse a Nextflow `.nf` file for includes, workflow blocks, and process blocks."""
    text = nf_path.read_text(encoding="utf-8")
    lines = text.splitlines()

    includes = []
    for match in INCLUDE_RE.finditer(text):
        raw_symbols = match.group(1)
        source = match.group(2)
        symbols = [" ".join(part.strip().split()) for part in raw_symbols.split(";") if part.strip()]
        includes.append({"symbols": symbols, "from": source})

    workflows = []
    for match in WORKFLOW_HEADER_RE.finditer(text):
        wf_name = match.group(1) or "<entry>"
        start_line = text[: match.start()].count("\n")
        brace_depth = 0
        wf_lines = []
        started = False
        for line in lines[start_line:]:
            wf_lines.append(line)
            brace_depth += line.count("{") - line.count("}")
            if "{" in line:
                started = True
            if started and brace_depth <= 0:
                break

        take_items = []
        emit_items = []
        stop_set = {"take:", "main:", "emit:", "publish:"}
        for idx, line in enumerate(wf_lines):
            stripped = _strip_nf_comments(line)
            if stripped == "take:":
                take_items = _extract_nf_block_statements(wf_lines, idx + 1, stop_set)
            elif stripped == "emit:":
                emit_items = [
                    s.split("=")[0].strip() if "=" in s else s
                    for s in _extract_nf_block_statements(wf_lines, idx + 1, stop_set)
                ]
        workflows.append({"name": wf_name, "take": take_items, "emit": emit_items})

    processes = []
    for match in PROCESS_HEADER_RE.finditer(text):
        proc_name = match.group(1)
        start_line = text[: match.start()].count("\n")
        brace_depth = 0
        proc_lines = []
        started = False
        for line in lines[start_line:]:
            proc_lines.append(line)
            brace_depth += line.count("{") - line.count("}")
            if "{" in line:
                started = True
            if started and brace_depth <= 0:
                break

        input_items = []
        output_items = []
        stop_set = {"input:", "output:", "when:", "script:", "shell:", "exec:", "stub:"}
        for idx, line in enumerate(proc_lines):
            stripped = _strip_nf_comments(line)
            if stripped == "input:":
                input_items = _extract_nf_block_statements(proc_lines, idx + 1, stop_set)
            elif stripped == "output:":
                output_items = _extract_nf_block_statements(proc_lines, idx + 1, stop_set)
        processes.append({"name": proc_name, "input": input_items, "output": output_items})

    return {
        "rel_path": str(nf_path.relative_to(REPO_ROOT)),
        "includes": includes,
        "workflows": workflows,
        "processes": processes,
    }


def _format_arg(arg: ast.arg, default_node: ast.AST | None = None) -> str:
    text = arg.arg
    if arg.annotation is not None:
        text += f": {ast.unparse(arg.annotation)}"
    if default_node is not None:
        text += f"={ast.unparse(default_node)}"
    return text


def format_function_signature(node: ast.FunctionDef | ast.AsyncFunctionDef) -> str:
    """Format an AST FunctionDef/AsyncFunctionDef into a deterministic signature string."""
    args = node.args
    parts: list[str] = []

    pos_args = list(args.posonlyargs) + list(args.args)
    num_defaults = len(args.defaults)
    default_offset = len(pos_args) - num_defaults

    for idx, arg in enumerate(pos_args):
        default_node = args.defaults[idx - default_offset] if idx >= default_offset else None
        parts.append(_format_arg(arg, default_node))
        if args.posonlyargs and idx == len(args.posonlyargs) - 1:
            parts.append("/")

    if args.vararg is not None:
        parts.append("*" + _format_arg(args.vararg))
    elif args.kwonlyargs:
        parts.append("*")

    for kw_arg, kw_default in zip(args.kwonlyargs, args.kw_defaults):
        parts.append(_format_arg(kw_arg, kw_default))

    if args.kwarg is not None:
        parts.append("**" + _format_arg(args.kwarg))

    ret = f" -> {ast.unparse(node.returns)}" if node.returns is not None else ""
    prefix = "async def " if isinstance(node, ast.AsyncFunctionDef) else "def "
    return f"{prefix}{node.name}({', '.join(parts)}){ret}"


def _one_line_doc(node: ast.AST) -> str:
    doc = ast.get_docstring(node)
    if not doc:
        return "—"
    first_line = doc.strip().splitlines()[0].strip()
    return first_line.replace("|", r"\|")


def extract_python_symbols(py_path: pathlib.Path) -> list[dict]:
    """Parse a Python file via `ast.parse` and return top-level classes, methods, and functions."""
    source = py_path.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(py_path))
    symbols: list[dict] = []

    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            symbols.append(
                {
                    "kind": "function",
                    "name": node.name,
                    "signature": format_function_signature(node),
                    "lineno": node.lineno,
                    "end_lineno": getattr(node, "end_lineno", node.lineno),
                    "doc": _one_line_doc(node),
                }
            )
        elif isinstance(node, ast.ClassDef):
            bases = [ast.unparse(b) for b in node.bases]
            base_str = f"({', '.join(bases)})" if bases else ""
            symbols.append(
                {
                    "kind": "class",
                    "name": node.name,
                    "signature": f"class {node.name}{base_str}",
                    "lineno": node.lineno,
                    "end_lineno": getattr(node, "end_lineno", node.lineno),
                    "doc": _one_line_doc(node),
                }
            )
            for item in node.body:
                if isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    sig = format_function_signature(item)
                    sig = re.sub(r"^(async def |def )", rf"\1{node.name}.", sig, count=1)
                    symbols.append(
                        {
                            "kind": "method",
                            "name": f"{node.name}.{item.name}",
                            "signature": sig,
                            "lineno": item.lineno,
                            "end_lineno": getattr(item, "end_lineno", item.lineno),
                            "doc": _one_line_doc(item),
                        }
                    )
    return symbols


def discover_python_files() -> list[pathlib.Path]:
    """Return deterministically sorted `.py` files under `modules/local/*/resources/usr/{bin,lib}/`."""
    modules_dir = REPO_ROOT / "modules" / "local"
    py_files: set[pathlib.Path] = set()
    for sub in ("bin", "lib"):
        for p in modules_dir.glob(f"*/resources/usr/{sub}/**/*.py"):
            if p.is_file():
                py_files.add(p)
        for p in modules_dir.glob(f"*/resources/usr/{sub}/*.py"):
            if p.is_file():
                py_files.add(p)
    return sorted(py_files, key=lambda p: str(p.relative_to(REPO_ROOT)))


def discover_nextflow_files() -> tuple[list[pathlib.Path], list[pathlib.Path], list[pathlib.Path]]:
    """Return sorted (entry_and_workflows, subworkflows, modules) Nextflow files."""
    entry_wf = [REPO_ROOT / "main.nf"]
    entry_wf.extend(sorted((REPO_ROOT / "workflows").glob("*.nf")))
    subworkflows = sorted((REPO_ROOT / "subworkflows" / "local").glob("*/main.nf"))
    mod_set = set((REPO_ROOT / "modules" / "local").glob("*/main.nf")) | set(
        (REPO_ROOT / "modules" / "local").glob("*/*/main.nf")
    )
    modules = sorted(mod_set, key=lambda p: str(p.relative_to(REPO_ROOT)))
    return entry_wf, subworkflows, modules


def render_code_contracts_markdown() -> str:
    """Generate the complete deterministic Markdown content for `docs/CODE_CONTRACTS.md`."""
    lines: list[str] = []
    lines.append("# Code Contracts, Nextflow Topology & Python AST Symbol Table")
    lines.append("")
    lines.append(
        "> **AUTO-GENERATED FILE — DO NOT EDIT MANUALLY.**  \n"
        "> Deterministically generated by `python3 scripts/extract_code_contracts.py --all` "
        "and verified in `scripts/pre-commit`."
    )
    lines.append("")
    lines.append("## Table of Contents")
    lines.append("1. [Nextflow Workflow, Subworkflow & Module Topology](#1-nextflow-workflow-subworkflow--module-topology)")
    lines.append("2. [Inter-Module Data Contracts (JSON & TSV Schemas)](#2-inter-module-data-contracts-json--tsv-schemas)")
    lines.append("3. [Python Module AST Symbol Table](#3-python-module-ast-symbol-table)")
    lines.append("")

    # -------------------------------------------------------------------------
    # SECTION 1: Nextflow Topology
    # -------------------------------------------------------------------------
    lines.append("## 1. Nextflow Workflow, Subworkflow & Module Topology")
    lines.append("")
    entry_wf, subworkflows, modules = discover_nextflow_files()

    lines.append("### 1.1 Entrypoint & Top-Level Workflows (`main.nf`, `workflows/*.nf`)")
    lines.append("")
    for nf_path in entry_wf:
        info = parse_nextflow_file(nf_path)
        lines.append(f"#### `{info['rel_path']}`")
        if info["includes"]:
            lines.append("- **Includes:**")
            for inc in info["includes"]:
                lines.append(f"  - `{', '.join(inc['symbols'])}` from `{inc['from']}`")
        for wf in info["workflows"]:
            lines.append(f"- **Workflow `{wf['name']}`:**")
            lines.append(f"  - `take:` {', '.join(f'`{t}`' for t in wf['take']) if wf['take'] else '*(none)*'}")
            lines.append(f"  - `emit:` {', '.join(f'`{e}`' for e in wf['emit']) if wf['emit'] else '*(none)*'}")
        lines.append("")

    lines.append("### 1.2 Local Subworkflows (`subworkflows/local/*/main.nf`)")
    lines.append("")
    if not subworkflows:
        lines.append("*(No local subworkflows defined.)*")
        lines.append("")
    else:
        for nf_path in subworkflows:
            info = parse_nextflow_file(nf_path)
            lines.append(f"#### `{info['rel_path']}`")
            if info["includes"]:
                lines.append("- **Includes:**")
                for inc in info["includes"]:
                    lines.append(f"  - `{', '.join(inc['symbols'])}` from `{inc['from']}`")
            for wf in info["workflows"]:
                lines.append(f"- **Workflow `{wf['name']}`:**")
                lines.append(f"  - `take:` {', '.join(f'`{t}`' for t in wf['take']) if wf['take'] else '*(none)*'}")
                lines.append(f"  - `emit:` {', '.join(f'`{e}`' for e in wf['emit']) if wf['emit'] else '*(none)*'}")
            lines.append("")

    lines.append("### 1.3 Local Modules (`modules/local/**/main.nf`)")
    lines.append("")
    lines.append("| Process | File | Inputs (`input:`) | Outputs (`output:`) |")
    lines.append("| :--- | :--- | :--- | :--- |")
    escaped_pipe = r"\|"
    for nf_path in modules:
        info = parse_nextflow_file(nf_path)
        for proc in info["processes"]:
            inp_str = "<br>".join(f"`{i.replace('|', escaped_pipe)}`" for i in proc["input"]) if proc["input"] else "—"
            out_str = "<br>".join(f"`{o.replace('|', escaped_pipe)}`" for o in proc["output"]) if proc["output"] else "—"
            lines.append(f"| `{proc['name']}` | `{info['rel_path']}` | {inp_str} | {out_str} |")
    lines.append("")

    # -------------------------------------------------------------------------
    # SECTION 2: Inter-Module Data Contracts
    # -------------------------------------------------------------------------
    lines.append("## 2. Inter-Module Data Contracts (JSON & TSV Schemas)")
    lines.append("")

    if CENSUS_SCHEMA_PATH.is_file():
        lines.append("### 2.1 CBM Census Schema (`assets/schemas/cbm_census.schema.json`)")
        lines.append("")
        try:
            census_schema = json.loads(CENSUS_SCHEMA_PATH.read_text(encoding="utf-8"))
            lines.append(
                f"- **Schema URI:** `{census_schema.get('$schema', 'https://json-schema.org/draft/2020-12/schema')}`"
            )
            if "title" in census_schema:
                lines.append(f"- **Title:** {census_schema['title']}")
            if "description" in census_schema:
                lines.append(f"- **Description:** {census_schema['description']}")
            lines.append(
                "- **Producer:** `BAITGERMLINEJUNCTIONS` (`cbm_bait_germline_junctions.py` -> `coverage_and_chromosome_census.json`)"
            )
            req = census_schema.get("required", [])
            if req:
                lines.append(f"- **Required Top-Level Keys:** {', '.join(f'`{k}`' for k in req)}")
            props = census_schema.get("properties", {})
            if props:
                lines.append("")
                lines.append("| Property | Type | Description |")
                lines.append("| :--- | :--- | :--- |")
                for prop_name, prop_spec in sorted(props.items()):
                    ptype = prop_spec.get("type", "—")
                    pdesc = prop_spec.get("description", "—").replace("|", r"\|")
                    lines.append(f"| `{prop_name}` | `{ptype}` | {pdesc} |")
        except Exception as exc:
            lines.append(f"*(Error parsing schema: {exc})*")
        lines.append("")
    else:
        lines.append("### 2.1 CBM Census Schema (`assets/schemas/cbm_census.schema.json`)")
        lines.append("")
        lines.append("*(Schema file `assets/schemas/cbm_census.schema.json` not yet present.)*")
        lines.append("")

    if TSV_CONTRACTS_PATH.is_file():
        lines.append("### 2.2 Inter-Module TSV Contracts (`assets/schemas/inter_module_tsv_contracts.json`)")
        lines.append("")
        try:
            tsv_contracts = json.loads(TSV_CONTRACTS_PATH.read_text(encoding="utf-8"))
            contracts = tsv_contracts.get("contracts", {})
            if contracts:
                lines.append(
                    "| Contract TSV | Producer Module | Consumer Module(s) | Ordered Columns | Status Enums / Details |"
                )
                lines.append("| :--- | :--- | :--- | :--- | :--- |")
                for tsv_name, spec in sorted(contracts.items()):
                    prod = f"`{spec.get('producer_module', '—')}`"
                    cons_val = spec.get("consumer_modules", spec.get("consumer_module", []))
                    if isinstance(cons_val, list):
                        cons = ", ".join(f"`{c}`" for c in cons_val) if cons_val else "—"
                    else:
                        cons = f"`{cons_val}`"
                    cols_val = spec.get("columns", [])
                    if isinstance(cols_val, list):
                        cols_list = []
                        for c in cols_val:
                            if isinstance(c, dict):
                                cols_list.append(f"`{c.get('name', '')}`")
                            else:
                                cols_list.append(f"`{c}`")
                        cols = ", ".join(cols_list) if cols_list else "—"
                    else:
                        cols = f"`{cols_val}`"

                    statuses = (
                        spec.get("status_enums")
                        or spec.get("allowed_statuses")
                        or spec.get("spliceable_statuses")
                        or spec.get("statuses")
                        or []
                    )
                    if isinstance(statuses, list) and statuses:
                        stat_str = ", ".join(f"`{s}`" for s in statuses)
                    elif isinstance(statuses, dict):
                        stat_str = ", ".join(f"`{k}`" for k in statuses)
                    else:
                        stat_str = "—"
                    lines.append(f"| `{tsv_name}` | {prod} | {cons} | {cols} | {stat_str} |")
            else:
                lines.append("*(No contracts declared in JSON.)*")
        except Exception as exc:
            lines.append(f"*(Error parsing TSV contracts: {exc})*")
        lines.append("")
    else:
        lines.append("### 2.2 Inter-Module TSV Contracts (`assets/schemas/inter_module_tsv_contracts.json`)")
        lines.append("")
        lines.append("*(Contracts file `assets/schemas/inter_module_tsv_contracts.json` not yet present.)*")
        lines.append("")

    # -------------------------------------------------------------------------
    # SECTION 3: Python Module AST Symbol Table
    # -------------------------------------------------------------------------
    lines.append("## 3. Python Module AST Symbol Table")
    lines.append("")
    for py_path in discover_python_files():
        rel_path = str(py_path.relative_to(REPO_ROOT))
        symbols = extract_python_symbols(py_path)
        lines.append(f"### `{rel_path}`")
        if not symbols:
            lines.append("*(No top-level functions or classes defined.)*")
            lines.append("")
            continue
        lines.append("")
        lines.append("| Symbol | Kind | Lines | Signature | Summary |")
        lines.append("| :--- | :--- | :--- | :--- | :--- |")
        for sym in symbols:
            line_range = f"`L{sym['lineno']}-L{sym['end_lineno']}`"
            sig_escaped = sym["signature"].replace("|", r"\|")
            lines.append(
                f"| `{sym['name']}` | `{sym['kind']}` | {line_range} | `{sig_escaped}` | {sym['doc']} |"
            )
        lines.append("")

    return "\n".join(lines).rstrip() + "\n"


def _has_relevant_staged_changes() -> bool:
    """Return True if staged git changes touch `.nf`, `.py`, or `assets/schemas/*.json`."""
    res = subprocess.run(
        ["git", "diff", "--cached", "--name-only", "--diff-filter=ACMRD"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=True,
    )
    for line in res.stdout.splitlines():
        if (
            line.endswith(".py")
            or line.endswith(".nf")
            or line.startswith("assets/schemas/")
            or line == "docs/CODE_CONTRACTS.md"
        ):
            return True
    return False


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--all", action="store_true", help="Regenerate docs/CODE_CONTRACTS.md")
    group.add_argument(
        "--pre-commit",
        action="store_true",
        help="Regenerate docs/CODE_CONTRACTS.md and git-add if relevant files are staged",
    )
    group.add_argument("--check", action="store_true", help="Check that docs/CODE_CONTRACTS.md has zero drift")
    args = parser.parse_args()

    expected_content = render_code_contracts_markdown()
    current_content = DOCS_OUT.read_text(encoding="utf-8") if DOCS_OUT.is_file() else ""

    if args.check:
        if current_content != expected_content:
            print(
                "docs/CODE_CONTRACTS.md is out of date or missing. "
                "Run: python3 scripts/extract_code_contracts.py --all",
                file=sys.stderr,
            )
            return 1
        print("docs/CODE_CONTRACTS.md is up to date (zero contract/AST drift).")
        return 0

    if args.pre_commit:
        if not _has_relevant_staged_changes() and DOCS_OUT.is_file() and current_content == expected_content:
            return 0
        if current_content != expected_content:
            DOCS_OUT.parent.mkdir(parents=True, exist_ok=True)
            DOCS_OUT.write_text(expected_content, encoding="utf-8")
            subprocess.run(
                ["git", "add", str(DOCS_OUT.relative_to(REPO_ROOT))],
                cwd=REPO_ROOT,
                check=True,
            )
            print("[pre-commit] Auto-updated docs/CODE_CONTRACTS.md")
        return 0

    # --all
    DOCS_OUT.parent.mkdir(parents=True, exist_ok=True)
    changed = current_content != expected_content
    if changed:
        DOCS_OUT.write_text(expected_content, encoding="utf-8")
    status = "updated" if changed else "unchanged"
    print(f"Generated {DOCS_OUT.relative_to(REPO_ROOT)} ({status}, {len(expected_content.splitlines())} lines).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
