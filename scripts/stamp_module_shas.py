#!/usr/bin/env python3
"""Automatically stamp literal Script-SHA comments into Nextflow module main.nf files.

Nextflow only invalidates its `-resume` task cache when a declared input channel
changes or when the literal text of the process definition in `main.nf` changes.
Dynamic Groovy expressions inside `script:` (such as `file(...).text.sha256()`)
do not trigger cache invalidation if `main.nf` itself is untouched.

This script computes a deterministic 16-character SHA-256 digest across all Python
files in `modules/local/<module>/resources/usr/{bin,lib}/` and writes a literal
comment line inside the `script:` block of `modules/local/<module>/main.nf`:

    # Tool: <script.py> | Script-SHA: <16-char-sha256>

Usage:
    python3 scripts/stamp_module_shas.py --all          # Stamp all Python modules
    python3 scripts/stamp_module_shas.py --pre-commit   # Stamp & git-add staged modules
    python3 scripts/stamp_module_shas.py --check        # Verify all stamps are current
"""

import argparse
import hashlib
import json
import pathlib
import pprint
import re
import subprocess
import sys


REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
MODULES_DIR = REPO_ROOT / "modules" / "local"
SSOT_JSON = REPO_ROOT / "conf" / "cbm_params.json"


def sync_module_ssot_defaults(mod_dir: pathlib.Path, check_only: bool = False) -> bool:
    """Sync conf/cbm_params.json into resources/usr/lib/cbm_defaults.py for container portability."""
    if not SSOT_JSON.is_file():
        return False
    params = json.loads(SSOT_JSON.read_text(encoding="utf-8"))
    lib_dir = mod_dir / "resources" / "usr" / "lib"
    lib_dir.mkdir(parents=True, exist_ok=True)
    target_py = lib_dir / "cbm_defaults.py"
    formatted = pprint.pformat(params, width=88, sort_dicts=True)
    content = (
        '"""Auto-synced Single Source of Truth (SSOT) parameters from conf/cbm_params.json."""\n\n'
        "import json\n"
        "from pathlib import Path\n\n"
        f"_FALLBACK_SSOT = {formatted}\n\n\n"
        "def load_ssot_params() -> dict:\n"
        '    """Load canonical parameters from conf/cbm_params.json if accessible, else synced SSOT."""\n'
        "    cur = Path(__file__).resolve().parent\n"
        "    for _ in range(8):\n"
        '        cand = cur / "conf" / "cbm_params.json"\n'
        "        if cand.is_file():\n"
        '            return json.loads(cand.read_text(encoding="utf-8"))\n'
        "        if cur.parent == cur:\n"
        "            break\n"
        "        cur = cur.parent\n"
        "    return dict(_FALLBACK_SSOT)\n\n\n"
        "SSOT = load_ssot_params()\n"
    )
    existing = target_py.read_text(encoding="utf-8") if target_py.is_file() else ""
    if existing != content:
        if not check_only:
            target_py.write_text(content, encoding="utf-8")
        return True
    return False

# Stale Groovy lines previously used for dynamic fingerprinting
DYNAMIC_GROOVY_RE = re.compile(
    r"^[ \t]*def (?:lib_files|combined_src|script_fingerprint)\b.*\n",
    re.MULTILINE,
)

# Existing # Tool: ... | Script-SHA / Script-Fingerprint comment line
TOOL_COMMENT_RE = re.compile(
    r"^([ \t]*)# Tool: [^\n]*(?:Script-SHA|Script-Fingerprint): [^\n]*\n",
    re.MULTILINE,
)


def discover_python_modules() -> list[pathlib.Path]:
    """Return sorted module directories under modules/local that have usr/bin/*.py."""
    modules = []
    if not MODULES_DIR.is_dir():
        return modules
    for main_nf in sorted(MODULES_DIR.rglob("main.nf")):
        mod_dir = main_nf.parent
        bin_dir = mod_dir / "resources" / "usr" / "bin"
        if bin_dir.is_dir() and list(bin_dir.glob("*.py")):
            modules.append(mod_dir)
    return modules


def compute_module_sha(mod_dir: pathlib.Path) -> tuple[str, str]:
    """Return (primary_bin_name, sha16) for all .py files in resources/usr/{bin,lib}."""
    usr_dir = mod_dir / "resources" / "usr"
    bin_files = sorted((usr_dir / "bin").glob("*.py"))
    if not bin_files:
        raise ValueError(f"No .py files in {usr_dir / 'bin'}")

    primary_bin = bin_files[0].name
    lib_dir = usr_dir / "lib"
    lib_files = sorted(lib_dir.rglob("*.py")) if lib_dir.is_dir() else []

    hasher = hashlib.sha256()
    for py_path in bin_files + lib_files:
        hasher.update(py_path.read_bytes())
    return primary_bin, hasher.hexdigest()[:16]


def stamp_main_nf(mod_dir: pathlib.Path, check_only: bool = False) -> bool:
    """Ensure cbm_defaults.py and main.nf contain the current SSOT and Script-SHA comment.

    Returns True if cbm_defaults.py or main.nf was modified (or needs modification when check_only=True).
    """
    ssot_changed = sync_module_ssot_defaults(mod_dir, check_only=check_only)
    main_nf = mod_dir / "main.nf"
    original_text = main_nf.read_text(encoding="utf-8")
    primary_bin, sha16 = compute_module_sha(mod_dir)
    target_comment_body = f"# Tool: {primary_bin} | Script-SHA: {sha16}"

    # 1. Strip any obsolete dynamic Groovy fingerprint definitions
    updated_text = DYNAMIC_GROOVY_RE.sub("", original_text)

    # 2. Split before `stub:` so we only inspect/modify the `script:` block
    stub_idx = updated_text.find("\n    stub:")
    if stub_idx == -1:
        script_part = updated_text
        stub_part = ""
    else:
        script_part = updated_text[:stub_idx]
        stub_part = updated_text[stub_idx:]

    if TOOL_COMMENT_RE.search(script_part):
        script_part = TOOL_COMMENT_RE.sub(
            rf"\1{target_comment_body}\n",
            script_part,
            count=1,
        )
    else:
        # Insert immediately after the opening `"""` of the script block
        triple_quote_match = re.search(r'([ \t]*)"""\n', script_part)
        if not triple_quote_match:
            raise RuntimeError(f"Could not find opening '\"\"\"' in {main_nf}")
        indent = triple_quote_match.group(1)
        insert_pos = triple_quote_match.end()
        script_part = (
            script_part[:insert_pos]
            + f"{indent}{target_comment_body}\n"
            + script_part[insert_pos:]
        )

    final_text = script_part + stub_part
    changed = (final_text != original_text) or ssot_changed
    if final_text != original_text and not check_only:
        main_nf.write_text(final_text, encoding="utf-8")
    return changed


def get_staged_module_dirs() -> list[pathlib.Path]:
    """Return module dirs that have staged changes in git (or all if conf/cbm_params.json changed)."""
    res = subprocess.run(
        ["git", "diff", "--cached", "--name-only"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=True,
    )
    staged_lines = res.stdout.splitlines()
    if any("cbm_params.json" in line for line in staged_lines):
        return discover_python_modules()

    staged_mods = set()
    for line in staged_lines:
        candidate = (REPO_ROOT / line).resolve()
        try:
            candidate.relative_to(MODULES_DIR)
        except ValueError:
            continue
        cur = candidate if candidate.is_dir() else candidate.parent
        while cur != MODULES_DIR and MODULES_DIR in cur.parents:
            if (cur / "main.nf").is_file() and (cur / "resources" / "usr" / "bin").is_dir():
                if list((cur / "resources" / "usr" / "bin").glob("*.py")):
                    staged_mods.add(cur)
                break
            cur = cur.parent
    return sorted(staged_mods)


def _mod_label(mod_dir: pathlib.Path) -> str:
    return str(mod_dir.relative_to(MODULES_DIR))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--all", action="store_true", help="Stamp all sanger-tol Python modules")
    group.add_argument("--pre-commit", action="store_true", help="Stamp modules with staged changes and git-add main.nf")
    group.add_argument("--check", action="store_true", help="Check that all Python modules have current Script-SHA stamps")
    args = parser.parse_args()

    if args.pre_commit:
        target_mods = get_staged_module_dirs()
        if not target_mods:
            return 0
        stamped = []
        for mod_dir in target_mods:
            if stamp_main_nf(mod_dir, check_only=False):
                main_nf_rel = (mod_dir / "main.nf").relative_to(REPO_ROOT)
                defaults_py = mod_dir / "resources" / "usr" / "lib" / "cbm_defaults.py"
                to_add = [str(main_nf_rel)]
                if defaults_py.is_file():
                    to_add.append(str(defaults_py.relative_to(REPO_ROOT)))
                subprocess.run(["git", "add"] + to_add, cwd=REPO_ROOT, check=True)
                stamped.append(_mod_label(mod_dir))
        if stamped:
            print(f"[pre-commit] Auto-stamped Script-SHA in main.nf for: {', '.join(stamped)}")
        return 0

    all_mods = discover_python_modules()
    if args.check:
        stale = [_mod_label(mod_dir) for mod_dir in all_mods if stamp_main_nf(mod_dir, check_only=True)]
        if stale:
            print(f"Stale or missing Script-SHA in main.nf for: {', '.join(stale)}", file=sys.stderr)
            return 1
        print(f"All {len(all_mods)} Python modules have up-to-date Script-SHA stamps.")
        return 0

    # --all
    stamped = []
    for mod_dir in all_mods:
        if stamp_main_nf(mod_dir, check_only=False):
            stamped.append(_mod_label(mod_dir))
    print(f"Scanned {len(all_mods)} Python modules; updated {len(stamped)}: {', '.join(stamped) if stamped else 'none'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
