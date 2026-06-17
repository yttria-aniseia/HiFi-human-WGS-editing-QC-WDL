#!/usr/bin/env python3
"""
Validate a CRISPR expected-edit JSON against the authoritative schema plus a few
semantic rules that the schema alone can't express.

Usage:
    python3 validate_edit_json.py <edit.json> [<edit2.json> ...]

Exit status is non-zero if any file fails. Schema validation uses `jsonschema` if
available; otherwise it falls back to a built-in structural check covering the same
required fields, so the script is useful even in a bare environment.
"""

import json
import sys
from pathlib import Path

# Authoritative schema lives in the workflow tree; locate it relative to repo root.
REPO_ROOT = Path(__file__).resolve().parents[3]
SCHEMA_PATH = REPO_ROOT / "workflows" / "edit-qc" / "crispr_edit.schema.json"

DNA = set("ACGTNacgtn")
SUPPORTED_TYPES = {"INS"}          # only large INS is exercised end-to-end
SCHEMA_TYPES = {"SNV", "INS", "DEL"}  # what the schema permits


def _structural_check(doc, errors):
    """Minimal required-field check used when jsonschema is unavailable."""
    if not isinstance(doc, dict):
        errors.append("top level must be an object")
        return
    if doc.get("ref") != "hg38":
        errors.append('"ref" must be "hg38"')
    edits = doc.get("edits")
    if not isinstance(edits, list) or not edits:
        errors.append('"edits" must be a non-empty array')
        return
    for i, e in enumerate(edits):
        p = f"edits[{i}]"
        for k in ("id", "target", "edit"):
            if k not in e:
                errors.append(f"{p}: missing required key '{k}'")
        t = e.get("target", {})
        for k in ("chr", "start", "end"):
            if k not in t:
                errors.append(f"{p}.target: missing required key '{k}'")
        ed = e.get("edit", {})
        for k in ("type", "strand", "guide", "left_ha", "payload", "right_ha"):
            if k not in ed:
                errors.append(f"{p}.edit: missing required key '{k}'")


def _schema_check(doc, errors):
    try:
        import jsonschema  # type: ignore
    except ImportError:
        _structural_check(doc, errors)
        return "structural (jsonschema not installed)"
    schema = json.loads(SCHEMA_PATH.read_text())
    validator = jsonschema.Draft202012Validator(schema)
    for err in sorted(validator.iter_errors(doc), key=lambda e: e.path):
        loc = "/".join(str(p) for p in err.path) or "<root>"
        errors.append(f"{loc}: {err.message}")
    return "jsonschema"


def _semantic_check(doc, errors, warnings):
    for i, e in enumerate(doc.get("edits", []) if isinstance(doc, dict) else []):
        p = f"edits[{i}]"
        ed = e.get("edit", {}) if isinstance(e, dict) else {}
        t = e.get("target", {}) if isinstance(e, dict) else {}

        etype = ed.get("type")
        if etype in SCHEMA_TYPES and etype not in SUPPORTED_TYPES:
            warnings.append(
                f"{p}.edit.type={etype!r} is schema-valid but NOT exercised "
                f"end-to-end; only large INS is currently supported.")

        for seq_key in ("guide", "left_ha", "payload", "right_ha",
                        "donor_context_before", "donor_context_after"):
            seq = ed.get(seq_key)
            if isinstance(seq, str) and seq and not set(seq) <= DNA:
                bad = sorted(set(seq) - DNA)
                errors.append(f"{p}.edit.{seq_key}: non-DNA characters {bad}")

        start, end = t.get("start"), t.get("end")
        if isinstance(start, (int, float)) and isinstance(end, (int, float)):
            if start < 1:
                errors.append(f"{p}.target.start must be >= 1 (1-based)")
            if end < start:
                errors.append(f"{p}.target.end ({end}) < start ({start})")
        chrom = t.get("chr")
        if isinstance(chrom, str) and not chrom.startswith("chr"):
            warnings.append(f"{p}.target.chr={chrom!r} lacks 'chr' prefix "
                            f"(reference is hg38/UCSC-style).")


def validate(path):
    errors, warnings = [], []
    try:
        doc = json.loads(Path(path).read_text())
    except (OSError, json.JSONDecodeError) as e:
        return [f"could not read/parse JSON: {e}"], []
    mode = _schema_check(doc, errors)
    _semantic_check(doc, errors, warnings)
    return errors, warnings, mode


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 2
    if not SCHEMA_PATH.exists():
        print(f"WARNING: schema not found at {SCHEMA_PATH}", file=sys.stderr)
    any_fail = False
    for path in argv[1:]:
        errors, warnings, mode = validate(path)
        print(f"== {path}  [{mode}] ==")
        for w in warnings:
            print(f"  WARN: {w}")
        if errors:
            any_fail = True
            for e in errors:
                print(f"  ERROR: {e}")
            print(f"  -> INVALID ({len(errors)} error(s))")
        else:
            print("  -> OK")
    return 1 if any_fail else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
