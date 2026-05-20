---
name: prepare-edit-inputs
description: Prepare and validate inputs for a HiFi WGS CRISPR edit-QC family run — build the input config JSON, generate expected-edit JSONs from Benchling GenBank files, fix ref-map placeholders, handle dbNSFP licensing, and reuse precomputed assemblies. Use when a user is onboarding, setting up a new run, or preparing input data.
---

# Prepare edit-QC inputs

Goal: produce a valid `input_config.json` (+ per-sample expected-edit JSONs) that
`scripts/launch.sh` can stage and run. Work through the steps below; do not skip validation.

## 0. Preflight: is the reference data present?

Before anything else, confirm `scripts/setup.sh` has been run (see AGENTS.md guardrail 1):
- Check `hifi-wdl-resources-v3.1.0/` exists and the paths inside
  `GRCh38.ref_map.v3p1p0.template.tsv` resolve to real files.
- **`<prefix>` placeholders:** unmodified templates contain `<prefix>/...` paths.
  `setup.sh` strips these (`sed "s/<prefix>\///g"`). If you still see `<prefix>` in a map
  file, the template was never populated — tell the user to run `scripts/setup.sh`. Do not
  hand-edit placeholders into guessed absolute paths.

## 1. dbNSFP licensing (do this once, flag every time)

`setup.sh` auto-downloads **dbNSFP v4.9a** as a *fallback only*. v4.x omits columns that are
license-gated for commercial use. **Academic users should register at
https://www.dbnsfp.org/download and obtain dbNSFP v5.3+**, then:
- place the indexed file alongside the other references,
- update `DBNSFP_VERSION` in `setup.sh` (or the download step) and
- update the `dbnsfp` entry in `GRCh38.somatic_map.v3p1p0.template.tsv` to point at it.

The agent cannot download a registered file — surface this to the user and confirm which
dbNSFP version their somatic map points at.

## 2. Expected-edit JSON from a Benchling GenBank file

Use `genbank_to_crispr_json.py`. The GenBank record must have features labeled:
- a **left homology arm** (label matching `5' Homology` / `Left HA` / `Left Homology Arm`),
- a **right homology arm** (`3' Homology` / `Right HA` / `Right Homology Arm`),
- a feature labeled exactly **`payload`**,
- payload sub-components as `misc_feature`s *inside* the payload (optional; become `payload_components`).
- Sequence flanking the HDR cassette (e.g. plasmid backbone) is captured as
  `donor_context_before`/`donor_context_after` automatically.

```bash
python3 genbank_to_crispr_json.py \
  -i plasmid.gb -o my_inputs/<edit>.json \
  -e <edit_id> -t <ENSEMBL_transcript_id> -g <gRNA_seq> \
  --type INS [--strand +|-]
```

Target coordinates come from the Ensembl REST API for the transcript ID (needs network).
**Only `--type INS` is exercised end-to-end** — warn if the user wants SNV/DEL.

See `my_inputs/genbank/` for example `.gb` files and `example_inputs/example_expected_edit.json`
for a known-good output.

## 3. Validate the expected-edit JSON

Always run the bundled validator against the authoritative schema before wiring it into a run:

```bash
python3 .claude/skills/prepare-edit-inputs/validate_edit_json.py my_inputs/<edit>.json
```

It checks `workflows/edit-qc/crispr_edit.schema.json` plus extra semantic rules (large-INS
only, DNA alphabet, sane coordinates). Fix every error it reports.

## 4. Build the input config JSON

Model after `my_inputs/all_KOLF21J.json` (and `README.md`). Shape:

```json
{
  "humanwgs_family.family": {
    "family_id": "my_experiment",
    "samples": [
      { "sample_id": "parent_wt", "hifi_reads": ["/abs/parent.bam"],
        "sex": "FEMALE", "affected": false },
      { "sample_id": "clone1", "hifi_reads": ["/abs/clone1.bam"],
        "sex": "FEMALE", "affected": true, "mother_id": "parent_wt",
        "expected_edits": "/abs/my_inputs/<edit>.json" }
    ]
  },
  "humanwgs_family.ref_map_file": "/abs/GRCh38.ref_map.v3p1p0.template.tsv",
  "humanwgs_family.tertiary_map_file": "/abs/GRCh38.tertiary_map.v3p1p0.template.tsv",
  "humanwgs_family.somatic_map_file": "/abs/GRCh38.somatic_map.v3p1p0.template.tsv"
}
```

`hifi_reads` may list multiple BAMs per sample; `process_input_config.py` merges + strips them.
You give raw source paths here — staging happens in `launch.sh`, not in this file.

## 5. Reuse precomputed assemblies (big time saver — always ask)

`hifiasm` de novo assembly is the longest step. When a parent line is reused or a family is
rerun, **prompt the user for archived per-sample assembly FASTAs** and add to the sample:

```json
"precomputed_asm_hap1": "/abs/<sample>.hap1.fa.gz",
"precomputed_asm_hap2": "/abs/<sample>.hap2.fa.gz"
```

(Fields consumed at `workflows/family.wdl:167-168`.) Check prior run output dirs / the archive
path for these before committing to a fresh assembly. If the user has an archived data path,
ask for it explicitly.

## 6. Hand off

Once the config validates, the run is launched with `scripts/launch.sh` — see the
**`run-and-monitor`** skill.
