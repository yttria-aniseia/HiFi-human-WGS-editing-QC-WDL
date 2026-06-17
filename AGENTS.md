# AGENTS.md — orientation for AI agents working in this repo

This file is loaded into context on every session. Keep it short and durable.
Procedures live in skills (see "Skills" below); this file is the map and the guardrails.

## What this pipeline is

WDL workflows for **QC of HiFi WGS CRISPR-editing experiments** (HDR knock-ins on
PacBio long reads). Primarily run via **miniwdl + apptainer/singularity on an HPC SLURM
cluster**.

- **Only the `family` workflow is supported** (`workflows/family.wdl`). Ignore singleton.
- Typical cohort: a parental (WT) line + one or more edited clones, optionally as a family.
- Edit-specific QC runs only when a sample has an `expected_edits` JSON.

## Hard guardrails (must hold even when no skill is invoked)

1. **Verify prerequisite reference data before doing anything that runs the workflow.**
   `scripts/setup.sh` downloads ~200 GB into `hifi-wdl-resources-v3.1.0/` and populates the
   map templates (`GRCh38.{ref,tertiary,somatic}_map.v3p1p0.template.tsv`). If the bundle dir
   or the paths referenced inside the map files are **missing, STOP and tell the user to run
   `scripts/setup.sh`** — never fabricate or guess reference paths.

2. **Never read pipeline logs in full.** `workflow.log` and task stderr are enormous and
   highly repetitive. Always `grep -E "ERROR|FAILED|failed"`, `tail`, or scope to one
   `call-<task>/` dir. Reading a full log wastes the entire context window.

   **Likewise, never run a broad `find`/`grep -r`.** This is an HPC tree: even a repo-scoped
   search descends into every existing pipeline work dir (`*_*/`, `workflow_run_*/`,
   `outputs/<run>/call-*/`), which hold hundreds of thousands of large intermediate files —
   such a search is slow, can hammer the shared filesystem, and floods context. Always scope
   to a known subdir (`workflows/`, `scripts/`, a single run's `out/`), prune work dirs
   (`-prune`), or use the project's structure (`outputs.json`, `out/`) instead of scanning.

3. **All workflow inputs must be staged under the run directory.** apptainer bind-mounts only
   paths under the cwd/run dir, which is why `scripts/launch.sh` → `process_input_config.py`
   copies, merges, and tag-strips input BAMs into `<work_dir>/inputs/`. **Do not hand
   `family.wdl` raw external BAM paths** — route new runs through `launch.sh`.

4. **Only large INS edits are currently supported** by the edit-QC analysis. The schema
   permits `SNV`/`DEL`, but those paths are not validated end-to-end — warn the user.

5. **`miniwdl.cfg` must exist before `launch.sh` will run.** `launch.sh` defaults to
   `<repo_root>/miniwdl.cfg` and aborts with "Default miniwdl.cfg not found" if it's missing —
   a very common first-run snag. There is **no repo-root cfg checked in**; the canonical
   SLURM/singularity template is `backends/hpc/miniwdl.cfg`. Recommend the user keep their
   own at `~/.config/miniwdl.cfg` (HPC-tuned), but an agent can always point launch.sh at the
   in-repo template directly: `./scripts/launch.sh <config> --miniwdl-cfg backends/hpc/miniwdl.cfg`
   (or `cp backends/hpc/miniwdl.cfg <repo>/miniwdl.cfg`). `launch.sh` copies the chosen cfg
   into the work dir and rewrites its cache paths — don't hand-edit the work-dir copy.

6. **Don't run the workflow directly in the foreground, and don't skip the container build.**
   `launch.sh` Phase 1 prepulls/builds all container images (via `create_image_manifest.sh` +
   `populate_miniwdl_singularity_cache.sh`) — never bypass it by calling `miniwdl run` on raw
   inputs. A family run takes many hours, so launch it inside a **`screen`** (or `tmux`)
   session so it survives logout; don't background it in a way that dies with the shell.

## Disk & paths

- **The work dir must live under the cloned repo** (e.g. `<repo>/workflow_run_*`).
  apptainer bind-mounts only paths under the cwd/run dir, so inputs and outputs both have
  to sit there (see guardrail 3). `launch.sh --work-dir` accepts a name (created under the
  repo) or a full path — keep it inside the repo tree.
- **Plan for large disk.** Staged inputs are ~50 GB per sample (`<work_dir>/inputs/`), plus
  ~3 GB × 2 more per sample when a precomputed assembly haplotype is supplied (see
  `publication_Ngn2-KOLF-may/inputs/`). Outputs scale with sample count — budget **~2 TB**
  for a family run. The cache dirs (`miniwdl_cache/`, `miniwdl_tmp/`) add more on top.
- **Share the container cache across runs.** `launch.sh --cache-dir` defaults to a *per-run*
  `<work_dir>/miniwdl_cache`, so each new run re-pulls all SIFs. Point `--cache-dir` at a
  single repo-root location (e.g. `<repo>/miniwdl_cache`) so the singularity image cache is
  shared and phase 1 is a near-no-op after the first run.
- **Archive, then free space — but only after the run completes.**
  - The run dir's cached task outputs let miniwdl skip completed work on resubmit, so
    **do not delete the work dir while a run is in progress or may be resumed.**
  - Once a family run finishes, run `scripts/archive.sh <run_dir> [archive_base_dir]` to
    copy the final files out (archive storage can be **anywhere** — it is not bind-mounted).
    Staged inputs and full outputs are too large to leave in the repo tree after completion.
  - The archive is smaller than the live run dir: it omits intermediate files that are only
    needed for partial/resumed runs. Use `archive.sh -d` to delete source outputs after the
    copy verifies.

## File map

- `workflows/family.wdl` — entrypoint. Edit-QC sub-workflows under `workflows/edit-qc/`.
- `workflows/edit-qc/crispr_edit.schema.json` — **authoritative** schema for expected-edit JSON.
- `GRCh38.{ref,tertiary,somatic}_map.v3p1p0.template.tsv` — reference/resource map templates.
  - Note: in `GRCh38.tertiary_map.v3p1p0.template.tsv` we set `slivar_max_af` to `1.00`
    (upstream default is `0.03`). This deliberately delays allele-frequency filtering until
    later, which gives cleaner variant-filtering statistics for publication. It should have
    no effect on the actual final reported variants.
- `genbank_to_crispr_json.py` — Benchling GenBank → expected-edit JSON converter.
- `scripts/setup.sh` — one-time reference + container download.
- `scripts/launch.sh` — stages inputs and (optionally `--run`) launches the workflow.
- `scripts/process_input_config.py` — BAM merge/strip/stage, called by launch.sh.
- `scripts/create_image_manifest.sh` / `populate_miniwdl_singularity_cache.sh` — container prepull.
- `my_inputs/` — real example input configs and expected-edit JSONs (use as references).
- `scripts/archive.sh` — copy a finished run's outputs to long-term archive storage.
- `docs/biohub-setup.md` — example HPC (SLURM) setup guide.

## Keeping the pipeline up to date

This repo uses git submodules (`workflows/wdl-common`). To check for and pull updates:

```bash
git -C <repo> fetch && git -C <repo> status   # see if behind the remote
git -C <repo> pull --recurse-submodules        # update, including submodules
git -C <repo> submodule update --init --recursive
```

After pulling, container image pins or reference versions may have changed — re-run
`scripts/create_image_manifest.sh` (launch.sh does this automatically) and check whether
`scripts/setup.sh` needs re-running for new reference files.

## Skills

- **`prepare-edit-inputs`** — build & validate a run's input config and expected-edit JSONs
  (GenBank conversion, ref-map placeholder fixes, dbNSFP licensing, assembly reuse). Invoke
  whenever a user is setting up a new run or preparing input data.
- **`run-and-monitor`** — launch via `launch.sh`, monitor logs sanely, and triage failures
  (when to auto-resubmit vs escalate). Invoke when running or debugging a workflow.
