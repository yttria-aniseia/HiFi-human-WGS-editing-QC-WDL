---
name: run-and-monitor
description: Launch a HiFi WGS CRISPR edit-QC family run via launch.sh, monitor logs without drowning in them, and triage failures — deciding when a failed task is safe to auto-resubmit (container-copy thundering herd) versus when it needs human reconfiguration (OOM, time limit, real command errors). Use when running or debugging a workflow.
---

# Run and monitor an edit-QC workflow

## Launching

Inputs must be staged under the run dir (apptainer only bind-mounts paths under cwd) — always
go through `launch.sh`, never hand `family.wdl` raw external BAM paths.

```bash
# stage inputs + run
./scripts/launch.sh my_input_config.json --work-dir <name> --run
# or stage only, then run manually
./scripts/launch.sh my_input_config.json --work-dir <name>
conda activate hifi-wgs && bash <name>/run_workflow.sh
```

`launch.sh` phases: (1) prepull container images to the cache, (2) SLURM job that merges/strips/
stages BAMs via `process_input_config.py`, (3) write `run_workflow.sh` and optionally launch.

## Output directory structure (know this to avoid expensive find/grep)

```
<work_dir>/outputs/
  _LAST -> <timestamp>_humanwgs_family      # symlink to the most recent run; use it
  <timestamp>_humanwgs_family/
    outputs.json     # maps each workflow output name -> ABSOLUTE path (under out/ below)
    inputs.json, error.json, workflow.log
    wdl/             # snapshot of the WDL source used for this run
    out/             # CURATED final outputs (symlinks), organized by output name -- look HERE
    call-<task>-NN/  # hundreds of intermediate task dirs; scratch, not final products
```

- **Final results live in `out/` and are enumerated in `outputs.json`** — start there. Do
  **not** `find`/`grep` across the whole run dir (hundreds of `call-*` dirs, huge + slow).
- `outputs/_LAST` always points at the newest run.

## Monitoring (do NOT read logs in full)

Logs are huge and repetitive. Scope every look:

```bash
RUN=<work_dir>/outputs/_LAST
grep -E "ERROR|failed" "$RUN/workflow.log" | tail -40   # recent failures only
ls "$RUN" | grep '^call-'                                # which tasks ran
# for one failing task, read its small per-attempt stderr, not workflow.log:
tail -50 "$RUN/call-<task>/stderr.txt"
```

Useful task-dir contents: `error.json`, `slurm-*.out`, `slurm_singularity.log.txt`,
`stderr.txt`/`stderr2.txt`… (one per retry attempt), `command`.

## Failure triage — auto-resubmit vs escalate

Decide by **signature**, not by the bare exit code (most failures surface as
`exit status 255`, which alone is uninformative).

### MOST COMMON container failure — build temp dir out of space

When apptainer pulls an image it builds/converts a large temporary SIF in
`APPTAINER_TMPDIR` / `SINGULARITY_TMPDIR`. If that dir is on a small filesystem (or the
default `/tmp`, or a nearly-full scratch), the build fails. Signature in
`slurm_singularity.log.txt` / `slurm-*.out` / `workflow.log`:

```
FATAL: ... While building SIF / creating squashfs / mksquashfs ... : No space left on device
# or: write error, disk quota exceeded, failed to create build directory
```

**This is a reconfiguration, not a retry** — rerunning with the same temp dir fails identically.
Point the build temp + cache at a roomy location (scratch), then rerun:

```bash
# launch.sh flags (preferred — keeps everything under one work dir on scratch):
./scripts/launch.sh <config> --work-dir <name> \
    --cache-dir /hpc/scratch/<user>/cache --tmp-dir /hpc/scratch/<user>/tmp
# or set env before a manual pull/run:
export APPTAINER_TMPDIR=/hpc/scratch/<user>/apptainer_tmp
export SINGULARITY_TMPDIR=$APPTAINER_TMPDIR
export APPTAINER_CACHEDIR=/hpc/scratch/<user>/cache/download_cache
export SINGULARITY_CACHEDIR=$APPTAINER_CACHEDIR
```

`launch.sh` defaults these to `<work_dir>/miniwdl_cache/...`, so the fix is usually putting the
**work dir on a filesystem with tens of GB free**, not on a small home/quota volume. Check free
space with `df -h <tmp-dir>` before resubmitting.

**Incremental workaround when temp space is tight:** completed builds land in the
**cache** (`singularity_cache/`, the keepers); the **temp dir** only holds in-progress build
scratch. So you can let it build as many containers as fit, then **clear the build temp dir
(NOT the cache) and rerun** — already-built images are skipped and the freed temp space is
reused for the next batch. Repeat until all images are cached:

```bash
rm -rf "$APPTAINER_TMPDIR"/* "$SINGULARITY_TMPDIR"/*   # safe: only build scratch
# leave singularity_cache/ alone — those are the finished SIFs
bash <work_dir>/run_workflow.sh                        # or rerun the prepull
```

### SAFE TO AUTO-RESUBMIT — container-copy "thundering herd"

Signature in `workflow.log` — many tasks failing within the same second with:

```
ERROR slurm_singularity pull failed :: stderr: ["FATAL: While making image from oci registry:
error fetching image to cache: failed to get checksum for docker://<img>:
pinging container registry registry-1.docker.io: ... context canceled"]
```

Cause: many tasks try to pull the **same uncached image** from docker.io at once; the registry
cancels/rate-limits. **Action:** simply rerun — miniwdl call-caching skips completed work, and
the image is usually cached after the first success. If it recurs, **prepull the image first**:

```bash
bash scripts/populate_miniwdl_singularity_cache.sh image_manifest.txt <cache>/singularity_cache
```

**Root-cause fix (preferred over retrying):** `scripts/create_image_manifest.sh` only captures
`@sha256`-pinned images. **Tag-pinned images are missing from the manifest**, so they're never
prepulled and always pull live → this failure. Current offenders (grep
`'docker:.*"' workflows/ | grep -v @sha256` to refresh this list):
`quay.io/biocontainers/bcftools:1.17--h3cc50cf_1`,
`quay.io/biocontainers/pybigwig:0.3.24--py311hd8c7dd8_0` (the `knock-knock` local image is
built separately by `setup.sh`). Add such tags to `image_manifest.txt` and prepull. Don't rely
on retries to paper over it.

### MUST ESCALATE / RECONFIGURE — do not blindly retry

- **Out of memory** — `slurmstepd: error: ... oom-kill`, `Killed`, `MemoryError`,
  `std::bad_alloc`, `Cannot allocate memory` in `slurm-*.out`/stderr, or a SLURM exit
  `OUT_OF_MEMORY`. Retrying as-is will fail identically. Requires raising the task's `mem_gb`
  in the WDL `runtime`/`RuntimeAttributes`, or running on a larger partition. Flag to the user.
- **Time limit** — `DUE TO TIME LIMIT`, `TIMEOUT`, `CANCELLED ... DUE TO TIME`. Needs a longer
  `--time` for that task. Escalate.
- **Disk / no space** — `No space left on device`, squashfs/extract failures. Cache or `--tmp`
  is full; escalate.
- **Real command error** — a non-empty `stderr.txt` with a tool-specific traceback/usage error
  (bad input, malformed JSON, missing reference). Read that stderr; fix the input/config. A
  retry won't help. This is where input-prep mistakes surface — re-check with `prepare-edit-inputs`.

### Decision rule

1. `grep ERROR workflow.log | tail` — are failures clustered in one second across many tasks
   with `slurm_singularity pull failed`? → thundering herd → **resubmit / prepull**.
2. Otherwise open the failing task's `stderr.txt` + `slurm-*.out` and match against the
   escalate signatures above. **If it's OOM / time / disk / real error, STOP and report to the
   user with the signature and the specific task** — do not loop retries.

Note: miniwdl already retries each task up to `max_retries` (2–3). If a task exhausted its
retries and the workflow aborted, an automatic in-band retry has *already failed* — only a
fresh `bash run_workflow.sh` (relying on call-caching) or a config change will help.

## Archiving a finished run (outputs get large)

Run dirs are large (hundreds of intermediate `call-*` dirs). Once a run is complete and
verified, archive only the curated outputs + metadata to long-term storage with
`scripts/archive.sh`:

```bash
# archives <run_dir>/outputs/_LAST  ->  <archive_base>/<run_dir>/
scripts/archive.sh <run_dir> <archive_base_dir>
scripts/archive.sh -i <run_dir> <archive_base_dir>   # -i: prompt before archive/delete
scripts/archive.sh -d <run_dir> <archive_base_dir>   # -d: delete source outputs/ after copy
```

It copies `_LAST/out` + `_LAST/wdl` and the `inputs.json`/`outputs.json`/`error.json`/
`workflow.log` metadata, and rewrites `expected_edit` paths in the copied `inputs.json` to
relative (`./<file>.json`), keeping the original as `inputs.orig.json`. Confirm the archive
base path with the user — there is no hardcoded default.

### Absolute paths break when outputs move — rebase them

`outputs.json` records every output as an **absolute** path under the original
`<work_dir>/outputs/<run>/out/...`. After archiving/moving, those paths no longer resolve.
The stable part is everything from `/out/...` onward; reinterpret it relative to the new
location. To consume a moved run, rebase the prefix:

```bash
# original prefix (up to and including the run dir) -> new archive location
sed 's|/orig/path/<work_dir>/outputs/_LAST|<archive_base>/<work_dir>|g' \
    <archive_base>/<work_dir>/outputs.json > outputs.rebased.json
# or, in code: strip everything through ".../out/" and join onto <archive>/.../out/
```

When pointing the downstream report repo (or any consumer) at an archived run, use the
rebased paths — the raw `outputs.json` from a moved run will have dangling absolute paths.
