# Installing and configuring for HPC backends

> **Note**: This fork only supports the `family` workflow. The instructions below are for running the family workflow on HPC using SLURM.

Either `miniwdl` or `Cromwell` can be used to run workflows on the HPC.

## Installing and configuring `miniwdl`

### Requirements

- [`miniwdl`](https://github.com/chanzuckerberg/miniwdl) >= 1.9.0
- [`miniwdl-slurm`](https://github.com/miniwdl-ext/miniwdl-slurm)

### Configuration

An example miniwdl.cfg file is provided in `backends/hpc/miniwdl.cfg`. This should be placed at `~/.config/miniwdl.cfg` and edited to match your SLURM configuration.

> [!IMPORTANT]
> In order to simplify workflow inputs, we make use of `map` files to specify the input data. This allows for a more concise input file, but requires changing a miniwdl configuration option to allow workflows to access files that are not expressly supplied with workflow inputs.  To enable this, add the following line to your `miniwdl.cfg` file:
>
> ```ini
> [file_io]
> allow_any_input = true
> ```
>
> This option is already included in the example miniwdl.cfg file described in this section.

## Installing and configuring `Cromwell`

Cromwell supports a number of different HPC backends; see [Cromwell's documentation](https://cromwell.readthedocs.io/en/stable/backends/HPC/) for more information on configuring each of the backends.  Cromwell can be used in a standalone "run" mode, or in "server" mode to allow for multiple users to submit workflows.  In the example below, we provide example commands for running Cromwell in "run" mode.

> [!NOTE]
> If running Cromwell on an HPC cluster using NFS for storage, you may encounter issues with NFS latency, which can cause Cromwell to fail to read files from the filesystem.  You can work around this by using `script-epilogue` to add a delay & sync to the end of each job.  This option is added to your backend provider config.
>
> ```bash
> script-epilogue = "sleep 60 && sync"
> ```

## Running the workflow

### Filling out workflow inputs

Create an input configuration JSON describing your samples and their relationships. See [example_input_config.json](../example_input_config.json) as a template. After downloading reference data with `./scripts/setup.sh`, the populated map files (`GRCh38.{ref,tertiary,somatic}_map.v3p1p0.template.tsv`) at the repository root will hold absolute local paths into the frozen bundle.

See [family.md](./family.md) for input structure details, or [biohub-setup.md](./biohub-setup.md) for biohub-specific instructions.

### Recommended: Using launch.sh

The automated launcher script handles file staging and setup:

```bash
./scripts/launch.sh my_inputs.json --work-dir my_analysis_name
conda activate hifi-wgs
bash my_analysis_name/run_workflow.sh
```

Or run everything at once:

```bash
./scripts/launch.sh my_inputs.json --work-dir my_analysis_name --run
```

See [scripts/README.md](../scripts/README.md) for details.

### Manual execution

#### Running via miniwdl

```bash
miniwdl run workflows/family.wdl --input <inputs_json_file>
```

If compute nodes cannot reach the internet, use `./scripts/populate_miniwdl_singularity_cache.sh` with `./image_manifest.txt` to pre-pull container images from a login node.

#### Running via Cromwell

```bash
cromwell run workflows/family.wdl --input <inputs_json_file>
```

## Reference data bundle

[<img src="https://zenodo.org/badge/DOI/10.5281/zenodo.20856447.svg" alt="10.5281/zenodo.20856447">](https://zenodo.org/records/20856447)

Reference data (~46 GB) is hosted on Zenodo. The download is resumable and checksum-verified; it typically takes 15–30 minutes depending on network speed.

Before running setup, collect two prerequisites:

**1. Set `SINGULARITY_CACHEDIR`** to a location with sufficient space (the knock-knock container
is several GB). If unset, `setup.sh` skips the knock-knock build with a warning. A good default
is a directory under the repo root on scratch:

```bash
export SINGULARITY_CACHEDIR="$(pwd)/miniwdl_cache/singularity_cache"
mkdir -p "${SINGULARITY_CACHEDIR}"
```

**2. Obtain dbNSFP** (license-gated, not in bundle). dbNSFP is required for somatic variant
annotation. Request a licensed copy from [https://www.dbnsfp.org/download](https://www.dbnsfp.org/download).
You need the **GRCh38 BGZF files** listed under _"dbNSFP variants in BGZF format for VEP and
SnpEff annotation programs (sorted by GRCh38 and GRCh37 coordinates)"_:

- `dbNSFP5.3.1a_grch38.gz`
- `dbNSFP5.3.1a_grch38.gz.tbi`

Have the path ready before running setup — passing it via `--dbnsfp` avoids re-running the
full download a second time. If you cannot obtain dbNSFP yet, setup will complete with
placeholders in the somatic map.

Then download the bundle (~46 GB, resumable, typically 15–30 min) and build the knock-knock container:

```bash
# Recommended: pass dbNSFP in the same invocation
./scripts/setup.sh --dbnsfp /path/to/dbNSFP5.3.1a_grch38.gz

# Without dbNSFP (somatic map will have placeholders):
./scripts/setup.sh

# Resources only, no knock-knock build:
./scripts/fetch_resources.sh --dbnsfp /path/to/dbNSFP5.3.1a_grch38.gz
```

All commands extract the bundle under the repo root and write ready-to-use map files (`GRCh38.{ref,tertiary,somatic}_map.v3p1p0.template.tsv`) with absolute paths. If a download is interrupted, re-run — it resumes from where it stopped.
