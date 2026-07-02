# Biohub HPC Setup Guide

Instructions for running the HiFi WGS CRISPR Edit QC pipeline on the Biohub HPC cluster.

> **Using an AI coding agent?** This repo ships agent-facing docs: `AGENTS.md` (orientation +
> guardrails, also imported by `CLAUDE.md`) and two skills under `.claude/skills/` —
> `prepare-edit-inputs` (build/validate input configs and expected-edit JSONs) and
> `run-and-monitor` (launch, monitor logs, triage failures, archive outputs). Point your agent
> at those when onboarding a new run.

## Initial Setup

### 1. Clone repository with submodules

```bash
git clone --recurse-submodules --depth=1 \
  https://github.com/yttria-aniseia/HiFi-human-WGS-editing-QC-WDL.git
cd HiFi-human-WGS-editing-QC-WDL
```

### 2. Create conda environment

```bash
conda env create -f environment.yml
conda activate hifi-wdl
```

### 3. Configure miniwdl

Create or edit `~/.config/miniwdl.cfg` for your HPC environment. See [docs/backend-hpc.md](backend-hpc.md) for SLURM-specific configuration.

### 4. Download reference data and build containers

Before running setup, collect two prerequisites:

#### 4a. Set container cache location

`setup.sh` builds the knock-knock Singularity image into `$SINGULARITY_CACHEDIR`. On Biohub
HPC the default home-directory cache is too small — point it at scratch storage first:

```bash
export SINGULARITY_CACHEDIR="$(pwd)/miniwdl_cache/singularity_cache"
mkdir -p "${SINGULARITY_CACHEDIR}"
```

If `SINGULARITY_CACHEDIR` is unset, `setup.sh` will skip the knock-knock build with a warning
and you will need to re-run it later with the variable set.

#### 4b. Obtain dbNSFP (license-gated)

dbNSFP is required for somatic variant annotation but is not in the bundle due to licensing.
Request a licensed copy from [https://www.dbnsfp.org/download](https://www.dbnsfp.org/download).
You need the **GRCh38 BGZF files** listed under _"dbNSFP variants in BGZF format for VEP and
SnpEff annotation programs (sorted by GRCh38 and GRCh37 coordinates)"_:

- `dbNSFP5.3.1a_grch38.gz`
- `dbNSFP5.3.1a_grch38.gz.tbi`

Have the path to `dbNSFP5.3.1a_grch38.gz` ready before running setup — passing it via
`--dbnsfp` avoids having to re-run the full download a second time. If you cannot obtain
dbNSFP yet, setup will complete with placeholders in the somatic map.

#### 4c. Run setup

Downloads ~46 GB from Zenodo (resumable; typically 15–30 min), writes map files with
absolute paths, and builds the knock-knock container:

```bash
# Recommended: pass dbNSFP in the same invocation
./scripts/setup.sh --dbnsfp /path/to/dbNSFP5.3.1a_grch38.gz

# Without dbNSFP (somatic map will have placeholders — patch later with fetch_resources.sh --dbnsfp):
./scripts/setup.sh
```

## Preparing Input Files

### Input configuration JSON

Create a JSON file describing your samples. Start from [example_input_config.json](../example_input_config.json):

```json
{
  "humanwgs_family.family": {
    "family_id": "my_experiment",
    "samples": [
      {
        "sample_id": "parent_wt",
        "hifi_reads": ["/path/to/parent.bam"],
        "sex": "FEMALE",
        "affected": false
      },
      {
        "sample_id": "clone_edited",
        "hifi_reads": ["/path/to/clone.bam"],
        "sex": "FEMALE",
        "affected": true,
        "mother_id": "parent_wt",
        "expected_edits": "/path/to/edits.json"
      }
    ]
  },
  "humanwgs_family.ref_map_file": "GRCh38.ref_map.v3p1p0.template.tsv",
  "humanwgs_family.tertiary_map_file": "GRCh38.tertiary_map.v3p1p0.template.tsv",
  "humanwgs_family.somatic_map_file": "GRCh38.somatic_map.v3p1p0.template.tsv"
}
```

### Expected edits JSON (for CRISPR experiments)

Specify expected genomic changes for edited samples. See [example_expected_edit.json](../example_expected_edit.json) for format.

For complex knock-ins, use `genbank_to_crispr_json.py` to generate from GenBank files:

```bash
python3 genbank_to_crispr_json.py \
  -i my_plasmid.gb \
  -o my_edits.json \
  -e my_edit_name \
  -t ENST00000123456 \
  -g GGGCCACATCAGCGCGATCC \
  --type INS
```

Arguments:
- `-i`: Input GenBank file
- `-o`: Output JSON file
- `-e`: Edit identifier
- `-t`: Target Ensembl transcript ID
- `-g`: Guide RNA sequence
- `--type`: Edit type (INS, DEL, or SUB)

### CNV ploidy configuration (optional)

For non-standard ploidy expectations, set `sawfish_expected_bed_male` and/or `sawfish_expected_bed_female` in the ref_map TSV. See [sawfish expected CN files](https://github.com/PacificBiosciences/sawfish/tree/main/data/expected_cn).

## Running the Pipeline

### Using launch.sh (recommended)

```bash
./scripts/launch.sh my_inputs.json --work-dir my_analysis_name
```

This stages inputs, merges BAM files, strips tags, and creates a run script. Then:

```bash
conda activate hifi-wgs
bash my_analysis_name/run_workflow.sh
```

Or run everything at once with `--run` flag:

```bash
./scripts/launch.sh my_inputs.json --work-dir my_analysis_name --run
```

The workflow runs for ~24 hours depending on sample count and relatedness.

### Manual execution

```bash
conda activate hifi-wgs
miniwdl run --verbose \
  --cfg miniwdl.cfg \
  --dir output_dir \
  workflows/family.wdl \
  -i my_inputs.json
```

## Generating Reports

After workflow completion, generate QC reports using the [hifi-wgs-qc-report](https://github.com/your-org/hifi-wgs-qc-report) repository:

```bash
cd /path/to/hifi-wgs-qc-report
source venv/bin/activate

RUNNAME="my_analysis"
RESULTS="../HiFi-human-WGS-editing-QC-WDL/${RUNNAME}/outputs/_LAST/"
REPORTDIR="/hpc/websites/onsite.czbiohub.org/compbio/hifi-wgs-qc-reports/${RUNNAME}_$(date +%Y%m%d)"

python3 generate_reports.py -o "${REPORTDIR}" "${RESULTS}"
```

## Archiving Results

Move completed analysis to archive storage:

```bash
./scripts/archive.sh <my_analysis_name> <archive_dir>
```

## Troubleshooting

**Environment issues**: Verify conda environment installation and that miniwdl is available.

**Container cache issues**: Use a shared cache directory to reuse containers across runs:

```bash
./scripts/launch.sh my_inputs.json --work-dir my_analysis --cache-dir /scratch/shared/containers
```

**BAM merging failures**: Check that samtools is installed and input BAM files are accessible.

**Workflow failures**: Check logs in `work_dir/logs/` and miniwdl output for error details.
