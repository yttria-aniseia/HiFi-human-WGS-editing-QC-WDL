# Biohub HPC Setup Guide

Instructions for running the HiFi WGS CRISPR Edit QC pipeline on the Biohub HPC cluster.

## Initial Setup

### 1. Clone repository with submodules

```bash
git clone --recurse-submodules --depth=1 \
  git@github.com:yttria-aniseia/HiFi-human-WGS-editing-QC-WDL.git
cd HiFi-human-WGS-editing-QC-WDL
```

### 2. Create conda environment

```bash
conda env create -f environment.yml
conda activate hifi-wgs
```

### 3. Configure miniwdl

Create or edit `~/.config/miniwdl.cfg` for your HPC environment. See [docs/backend-hpc.md](backend-hpc.md) for SLURM-specific configuration.

### 4. Download reference data

This downloads ~200GB of reference files and will take several hours:

```bash
./scripts/setup.sh
```

This populates the reference map template files at the repository root with local paths.

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
