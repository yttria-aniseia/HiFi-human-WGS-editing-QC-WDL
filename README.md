<h1 align="center"><img width="300px" src="https://github.com/PacificBiosciences/HiFi-human-WGS-WDL/blob/main/images/logo_wdl_workflows.svg" alt="PacBio WGS Variant Pipeline"/></h1>

<h1 align="center">HiFi Human WGS CRISPR Edit QC Pipeline</h1>

Workflow for analyzing human PacBio whole genome sequencing (WGS) data with CRISPR editing quality control using [Workflow Description Language (WDL)](https://openwdl.org/).

**This is a fork of the [PacBio HiFi-human-WGS-WDL pipeline](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL) customized for CRISPR editing experiments.** This fork adds edit-specific QC analysis while maintaining the core variant calling and analysis features.

- Docker images used by this workflow are defined in [the wdl-dockerfiles repo](../../../wdl-dockerfiles). Images are hosted in PacBio's [quay.io repo](https://quay.io/organization/pacbio).
- Common tasks that may be reused within or between workflows are defined in [the wdl-common repo](../../../wdl-common). Note: In this fork, `wdl-common` is not actively maintained.

## Workflow

**Only the `family` workflow is supported in this fork.** The family workflow is designed to analyze cohorts of related samples, which is ideal for CRISPR editing experiments where you typically have:
- Parental (wild-type) samples
- Edited clones derived from the parent
- Optional additional family relationships

The workflow analyzes human PacBio HiFi whole genome sequencing data and includes specialized analysis for CRISPR edits when `expected_edits` are provided in the input configuration.

**Workflow entrypoint**:
- [workflows/family.wdl](workflows/family.wdl)

**CRISPR Edit QC Features**:
- Detection and validation of expected edits (insertions, deletions, substitutions)
- Comparison of edited samples against parental baseline
- Integration with standard germline and somatic variant calling
- Support for complex edits including knock-ins with payload sequences

## Setup

This fork uses git submodules for common tasks. To clone the repository with submodules:

```bash
git clone --recurse-submodules --depth=1 \
  https://github.com/yttria-aniseia/HiFi-human-WGS-editing-QC-WDL.git
```

**For biohub-specific setup instructions** including conda environment setup, reference data download, and running the pipeline, see [docs/biohub-setup.md](docs/biohub-setup.md).

## Resource requirements

The most resource-heavy step in the workflow requires 64 cpu cores and 256 GB of RAM. Ensure that the backend environment you're using has enough quota to run the workflow.

On some backends, you may be able to make use of a GPU to accelerate the DeepVariant step.  The GPU is not required, but it can significantly speed up the workflow.  If you have access to a GPU, you can set the `gpu` parameter to `true` in the inputs JSON file.

## Reference datasets and associated workflow files

Reference datasets are hosted publicly for use in the pipeline. For data locations, see the backend-specific documentation and template inputs files for each backend with paths to publicly hosted reference files filled out.

## Setting up and executing the workflow

1. [Select a backend environment](#selecting-a-backend)
2. [Configure a workflow execution engine in the chosen environment](#configuring-a-workflow-engine-and-container-runtime)
3. [Fill out the inputs JSON file for your cohort](#filling-out-the-inputs-json)
4. [Run the workflow](#running-the-workflow)

### Selecting a backend

The workflow can be run on Azure, AWS, GCP, or HPC. Your choice of backend will largely be determined by the location of your data.

For backend-specific configuration, see the relevant documentation:

- [Azure](./docs/backend-azure.md)
- [AWS](./docs/backend-aws-healthomics.md)
- [GCP](./docs/backend-gcp.md)
- [HPC](./docs/backend-hpc.md)

### Configuring a workflow engine and container runtime

An execution engine is required to run workflows. Two popular engines for running WDL-based workflows are [`miniwdl`](https://miniwdl.readthedocs.io/en/latest/getting_started.html) and [`Cromwell`](https://cromwell.readthedocs.io/en/stable/tutorials/FiveMinuteIntro/).

Because workflow dependencies are containerized, a container runtime is required. This workflow has been tested with [Docker](https://docs.docker.com/get-docker/) and [Singularity](https://docs.sylabs.io/guides/3.10/user-guide/) container runtimes.

See the backend-specific documentation for details on setting up an engine.

| Engine | [Azure](./docs/backend-azure.md) | [AWS](./docs/backend-aws-healthomics.md) | [GCP](./docs/backend-gcp.md) | [HPC](./docs/backend-hpc.md) |
| :- | :- | :- | :- | :- |
| [**miniwdl**](https://github.com/chanzuckerberg/miniwdl#scaling-up) | _Unsupported_ | Supported via [AWS HealthOmics](https://aws.amazon.com/healthomics/) | _Unsupported_ | (SLURM only) Supported via the [`miniwdl-slurm`](https://github.com/miniwdl-ext/miniwdl-slurm) plugin |
| [**Cromwell**](https://cromwell.readthedocs.io/en/stable/backends/Backends/) | Supported via [Cromwell on Azure](https://github.com/microsoft/CromwellOnAzure) | _Unsupported_ | Supported via Google's [Pipelines API](https://cromwell.readthedocs.io/en/stable/backends/Google/) | Supported - [Configuration varies depending on HPC infrastructure](https://cromwell.readthedocs.io/en/stable/tutorials/HPCIntro/) |

### Filling out the inputs JSON

The input to a workflow run is defined in JSON format. Use [example_input_config.json](example_inputs/example_input_config.json) as a template.

Key steps:
1. Define your samples with their HiFi read BAM files
2. Specify family relationships (parent/child)
3. For CRISPR-edited samples, provide an `expected_edits` JSON file path
4. Point to your local reference map files (after running `./scripts/setup.sh`)

See [docs/biohub-setup.md](docs/biohub-setup.md) for detailed instructions.

**Automated workflow setup**: This fork includes a `launch.sh` script that automates file staging and workflow setup. See [scripts/README.md](scripts/README.md) for details.

### Running the workflow

**Recommended approach** using the automated launcher script:

```bash
# 1. Setup and stage inputs (with automatic workflow execution)
./scripts/launch.sh my_input_config.json --work-dir my_analysis_name --run

# Or setup only, run workflow manually later
./scripts/launch.sh my_input_config.json --work-dir my_analysis_name
conda activate hifi-wgs
bash my_analysis_name/run_workflow.sh
```

See [scripts/README.md](scripts/README.md) for detailed launcher script documentation.

#### Run directly using miniwdl (HPC with SLURM)

If not using the launcher script, you can run miniwdl directly:

```bash
miniwdl run --verbose \
  --cfg miniwdl.cfg \
  --dir output_directory \
  workflows/family.wdl \
  -i input_config.json
```

**Note**: This fork is primarily tested with miniwdl on HPC/SLURM environments. Support for other backends (Azure, GCP, AWS) may be limited.

## Workflow inputs

Workflow inputs for the family entrypoint are described in [family](./docs/family.md) documentation.

At a high level, we have two types of input files:

- **Map files** (TSV format) describe reference data and resources used for every workflow execution:
  - `ref_map_file`: Reference genome FASTA, indices, and core annotation files
  - `tertiary_map_file`: Population VCFs, SV databases, and tertiary analysis resources
  - `somatic_map_file`: Somatic variant calling resources

- **Input configuration JSON** describes the samples to analyze and their relationships:
  - Sample metadata (ID, sex, affected status, family relationships)
  - Paths to HiFi read BAM files
  - **For CRISPR editing experiments**: `expected_edits` field defining anticipated genomic changes

**Example input configuration**:
```json
{
  "humanwgs_family.family": {
    "family_id": "EXAMPLE_FAM",
    "samples": [
      {
        "sample_id": "parent",
        "hifi_reads": ["/path/to/parent.bam"],
        "sex": "FEMALE",
        "affected": false
      },
      {
        "sample_id": "edited_clone1",
        "hifi_reads": ["/path/to/clone1.bam"],
        "sex": "FEMALE",
        "affected": true,
        "mother_id": "parent",
        "expected_edits": "/path/to/edits.json"
      }
    ]
  },
  "humanwgs_family.ref_map_file": "/path/to/GRCh38.ref_map.tsv",
  "humanwgs_family.tertiary_map_file": "/path/to/GRCh38.tertiary_map.tsv",
  "humanwgs_family.somatic_map_file": "/path/to/GRCh38.somatic_map.tsv"
}
```

See [example_expected_edit.json](example_inputs/example_expected_edit.json) for the expected edit file format, or use the provided `genbank_to_crispr_json.py` helper script to generate edit descriptions from GenBank files.

The resource bundle containing the GRCh38 reference and other files used in this workflow can be downloaded from Zenodo:

[<img src="https://zenodo.org/badge/DOI/10.5281/zenodo.17086906.svg" alt="10.5281/zenodo.17086906">](https://zenodo.org/records/17086906)

Template map files are provided at the repository root: `GRCh38.ref_map.v3p1p0.template.tsv`, `GRCh38.tertiary_map.v3p1p0.template.tsv`, and `GRCh38.somatic_map.v3p1p0.template.tsv`. After downloading the reference bundle, update the paths in these templates to point to your local copies.

# Tool versions and Docker images

Docker images definitions used by this workflow can be found in [the wdl-dockerfiles repository](../../../wdl-dockerfiles/). Images are hosted in PacBio's [quay.io repo](https://quay.io/organization/pacbio). Docker images used in the workflow are pinned to specific versions by referring to their digests rather than tags.

The Docker image used by a particular step of the workflow can be identified by looking at the `docker` key in the `runtime` block for the given task. Images can be referenced in the following table by looking for the name after the final `/` character and before the `@sha256:...`. For example, the image referred to here is "align_hifiasm":
> ~{runtime_attributes.container_registry}/pb_wdl_base@sha256:4b889a1f ... b70a8e87

Tool versions and Docker images used in these workflows can be found in the [tools and containers](./docs/tools_containers.md) documentation.

---

## DISCLAIMER

TO THE GREATEST EXTENT PERMITTED BY APPLICABLE LAW, THIS WEBSITE AND ITS CONTENT, INCLUDING ALL SOFTWARE, SOFTWARE CODE, SITE-RELATED SERVICES, AND DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. ALL WARRANTIES ARE REJECTED AND DISCLAIMED. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THE FOREGOING. PACBIO IS NOT OBLIGATED TO PROVIDE ANY SUPPORT FOR ANY OF THE FOREGOING, AND ANY SUPPORT PACBIO DOES PROVIDE IS SIMILARLY PROVIDED WITHOUT REPRESENTATION OR WARRANTY OF ANY KIND. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A REPRESENTATION OR WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACBIO.
