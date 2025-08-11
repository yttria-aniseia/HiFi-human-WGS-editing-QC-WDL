# HiFi Human WGS Workflow Launch Scripts

This directory contains scripts to help launch and manage the HiFi human WGS family workflow.

## launch.sh - Automated Workflow Launcher

The `launch.sh` script automates the entire setup process for running the HiFi human WGS family workflow. It eliminates the need for manual file copying, JSON editing, and directory setup.

### Quick Start

```bash
# Basic usage with default settings
./scripts/launch.sh your_input_config.json

# Custom workflow directory name
./scripts/launch.sh your_input_config.json my_analysis_name

# Full customization
./scripts/launch.sh your_input_config.json my_analysis hifi-wgs-env
```

### Usage

```bash
./scripts/launch.sh <input_config_json> [work_dir_name] [conda_env]
```

**Arguments:**
- `input_config_json`: JSON configuration file with sample information and file paths
- `work_dir_name`: Optional name for the working directory (default: `workflow_run_YYYYMMDD_HHMMSS`)
- `conda_env`: Optional conda environment name (default: `hifi-wgs`)

### What the Script Does

1. **Environment Setup**: Detects and activates conda/mamba environment
2. **Directory Creation**: Creates organized work directory structure:
   ```
   work_directory/
   ├── inputs/           # Copied input files
   ├── outputs/          # Workflow outputs
   ├── logs/             # Execution logs
   ├── workflows/        # WDL workflow files
   ├── backends/         # Backend configurations
   ├── miniwdl.cfg       # MiniWDL configuration
   ├── family.hpc.inputs.json  # Generated inputs file
   └── run_workflow.sh   # Executable run script
   ```
3. **File Management**: Copies all input files (HiFi reads, reference files) to local `inputs/` directory
4. **JSON Generation**: Creates new `family.hpc.inputs.json` with updated local file paths
5. **Container Setup**: Configures Singularity/Apptainer cache directories
6. **Run Script Creation**: Generates executable `run_workflow.sh` for easy workflow execution

### Input Configuration File

Create a JSON file similar to the existing `backends/hpc/family.hpc.inputs.v2.json`. Use `example_input_config.json` as a template:

```json
{
  "humanwgs_family.family": {
    "family_id": "MY_FAMILY",
    "samples": [
      {
        "sample_id": "sample1",
        "hifi_reads": [
          "/path/to/your/sample1_reads.bam"
        ],
        "sex": "MALE",
        "affected": false,
        "father_id": "NA",
        "mother_id": "NA"
      },
      {
        "sample_id": "sample2", 
        "hifi_reads": [
          "/path/to/your/sample2_reads.bam"
        ],
        "sex": "FEMALE",
        "affected": true,
        "father_id": "sample1",
        "mother_id": "NA"
      }
    ]
  },
  "humanwgs_family.phenotypes": "HP:0000001,HP:0000002",
  "humanwgs_family.ref_map_file": "/path/to/GRCh38.ref_map.v2p0p0.template.tsv",
  "humanwgs_family.tertiary_map_file": "/path/to/GRCh38.tertiary_map.v2p0p0.template.tsv",
  "humanwgs_family.pharmcat_min_coverage": 1,
  "humanwgs_family.backend": "HPC",
  "humanwgs_family.preemptible": true
}
```

### Work Directory Naming Examples

**Default timestamped name:**
```bash
./scripts/launch.sh config.json
# Creates: ./workflow_run_20241211_143022
```

**Simple descriptive name:**
```bash
./scripts/launch.sh config.json trio_analysis
# Creates: ./trio_analysis
```

**Project-specific naming:**
```bash
./scripts/launch.sh config.json KOLF2_family_Dec2024
# Creates: ./KOLF2_family_Dec2024
```

**Full path specification:**
```bash
./scripts/launch.sh config.json /scratch/user/my_project/analysis1
# Creates: /scratch/user/my_project/analysis1
```

### Running the Workflow

After the launch script completes, you have two options to run the workflow:

**Option 1: Use the generated run script (recommended)**
```bash
cd your_work_directory
./run_workflow.sh
```

**Option 2: Run directly with miniwdl**
```bash
cd your_work_directory
miniwdl run --cfg miniwdl.cfg --dir outputs workflows/family.wdl -i family.hpc.inputs.json --verbose
```

### Monitoring Progress

Monitor workflow execution:
```bash
# Follow the log file
tail -f logs/workflow_*.log

# Check miniwdl status
miniwdl status outputs/

# View outputs
ls -la outputs/
```

### Prerequisites

1. **Conda Environment**: Ensure you have a conda environment with miniwdl installed
   ```bash
   conda create -n hifi-wgs miniwdl
   conda activate hifi-wgs
   ```

2. **Container Runtime**: Singularity/Apptainer must be available on your system

3. **Input Files**: Ensure all file paths in your input configuration JSON are accessible

4. **Reference Data**: Download the reference bundle from Zenodo if using local paths

### Troubleshooting

**Environment not found:**
```bash
# List available environments
conda env list

# Create environment if needed
conda create -n hifi-wgs miniwdl
```

**File not found errors:**
- Verify all file paths in your input JSON exist and are accessible
- Check file permissions for read access

**Container cache issues:**
- Ensure sufficient disk space for Singularity cache
- Check that Singularity/Apptainer is properly installed

**Workflow failures:**
- Check the log files in `logs/` directory
- Review miniwdl output for specific error messages
- Verify resource requirements are met (64 CPU cores, 256 GB RAM for peak usage)

### Example Complete Workflow

```bash
# 1. Create your input configuration
cp example_input_config.json my_trio_config.json
# Edit my_trio_config.json with your file paths

# 2. Launch the workflow setup
./scripts/launch.sh my_trio_config.json trio_analysis_2024

# 3. Run the workflow
cd trio_analysis_2024
./run_workflow.sh

# 4. Monitor progress
tail -f logs/workflow_*.log
```

## Other Scripts

### test_run.sh

Example script showing direct miniwdl execution. Used for testing and development.

### sample.sh

Helper script for sample-specific operations (if present).

---

For more information about the workflow itself, see the main [README.md](../README.md) and documentation in the [docs/](../docs/) directory.