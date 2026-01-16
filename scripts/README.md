# HiFi Human WGS Workflow Launch Scripts

This directory contains scripts to help launch and manage the HiFi human WGS family workflow.

## launch.sh - Automated Workflow Launcher

The `launch.sh` script automates the entire setup process for running the HiFi human WGS family workflow. It eliminates the need for manual file copying, JSON editing, and directory setup while providing intelligent file management and flexible configuration options.

### Quick Start

```bash
# Basic usage with default settings
./scripts/launch.sh your_input_config.json

# Custom workflow directory name
./scripts/launch.sh your_input_config.json --work-dir my_analysis_name

# With custom cache locations (recommended for HPC)
./scripts/launch.sh your_input_config.json --work-dir my_analysis --cache-dir /scratch/user/cache --tmp-dir /scratch/user/tmp

# With custom miniwdl.cfg
./scripts/launch.sh your_input_config.json --work-dir my_analysis --miniwdl-cfg /path/to/custom/miniwdl.cfg

# Auto-start workflow after setup
./scripts/launch.sh your_input_config.json --work-dir my_analysis --run

# Full customization
./scripts/launch.sh your_input_config.json --work-dir my_analysis --cache-dir /scratch/cache --tmp-dir /scratch/tmp --miniwdl-cfg /custom/miniwdl.cfg --run
```

### Usage

```bash
./scripts/launch.sh <input_config_json> [options]
```

**Required Arguments:**
- `input_config_json`: JSON configuration file with sample information and file paths

**Options:**
- `--work-dir <path>`: Working directory path (default: `workflow_run_YYYYMMDD_HHMMSS`)
- `--cache-dir <path>`: Apptainer cache directory (default: `<work_dir>/miniwdl_cache`)
- `--tmp-dir <path>`: Apptainer temp directory (default: `<work_dir>/miniwdl_tmp`)
- `--miniwdl-cfg <path>`: Path to custom miniwdl.cfg (default: `<repo_root>/miniwdl.cfg`)
- `--run`: Automatically start workflow execution after setup
- `--help`: Show help message

### What the Script Does

The script automates workflow setup in three phases:

**Phase 1: Setup (runs in current shell)**
1. Parse arguments and validate inputs
2. Create cache directories and export Apptainer environment variables
3. Pre-pull container images to the Singularity/Apptainer cache
4. Create work directory with `inputs/`, `outputs/`, and `logs/` subdirectories
5. Copy and update miniwdl.cfg with cache paths

**Phase 2: Input Processing (runs in SLURM job - 8 CPUs, 16GB, 4 hours)**
1. Merge multiple HiFi BAM files per sample using `samtools cat`
2. Strip phasing tags (fi, ri, fp, rp, ip, pw, HP, PS, PC) using `samtools reset`
3. Copy reference map files (ref_map, somatic_map, tertiary_map) to inputs directory
4. Copy expected_edits JSON file to inputs directory
5. Generate `family.hpc.inputs.json` with local file paths

**Phase 3: Workflow Execution (optional, if --run flag provided)**
1. Start miniwdl workflow automatically
2. Or create `run_workflow.sh` script for manual execution

The script skips BAM processing if the output file exists and passes `samtools quickcheck`.

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
./scripts/launch.sh config.json --work-dir trio_analysis
# Creates: ./trio_analysis
```

**Project-specific naming:**
```bash
./scripts/launch.sh config.json --work-dir KOLF2_family_Dec2024
# Creates: ./KOLF2_family_Dec2024
```

**Full path specification:**
```bash
./scripts/launch.sh config.json --work-dir /scratch/user/my_project/analysis1
# Creates: /scratch/user/my_project/analysis1
```

### Running the Workflow

After the launch script completes, you need to activate your conda environment and run the workflow **from the repository root directory** to avoid miniwdl I/O issues:

**Option 1: Use the generated run script (recommended)**
```bash
conda activate hifi-wgs  # or your environment name
cd /path/to/repository/root
bash your_work_directory/run_workflow.sh
```

**Option 2: Run directly with miniwdl**
```bash
conda activate hifi-wgs  # or your environment name
cd /path/to/repository/root
miniwdl run --cfg your_work_directory/miniwdl.cfg --dir your_work_directory/outputs workflows/family.wdl -i your_work_directory/family.hpc.inputs.json --verbose
```

**Important Notes:**
- The workflow must be run from the **repository root directory**, not the work directory
- The run script automatically handles all path prefixes for you
- This approach prevents miniwdl root directory I/O errors

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

1. **Conda Environment**: Ensure you have a conda environment with miniwdl and samtools installed
   ```bash
   conda env create -f environment.yml
   conda activate hifi-wgs
   ```

2. **Container Runtime**: Singularity/Apptainer must be available on your system

3. **Samtools**: Required for BAM merging and tag stripping (included in conda environment)

4. **Input Files**: Ensure all file paths in your input configuration JSON are accessible

5. **Reference Data**: Download the reference bundle from Zenodo if using local paths
   ```bash
   ./scripts/setup.sh
   ```

### Advanced Usage Examples

**Shared cache directory for multiple runs:**
```bash
# Use shared cache to avoid re-downloading containers across runs
./scripts/launch.sh config.json --work-dir analysis --cache-dir /scratch/shared/container_cache --tmp-dir /scratch/user/tmp
```

**Custom miniwdl configuration:**
```bash
# Use your own miniwdl.cfg with specific cluster settings
./scripts/launch.sh config.json --work-dir analysis --miniwdl-cfg /path/to/my/custom.cfg
```

**Auto-start workflow after setup:**
```bash
# Setup and immediately start workflow execution
./scripts/launch.sh config.json --work-dir analysis --run
```

**Re-running with same data:**
```bash
# Second run skips file copying (files pass samtools quickcheck)
./scripts/launch.sh same_config.json --work-dir analysis_v2
# Output: "Skipping merge for sample1 -> inputs/sample1_hifi_reads_merged.bam (already exists)"
```

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
- Ensure sufficient disk space for Apptainer cache (can be very large)
- Check that Singularity/Apptainer is properly installed
- Consider using custom cache location on high-capacity storage

**Workflow failures:**
- Check the log files in `work_directory/logs/` directory
- Review miniwdl output for specific error messages
- Verify resource requirements are met (64 CPU cores, 256 GB RAM for peak usage)
- Ensure you're running from repository root directory

**miniwdl I/O errors:**
- Make sure to run the workflow from the repository root directory
- Use `bash work_directory/run_workflow.sh` instead of running from within work directory
- Check that all paths in the generated run script are correct

**Cache path issues:**
- Verify cache directories have write permissions
- Check disk space in cache locations
- Ensure cache paths are accessible from compute nodes (in cluster environments)

### Cache Directory Configuration

The script supports flexible cache directory configuration for reusing container images across runs:

**Default (per-run cache):**
```bash
./scripts/launch.sh config.json --work-dir my_analysis
# Creates cache at: my_analysis/miniwdl_cache/
```

**Shared cache location (recommended for multiple runs):**
```bash
./scripts/launch.sh config.json --work-dir my_analysis --cache-dir /scratch/shared/containers --tmp-dir /scratch/user/tmp
# Reuses container images across different workflow runs
```

**Cache Structure Created:**
```
cache_location/
├── call_cache/           # miniwdl call caching (workflow execution)
├── cromwell_db/          # Cromwell database
├── download_cache/       # Apptainer image downloads (APPTAINER_CACHEDIR)
├── singularity_cache/    # Singularity/Apptainer images (reused across runs)
└── tmp/                  # Apptainer temp operations (APPTAINER_TMP_DIR)
```

### Smart File Copying

The script includes intelligent file management to avoid redundant operations:

**Features:**
- **Duplicate Detection**: Compares file size and MD5 hash
- **Skip Identical Files**: Won't re-copy files that already exist and are identical
- **Progress Feedback**: Clear messages about copied vs skipped files
- **Time & Space Savings**: Especially beneficial for large BAM files

**Example Output:**
```
Copying /data/sample1.bam -> inputs/sample1_hifi_reads_0.bam
Skipping /data/sample2.bam -> inputs/sample2_hifi_reads_0.bam (already exists and identical)
Copying /data/ref_map.tsv -> inputs/ref_map_file.tsv
```

### Example Complete Workflow

```bash
# 1. Create your input configuration
cp example_input_config.json my_trio_config.json
# Edit my_trio_config.json with your file paths

# 2. Launch the workflow setup with shared cache location
./scripts/launch.sh my_trio_config.json --work-dir trio_analysis_2024 --cache-dir /scratch/shared/containers

# 3. Run the workflow from repository root
conda activate hifi-wgs
cd /path/to/HiFi-human-WGS-editing-QC-WDL  # Repository root
bash trio_analysis_2024/run_workflow.sh

# 4. Monitor progress
tail -f trio_analysis_2024/logs/workflow_*.log
```

## Other Scripts

### test_run.sh

Example script showing direct miniwdl execution. Used for testing and development.

### sample.sh

Helper script for sample-specific operations (if present).

---

For more information about the workflow itself, see the main [README.md](../README.md) and documentation in the [docs/](../docs/) directory.