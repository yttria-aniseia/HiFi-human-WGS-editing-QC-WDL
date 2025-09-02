# HiFi Human WGS Workflow Launch Scripts

This directory contains scripts to help launch and manage the HiFi human WGS family workflow.

## launch.sh - Automated Workflow Launcher

The `launch.sh` script automates the entire setup process for running the HiFi human WGS family workflow. It eliminates the need for manual file copying, JSON editing, and directory setup while providing intelligent file management and flexible configuration options.

### Quick Start

```bash
# Basic usage with default settings
./scripts/launch.sh your_input_config.json

# Custom workflow directory name
./scripts/launch.sh your_input_config.json my_analysis_name

# With custom cache locations (recommended for HPC)
./scripts/launch.sh your_input_config.json my_analysis /scratch/user/cache /scratch/user/tmp

# With custom miniwdl.cfg
./scripts/launch.sh your_input_config.json my_analysis "" "" /path/to/custom/miniwdl.cfg

# Full customization
./scripts/launch.sh your_input_config.json my_analysis /scratch/cache /scratch/tmp /custom/miniwdl.cfg
```

### Usage

```bash
./scripts/launch.sh <input_config_json> [work_dir_name] [cache_dir] [tmp_dir] [miniwdl_cfg]
```

**Arguments:**
- `input_config_json`: JSON configuration file with sample information and file paths
- `work_dir_name`: Optional name for the working directory (default: `workflow_run_YYYYMMDD_HHMMSS`)
- `cache_dir`: Optional Apptainer cache directory (default: `<work_dir>/miniwdl_cache`)
- `tmp_dir`: Optional Apptainer temp directory (default: `<work_dir>/miniwdl_tmp`)
- `miniwdl_cfg`: Optional path to custom miniwdl.cfg (default: `<repo_root>/miniwdl.cfg`)

### What the Script Does

1. **Setup Validation**: Checks that required files, directories, and configurations exist
2. **Directory Creation**: Creates organized work directory structure:
   ```
   work_directory/
   ├── inputs/                    # Intelligently copied input files
   ├── outputs/                   # Workflow outputs (created during execution)
   ├── logs/                      # Execution logs
   ├── miniwdl.cfg               # MiniWDL configuration (from repo root or custom)
   ├── family.hpc.inputs.json    # Generated inputs file with local paths
   └── run_workflow.sh           # Executable run script (runs from repo root)
   ```
3. **Smart File Management**: Intelligently copies input files with duplicate detection
   - Skips files that already exist and are identical (MD5 hash comparison)
   - Provides clear feedback on copied vs skipped files
   - Saves time and storage space for repeated setup runs
4. **Cache Management**: Sets up Apptainer/Singularity cache structure
   - Creates cache subdirectories: `call_cache/`, `download_cache/`, `singularity_cache/`, `tmp/`
   - Supports custom cache locations (recommended for HPC environments)
   - Updates miniwdl.cfg with correct cache paths
5. **Configuration Management**: 
   - Uses miniwdl.cfg from repository root by default
   - Supports custom miniwdl.cfg paths
   - Automatically updates cache paths in configuration
6. **JSON Generation**: Creates new `family.hpc.inputs.json` with updated local file paths
7. **Run Script Creation**: Generates `run_workflow.sh` designed to be executed from repository root to avoid miniwdl I/O issues

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

1. **Conda Environment**: Ensure you have a conda environment with miniwdl installed
   ```bash
   conda create -n hifi-wgs miniwdl
   conda activate hifi-wgs
   ```

2. **Container Runtime**: Singularity/Apptainer must be available on your system

3. **Input Files**: Ensure all file paths in your input configuration JSON are accessible

4. **Reference Data**: Download the reference bundle from Zenodo if using local paths

### Advanced Usage Examples

**Multiple cache locations:**
```bash
# Use different locations for cache and temp
./scripts/launch.sh config.json analysis /fast/cache /local/tmp
```

**Custom miniwdl configuration:**
```bash
# Use your own miniwdl.cfg with specific cluster settings
./scripts/launch.sh config.json analysis "" "" /path/to/my/custom.cfg
```

**Re-running with same data:**
```bash
# Second run skips file copying (files are identical)
./scripts/launch.sh same_config.json analysis_v2
# Output: "Skipping file.bam -> inputs/file.bam (already exists and identical)"
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

The script supports flexible cache directory configuration, which is especially important for HPC environments:

**Default (local cache):**
```bash
./scripts/launch.sh config.json my_analysis
# Creates cache at: my_analysis/miniwdl_cache/
```

**Custom cache location (recommended for HPC):**
```bash
./scripts/launch.sh config.json my_analysis /scratch/user/cache /scratch/user/tmp
# Uses high-performance storage for cache operations
```

**Cache Structure Created:**
```
cache_location/
├── call_cache/           # miniwdl call caching
├── cromwell_db/          # Cromwell database  
├── download_cache/       # Apptainer image downloads (APPTAINER_CACHEDIR)
├── singularity_cache/    # Singularity/Apptainer images
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

# 2. Launch the workflow setup with custom cache location
./scripts/launch.sh my_trio_config.json trio_analysis_2024 /scratch/user/cache

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