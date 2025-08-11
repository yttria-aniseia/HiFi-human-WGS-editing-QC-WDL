#!/bin/bash
# launch.sh: Automated launcher for HiFi human WGS family workflow
# 
# This script automates the setup and execution of the HiFi human WGS workflow
# by copying input files to a local directory structure and generating the 
# appropriate inputs JSON file.

# Usage:
#   ./launch.sh <input_config_json> [work_dir_name] [conda_env]
#
# Arguments:
#   input_config_json: JSON file with sample information and file paths (similar to family.hpc.inputs.v2.json)
#   work_dir_name: Optional name for working directory (default: workflow_run_$(date +%Y%m%d_%H%M%S))
#                  If provided without path separators, creates ./[work_dir_name]
#                  If provided with path, uses as full path
#   conda_env: Optional conda environment name (default: hifi-wgs)

set -euo pipefail

# Check for required arguments
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <input_config_json> [work_dir_name] [conda_env]" >&2
    echo "" >&2
    echo "Arguments:" >&2
    echo "  input_config_json: JSON file with sample information and file paths" >&2
    echo "  work_dir_name: Optional name for working directory (default: workflow_run_YYYYMMDD_HHMMSS)" >&2
    echo "                 Examples: 'my_analysis', 'family_trio_2024', '/full/path/to/workdir'" >&2
    echo "  conda_env: Optional conda environment name (default: hifi-wgs)" >&2
    echo "" >&2
    echo "Example input_config_json format:" >&2
    echo '{' >&2
    echo '  "humanwgs_family.family": {' >&2
    echo '    "family_id": "FAM1",' >&2
    echo '    "samples": [' >&2
    echo '      {' >&2
    echo '        "sample_id": "sample1",' >&2
    echo '        "hifi_reads": ["/path/to/reads.bam"],' >&2
    echo '        "sex": "MALE",' >&2
    echo '        "affected": false' >&2
    echo '      }' >&2
    echo '    ]' >&2
    echo '  },' >&2
    echo '  "humanwgs_family.ref_map_file": "/path/to/ref_map.tsv",' >&2
    echo '  "humanwgs_family.tertiary_map_file": "/path/to/tertiary_map.tsv"' >&2
    echo '}' >&2
    exit 1
fi

# Input arguments
INPUT_CONFIG_JSON="$1"
WORK_DIR_NAME="${2:-workflow_run_$(date +%Y%m%d_%H%M%S)}"
CONDA_ENV="${3:-hifi-wgs}"

# Handle work directory path
if [[ "$WORK_DIR_NAME" == *"/"* ]]; then
    # Contains path separators, use as full path
    WORK_DIR="$WORK_DIR_NAME"
else
    # No path separators, create in current directory
    WORK_DIR="./$WORK_DIR_NAME"
fi

# Validate input file exists
if [[ ! -f "$INPUT_CONFIG_JSON" ]]; then
    echo "Error: Input config file '$INPUT_CONFIG_JSON' does not exist." >&2
    exit 1
fi

echo "=== HiFi Human WGS Family Workflow Launcher ==="
echo "Input config: $INPUT_CONFIG_JSON"
echo "Work directory: $WORK_DIR"
echo "Conda environment: $CONDA_ENV"
echo ""

# --- Conda/Mamba detection and activation ---
echo "Setting up conda environment..."
if [[ -f "$HOME/micromamba/bin/micromamba" ]]; then
    eval "$("$HOME/micromamba/bin/micromamba" shell hook)"
    CONDA_CMD="micromamba"
elif command -v mamba &>/dev/null; then
    eval "$(mamba shell hook)"
    CONDA_CMD="mamba"
elif command -v conda &>/dev/null; then
    eval "$(conda shell hook)"
    CONDA_CMD="conda"
else
    echo "Error: conda/mamba not found. Please install conda or mamba." >&2
    exit 1
fi

# Activate conda environment
if conda env list | awk '{print $1}' | grep -qx "$CONDA_ENV"; then
    echo "Activating conda environment: $CONDA_ENV"
    conda activate "$CONDA_ENV"
else
    echo "Error: Conda environment '$CONDA_ENV' not found." >&2
    echo "Available environments:" >&2
    conda env list >&2
    exit 1
fi

# Create work directory structure
echo "Creating work directory structure..."
mkdir -p "$WORK_DIR"/{inputs,outputs,logs}
cd "$WORK_DIR"

# Copy workflow files to work directory
echo "Copying workflow files..."
SCRIPT_DIR="$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

cp -r "$REPO_ROOT/workflows" .
cp -r "$REPO_ROOT/backends" .
cp "$REPO_ROOT/backends/hpc/miniwdl.cfg" .

# Parse input JSON and copy files
echo "Processing input files..."
python3 << 'EOF'
import json
import os
import shutil
import sys
from pathlib import Path

# Read input config
with open(sys.argv[1], 'r') as f:
    config = json.load(f)

# Create new config with local paths
new_config = {}
inputs_dir = Path("inputs")

# Process family samples
if "humanwgs_family.family" in config:
    family = config["humanwgs_family.family"]
    new_family = {
        "family_id": family["family_id"],
        "samples": []
    }
    
    for sample in family["samples"]:
        new_sample = sample.copy()
        new_hifi_reads = []
        
        # Copy HiFi reads files
        for i, hifi_read_path in enumerate(sample["hifi_reads"]):
            if os.path.exists(hifi_read_path):
                filename = f"{sample['sample_id']}_hifi_reads_{i}.bam"
                dest_path = inputs_dir / filename
                print(f"Copying {hifi_read_path} -> {dest_path}")
                shutil.copy2(hifi_read_path, dest_path)
                new_hifi_reads.append(str(dest_path.absolute()))
            else:
                print(f"Warning: HiFi reads file not found: {hifi_read_path}")
                new_hifi_reads.append(hifi_read_path)  # Keep original path
        
        new_sample["hifi_reads"] = new_hifi_reads
        new_family["samples"].append(new_sample)
    
    new_config["humanwgs_family.family"] = new_family

# Copy reference files if they exist locally
ref_files = ["ref_map_file", "tertiary_map_file"]
for ref_key in ref_files:
    full_key = f"humanwgs_family.{ref_key}"
    if full_key in config:
        ref_path = config[full_key]
        if os.path.exists(ref_path):
            filename = f"{ref_key}.tsv"
            dest_path = inputs_dir / filename
            print(f"Copying {ref_path} -> {dest_path}")
            shutil.copy2(ref_path, dest_path)
            new_config[full_key] = str(dest_path.absolute())
        else:
            print(f"Warning: Reference file not found: {ref_path}")
            new_config[full_key] = ref_path  # Keep original path

# Copy other config values
for key, value in config.items():
    if key not in new_config:
        new_config[key] = value

# Write new inputs JSON
with open("family.hpc.inputs.json", "w") as f:
    json.dump(new_config, f, indent=2)

print("Generated family.hpc.inputs.json with local file paths")
EOF python3 "$INPUT_CONFIG_JSON"

# Set up Singularity/Apptainer cache
echo "Setting up container cache..."
export APPTAINER_CACHEDIR="$(pwd)/miniwdl_cache/singularity_cache"
export APPTAINER_TMP_DIR="$(pwd)/miniwdl_cache/tmp"
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMP_DIR"

# Create run script
echo "Creating run script..."
cat > run_workflow.sh << 'EOF'
#!/bin/bash
set -euo pipefail

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate hifi-wgs

# Set container cache
export APPTAINER_CACHEDIR="$(pwd)/miniwdl_cache/singularity_cache"
export APPTAINER_TMP_DIR="$(pwd)/miniwdl_cache/tmp"

# Run workflow
echo "Starting HiFi human WGS family workflow..."
echo "Timestamp: $(date)"
echo "Working directory: $(pwd)"
echo ""

miniwdl run \
    --cfg miniwdl.cfg \
    --dir outputs \
    workflows/family.wdl \
    -i family.hpc.inputs.json \
    --verbose 2>&1 | tee logs/workflow_$(date +%Y%m%d_%H%M%S).log

echo ""
echo "Workflow completed at: $(date)"
echo "Check outputs in: $(pwd)/outputs"
EOF

chmod +x run_workflow.sh

# Summary
echo ""
echo "=== Setup Complete ==="
echo "Work directory: $(pwd)"
echo "Input files copied to: $(pwd)/inputs/"
echo "Generated config: $(pwd)/family.hpc.inputs.json"
echo "Workflow files: $(pwd)/workflows/"
echo ""
echo "To run the workflow:"
echo "  cd $(pwd)"
echo "  ./run_workflow.sh"
echo ""
echo "Or run directly with:"
echo "  miniwdl run --cfg miniwdl.cfg --dir outputs workflows/family.wdl -i family.hpc.inputs.json --verbose"
echo ""
echo "Monitor progress with:"
echo "  tail -f logs/workflow_*.log"

