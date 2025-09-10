#!/bin/bash
# launch.sh: Automated launcher for HiFi human WGS family workflow
# 
# This script automates the setup and execution of the HiFi human WGS workflow
# by copying input files to a local directory structure and generating the
# appropriate inputs JSON file.

# Usage:
#   ./launch.sh <input_config_json> [work_dir_name] [cache_dir] [tmp_dir] [miniwdl_cfg]
#
# Arguments:
#   input_config_json: JSON file with sample information and file paths (similar to family.hpc.inputs.v2.json)
#   work_dir_name: Optional name for working directory (default: workflow_run_$(date +%Y%m%d_%H%M%S))
#                  If provided without path separators, creates ./[work_dir_name]
#                  If provided with path, uses as full path
#   cache_dir: Optional Apptainer/Singularity cache directory (default: <work_dir>/miniwdl_cache)
#   tmp_dir: Optional Apptainer/Singularity temp directory (default: <work_dir>/miniwdl_tmp)
#   miniwdl_cfg: Optional path to custom miniwdl.cfg file (default: <repo_root>/miniwdl.cfg)

set -euo pipefail

# Check for required arguments
if [[ $# -lt 1 ]]; then
	# <<- here-document ignores leading tab, preserves spaces
    cat >&2 <<- EOF
	Usage: $0 <input_config_json> [work_dir_name] [cache_dir] [tmp_dir] [miniwdl_cfg]
	Arguments:
	  input_config_json: JSON file with sample information and file paths
	  work_dir_name: Optional name for working directory (default: workflow_run_YYYYMMDD_HHMMSS)
	                 Examples: 'my_analysis', 'family_trio_2024', '/full/path/to/workdir'
	  cache_dir: Optional Apptainer cache directory (default: <work_dir>/miniwdl_cache)
	  tmp_dir: Optional Apptainer temp directory (default: <work_dir>/miniwdl_tmp)
	  miniwdl_cfg: Optional path to custom miniwdl.cfg (default: <repo_root>/miniwdl.cfg)
	Example input_config_json format:
	{
	  "humanwgs_family.family": {
	    "family_id": "FAM1",
	    "samples": [
	      {
	        "sample_id": "sample1",
	        "hifi_reads": ["/path/to/reads.bam"],
	        "sex": "MALE",
	        "affected": false
	      }
	    ]
	  },
	  "humanwgs_family.ref_map_file": "/path/to/ref_map.tsv",
	  "humanwgs_family.somatic_map_file": "/path/to/somatic_map.tsv",
	  "humanwgs_family.tertiary_map_file": "/path/to/tertiary_map.tsv"
	}
	EOF
fi

# Get absolute paths before changing directories
INPUT_CONFIG_JSON="$(readlink -f "$1")"
WORK_DIR_NAME="${2:-workflow_run_$(date +%Y%m%d_%H%M%S)}"
CACHE_DIR="${3:-}"
TMP_DIR="${4:-}"
MINIWDL_CFG="${5:-}"

#prepull images
bash scripts/create_image_manifest.sh
bash scripts/populate_miniwdl_singularity_cache.sh image_manifest.txt $CACHE_DIR/singularity_cache


# Handle work directory path
if [[ "$WORK_DIR_NAME" == *"/"* ]]; then
    # Contains path separators, use as full path
    WORK_DIR="$WORK_DIR_NAME"
else
    # No path separators, create in current directory
    WORK_DIR="$(pwd)/$WORK_DIR_NAME"
fi

# Get absolute path for work directory
WORK_DIR="$(readlink -f "$WORK_DIR" 2>/dev/null || echo "$WORK_DIR")"

# Validate input file exists
if [[ ! -f "$INPUT_CONFIG_JSON" ]]; then
    echo "Error: Input config file '$1' does not exist." >&2
    exit 1
fi

# Get repository root path
SCRIPT_DIR="$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

# Validate repository structure
if [[ ! -d "$REPO_ROOT/workflows" ]]; then
    echo "Error: workflows directory not found at $REPO_ROOT/workflows" >&2
    exit 1
fi

# Determine miniwdl.cfg path
if [[ -n "$MINIWDL_CFG" ]]; then
    # User specified a custom config
    MINIWDL_CFG_PATH="$(readlink -f "$MINIWDL_CFG" 2>/dev/null || echo "$MINIWDL_CFG")"
    if [[ ! -f "$MINIWDL_CFG_PATH" ]]; then
        echo "Error: Custom miniwdl.cfg not found at $MINIWDL_CFG_PATH" >&2
        exit 1
    fi
else
    # Use default from repository root
    MINIWDL_CFG_PATH="$REPO_ROOT/miniwdl.cfg"
    if [[ ! -f "$MINIWDL_CFG_PATH" ]]; then
        echo "Error: Default miniwdl.cfg not found at $MINIWDL_CFG_PATH" >&2
        echo "Please specify a custom miniwdl.cfg path or ensure miniwdl.cfg exists in repository root" >&2
        exit 1
    fi
fi

# Set up cache directories
if [[ -n "$CACHE_DIR" ]]; then
    ACTUAL_CACHE_DIR="$(readlink -f "$CACHE_DIR" 2>/dev/null || echo "$CACHE_DIR")"
else
    ACTUAL_CACHE_DIR="$WORK_DIR/miniwdl_cache"
fi

if [[ -n "$TMP_DIR" ]]; then
    ACTUAL_TMP_DIR="$(readlink -f "$TMP_DIR" 2>/dev/null || echo "$TMP_DIR")"
else
    ACTUAL_TMP_DIR="$WORK_DIR/miniwdl_tmp"
fi

# Set symbolic link paths in work directory
CACHE_SYMLINK="$WORK_DIR/miniwdl_cache"
TMP_SYMLINK="$WORK_DIR/miniwdl_tmp"

echo "=== HiFi Human WGS Family Workflow Launcher ==="
echo "Input config: $INPUT_CONFIG_JSON"
echo "Work directory: $WORK_DIR"
echo "Repository root: $REPO_ROOT"
echo "MiniWDL config: $MINIWDL_CFG_PATH"
echo "Apptainer cache: $ACTUAL_CACHE_DIR"
echo "Apptainer temp: $ACTUAL_TMP_DIR"
echo ""

# Create work directory structure
echo "Creating work directory structure..."
mkdir -p "$WORK_DIR"/{inputs,outputs,logs}
cd "$WORK_DIR"

# Copy the miniwdl.cfg file and update cache paths
echo "Copying and updating configuration files..."
cp "$MINIWDL_CFG_PATH" miniwdl.cfg

# Update cache paths in miniwdl.cfg to point to actual cache locations
echo "Updating cache paths in miniwdl.cfg..."
sed -i.bak \
    -e "s|dir = \"[^\"]*download_cache[^\"]*\"|dir = \"$ACTUAL_CACHE_DIR/download_cache\"|g" \
    -e "s|dir = \"[^\"]*call_cache[^\"]*\"|dir = \"$ACTUAL_CACHE_DIR/call_cache\"|g" \
    -e "s|image_cache = \"[^\"]*singularity_cache[^\"]*\"|image_cache = \"$ACTUAL_CACHE_DIR/singularity_cache\"|g" \
    miniwdl.cfg

# Remove backup file
rm -f miniwdl.cfg.bak

# Parse input JSON and copy files
echo "Processing input files..."

# Create a temporary Python script to process the JSON
cat > process_config.py << 'EOF'
import json
import os
import shutil
import sys
import filecmp
from pathlib import Path

def should_copy_file(src_path, dest_path):
    """Check if file should be copied (doesn't exist or is different)"""
    if not dest_path.exists():
        return True

    return not filecmp.cmp(src_path, dest_path, shallow=True)

if len(sys.argv) != 2:
    print("Usage: python3 process_config.py <input_config_json>")
    sys.exit(1)

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
                
                if should_copy_file(hifi_read_path, dest_path):
                    print(f"Copying {hifi_read_path} -> {dest_path}")
                    shutil.copy2(hifi_read_path, dest_path)
                else:
                    print(f"Skipping {hifi_read_path} -> {dest_path} (already exists and identical)")
                
                new_hifi_reads.append(str(dest_path.absolute()))
            else:
                print(f"Warning: HiFi reads file not found: {hifi_read_path}")
                new_hifi_reads.append(hifi_read_path)  # Keep original path
        
        new_sample["hifi_reads"] = new_hifi_reads
        new_family["samples"].append(new_sample)
    
    new_config["humanwgs_family.family"] = new_family

# Copy reference files if they exist locally
ref_files = ["ref_map_file", "somatic_map_file", "tertiary_map_file"]
for ref_key in ref_files:
    full_key = f"humanwgs_family.{ref_key}"
    if full_key in config:
        ref_path = config[full_key]
        if os.path.exists(ref_path):
            filename = f"{ref_key}.tsv"
            dest_path = inputs_dir / filename
            
            if should_copy_file(ref_path, dest_path):
                print(f"Copying {ref_path} -> {dest_path}")
                shutil.copy2(ref_path, dest_path)
            else:
                print(f"Skipping {ref_path} -> {dest_path} (already exists and identical)")
            
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
EOF

# Run the Python script
python3 process_config.py "$INPUT_CONFIG_JSON"

# Clean up the temporary script
rm process_config.py

# Set up Singularity/Apptainer cache directories
echo "Setting up container cache..."

# Create the actual cache directory structure
mkdir -p "$ACTUAL_CACHE_DIR"/{call_cache,cromwell_db,download_cache,singularity_cache,tmp}
mkdir -p "$ACTUAL_TMP_DIR"

# Create run script
echo "Creating run script..."

cat > run_workflow.sh << EOF
#!/bin/bash
set -euo pipefail

# NOTE: This script is designed to be run from the repository root directory
# Usage: bash $WORK_DIR/run_workflow.sh
# Make sure to activate your conda environment with miniwdl before running this script
# Example: conda activate hifi-wgs

# Check if we're in the repository root
if [[ ! -d "workflows" ]] || [[ ! -f "workflows/family.wdl" ]]; then
    echo "Error: This script must be run from the repository root directory" >&2
    echo "Current directory: \$(pwd)" >&2
    echo "Expected to find: workflows/family.wdl" >&2
    exit 1
fi

# Set container cache directories
export APPTAINER_CACHEDIR="$ACTUAL_CACHE_DIR/download_cache"
export APPTAINER_TMP_DIR="$ACTUAL_CACHE_DIR/tmp"

echo "Starting HiFi human WGS family workflow..."
echo "Timestamp: \$(date)"
echo "Repository root: \$(pwd)"
echo "Work directory: $WORK_DIR"
echo "Cache directory: \$APPTAINER_CACHEDIR"
echo ""

miniwdl run \\
    --cfg $WORK_DIR/miniwdl.cfg \\
    --dir $WORK_DIR/outputs \\
    workflows/family.wdl \\
    -i $WORK_DIR/family.hpc.inputs.json \\
    --verbose 2>&1 | tee $WORK_DIR/logs/workflow_\$(date +%Y%m%d_%H%M%S).log

echo ""
echo "Workflow completed at: \$(date)"
echo "Check outputs in: $WORK_DIR/outputs"
EOF

chmod +x run_workflow.sh

# Summary
echo ""
echo "=== Setup Complete ==="
echo "Work directory: $WORK_DIR"
echo "Input files copied to: $WORK_DIR/inputs/"
echo "Generated config: $WORK_DIR/family.hpc.inputs.json"
echo "Using workflow from: $REPO_ROOT/workflows/"
echo "MiniWDL config: $WORK_DIR/miniwdl.cfg"
echo "Apptainer cache: $ACTUAL_CACHE_DIR"
echo "Run script: $WORK_DIR/run_workflow.sh"
echo ""
echo "To run the workflow:"
echo "  1. Activate your conda environment: conda activate hifi-wgs"
echo "  2. cd $REPO_ROOT"
echo "  3. bash $WORK_DIR/run_workflow.sh"
echo ""
echo "Or run directly from repository root:"
echo "  conda activate hifi-wgs"
echo "  cd $REPO_ROOT"
echo "  miniwdl run --cfg $WORK_DIR/miniwdl.cfg --dir $WORK_DIR/outputs workflows/family.wdl -i $WORK_DIR/family.hpc.inputs.json --verbose"
echo ""
echo "Monitor progress with:"
echo "  tail -f $WORK_DIR/logs/workflow_*.log"

