#!/bin/bash
# launch.sh: Automated launcher for HiFi human WGS family workflow
#
# This script automates the setup and execution of the HiFi human WGS workflow
# in multiple phases:
#   Phase 1: Setup scripts in current shell (prepull images, create work dir)
#   Phase 2: Input processing in SLURM job (copy/merge/strip files)
#   Phase 3: Optional workflow execution (if --run flag provided)
#
# TODO: Phase 4: Automatic report generation after workflow completion

# Usage:
#   ./launch.sh <input_config_json> [options]
#
# Arguments:
#   input_config_json: JSON file with sample information and file paths
#
# Options:
#   --work-dir <path>    Work directory path (default: workflow_run_YYYYMMDD_HHMMSS)
#   --cache-dir <path>   Apptainer cache directory (default: <work_dir>/miniwdl_cache)
#   --tmp-dir <path>     Apptainer temp directory (default: <work_dir>/miniwdl_tmp)
#   --miniwdl-cfg <path> Custom miniwdl.cfg file (default: <repo_root>/miniwdl.cfg)
#   --run                Start miniwdl workflow after input processing
#   --help               Show this help message

set -euo pipefail

# Function to show help message
show_help() {
    cat <<- EOF
	Usage: $0 <input_config_json> [options]

	Arguments:
	  input_config_json: JSON file with sample information and file paths

	Options:
	  --work-dir <path>    Work directory path (default: workflow_run_YYYYMMDD_HHMMSS)
	  --cache-dir <path>   Apptainer cache directory (default: <work_dir>/miniwdl_cache)
	  --tmp-dir <path>     Apptainer temp directory (default: <work_dir>/miniwdl_tmp)
	  --miniwdl-cfg <path> Custom miniwdl.cfg file (default: <repo_root>/miniwdl.cfg)
	  --run                Start miniwdl workflow after input processing
	  --help               Show this help message

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
	  "humanwgs_family.tertiary_map_file": "/path/to/tertiary_map.tsv",
	  "humanwgs_family.expected_edits": "/path/to/expected_edits.json"
	}
	EOF
}

# Parse arguments
INPUT_CONFIG_JSON=""
WORK_DIR_NAME="workflow_run_$(date +%Y%m%d_%H%M%S)"
CACHE_DIR=""
TMP_DIR=""
MINIWDL_CFG=""
RUN_WORKFLOW=false

if [[ $# -lt 1 ]]; then
    show_help >&2
    exit 1
fi

# First positional argument is input config
INPUT_CONFIG_JSON="$1"
shift

# Parse optional flags
while [[ $# -gt 0 ]]; do
    case $1 in
        --work-dir)
            WORK_DIR_NAME="$2"
            shift 2
            ;;
        --cache-dir)
            CACHE_DIR="$2"
            shift 2
            ;;
        --tmp-dir)
            TMP_DIR="$2"
            shift 2
            ;;
        --miniwdl-cfg)
            MINIWDL_CFG="$2"
            shift 2
            ;;
        --run)
            RUN_WORKFLOW=true
            shift
            ;;
        --help)
            show_help
            exit 0
            ;;
        *)
            echo "Error: Unknown option: $1" >&2
            show_help >&2
            exit 1
            ;;
    esac
done

# Get absolute path for input config
INPUT_CONFIG_JSON="$(readlink -f "$INPUT_CONFIG_JSON")"

# Validate input file exists
if [[ ! -f "$INPUT_CONFIG_JSON" ]]; then
    echo "Error: Input config file '$INPUT_CONFIG_JSON' does not exist." >&2
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

# Handle work directory path
if [[ "$WORK_DIR_NAME" == *"/"* ]]; then
    # Contains path separators, use as full path
    WORK_DIR="$WORK_DIR_NAME"
else
    # No path separators, create in current directory
    WORK_DIR="$(pwd)/$WORK_DIR_NAME"
fi

# Get absolute path for work directory (create if needed for readlink to work)
mkdir -p "$WORK_DIR"
WORK_DIR="$(readlink -f "$WORK_DIR")"

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

echo "=== HiFi Human WGS Family Workflow Launcher ==="
echo "Input config: $INPUT_CONFIG_JSON"
echo "Work directory: $WORK_DIR"
echo "Repository root: $REPO_ROOT"
echo "MiniWDL config: $MINIWDL_CFG_PATH"
echo "Apptainer cache: $ACTUAL_CACHE_DIR"
echo "Apptainer temp: $ACTUAL_TMP_DIR"
echo "Run workflow: $RUN_WORKFLOW"
echo ""

################################################################################
# PHASE 1: Setup in current shell
################################################################################
echo "=== Phase 1: Setup ==="

# Create cache directory structure first (needed for prepull)
echo "Creating cache directory structure..."
mkdir -p "$ACTUAL_CACHE_DIR"/{call_cache,cromwell_db,download_cache,singularity_cache,tmp}
mkdir -p "$ACTUAL_TMP_DIR"

# Export environment variables for Apptainer (needed for prepull)
export APPTAINER_CACHEDIR="$ACTUAL_CACHE_DIR/download_cache"
export APPTAINER_TMPDIR="$ACTUAL_CACHE_DIR/tmp"
export SINGULARITY_CACHEDIR="$ACTUAL_CACHE_DIR/download_cache"
export SINGULARITY_TMPDIR="$ACTUAL_CACHE_DIR/tmp"

# Prepull container images
echo "Creating image manifest..."
cd "$REPO_ROOT"
bash scripts/create_image_manifest.sh

echo "Prepulling container images to cache..."
bash scripts/populate_miniwdl_singularity_cache.sh image_manifest.txt "$ACTUAL_CACHE_DIR/singularity_cache"

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

echo "Phase 1 setup complete."
echo ""

################################################################################
# PHASE 2: Input processing (SLURM job)
################################################################################
echo "=== Phase 2: Input Processing ==="

# Copy Python script for input processing
cp "$REPO_ROOT/scripts/process_input_config.py" "$WORK_DIR/process_input_config.py"
chmod +x "$WORK_DIR/process_input_config.py"

# Submit the SLURM job
echo "Submitting input processing job..."
srun -u --cpus-per-task=16 --mem=16G --time=4:00:00 python3 process_input_config.py "$INPUT_CONFIG_JSON"

echo ""
echo "Phase 2 input processing complete."
echo ""

################################################################################
# PHASE 3: Optional workflow execution
################################################################################
echo "=== Phase 3: Workflow Execution ==="

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
export APPTAINER_TMPDIR="$ACTUAL_CACHE_DIR/tmp"
export SINGULARITY_CACHEDIR="$ACTUAL_CACHE_DIR/download_cache"
export SINGULARITY_TMPDIR="$ACTUAL_CACHE_DIR/tmp"

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
    --verbose

echo ""
echo "Workflow completed at: \$(date)"
echo "Check outputs in: $WORK_DIR/outputs"
EOF

chmod +x run_workflow.sh

# Optionally start the workflow
if [[ "$RUN_WORKFLOW" == "true" ]]; then
    echo "Starting workflow execution (--run flag provided)..."
    echo "NOTE: Make sure your conda environment with miniwdl is activated"
    echo ""
    cd "$REPO_ROOT"
    bash "$WORK_DIR/run_workflow.sh"
else
    echo "Workflow execution skipped (use --run flag to auto-start)"
    echo ""
fi

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

if [[ "$RUN_WORKFLOW" != "true" ]]; then
    echo "To run the workflow manually:"
    echo "  1. Activate your conda environment: conda activate hifi-wgs"
    echo "  2. cd $REPO_ROOT"
    echo "  3. bash $WORK_DIR/run_workflow.sh"
    echo ""
    echo "Or run directly from repository root:"
    echo "  conda activate hifi-wgs"
    echo "  cd $REPO_ROOT"
    echo "  miniwdl run --cfg $WORK_DIR/miniwdl.cfg --dir $WORK_DIR/outputs workflows/family.wdl -i $WORK_DIR/family.hpc.inputs.json --verbose"
    echo ""
fi

echo "Monitor progress with:"
echo "  tail -f $WORK_DIR/logs/workflow_*.log"

