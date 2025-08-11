#!/bin/bash
# launch.sh: helper script to run a WGS-human-edit pipeline on a single dataset via Slurm

# Usage:
#   ./launch.sh <hifi_reads_path> <ccs_stats_json> [parent_reads_path] [sample_name]

set -euo pipefail

# Input arguments
HIFI_READS_PATH="$1"
CCS_STATS_JSON="$2"
# Parent reads optional
PARENT_READS_PATH="${3:-}"
# Sample name optional
SAMPLE_NAME="${4:-sample}"
ENV_YML="${5:-environment.yml}"

# --- Conda/Mamba detection ---
if [[ -f "$HOME/micromamba/bin/micromamba" ]]; then
    eval "$("$HOME/micromamba/bin/micromamba" shell hook)"
    CONDA_CMD="micromamba"
elif command -v mamba &>/dev/null; then
    eval "$(mamba shell hook)"
    CONDA_CMD="mamba"
elif command -v conda &>/dev/null; then
    eval "$(conda shell hook)"
    CONDA_CMD="conda"
elif module avail anaconda &>/dev/null 2>&1; then
    module load anaconda
    eval "$(conda shell hook)"
    CONDA_CMD="conda"
else
    echo "Error: micromamba or conda not found. Please install micromamba or Anaconda." >&2
    exit 1
fi

# --- Environment activation ---
if [[ ! -f "$ENV_YML" ]]; then
    echo "Error: environment file '$ENV_YML' does not exist." >&2
    exit 1
fi
ENV_NAME=$(grep -E '^name:' "$ENV_YML" | head -1 | awk '{print $2}')
if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
    echo "Activating existing environment: $ENV_NAME"
    "$CONDA_CMD" activate "$ENV_NAME"
else
    echo "Error: Conda environment '$ENV_NAME' not found."
    echo "Please create it with: $CONDA_CMD env create -f $ENV_YML" >&2
    exit 1
fi

# Directories
WORKDIR="$(pwd)/work/${SAMPLE_NAME}"
RES_DIR="$(pwd)/results/${SAMPLE_NAME}"
INPUT_JSON="${WORKDIR}/inputs.json"
MINIWDL_CFG="${WORKDIR}/miniwdl.cfg"
OUTPUTS_DIR="${WORKDIR}/outputs"
PIPELINE="${WORKDIR}/workflows/family.wdl"

# Set Singularity cache to high-capacity storage
export APPTAINER_CACHEDIR="${WORKDIR}/miniwdl_cache/singularity_cache/"
export APPTAINER_TMP_DIR="${WORKDIR}/miniwdl_cache/tmp"
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMP_DIR"

# Pipeline command using MiniWDL
PIPELINE_CMD="miniwdl run --cfg ${MINIWDL_CFG} ${PIPELINE} -i ${INPUT_JSON} --verbose"

# Prepare Slurm submission 
SBATCH_SCRIPT="${WORKDIR}/slurm_${SAMPLE_NAME}.sh"
mkdir -p "${WORKDIR}" "${RES_DIR}"

cat > "${SBATCH_SCRIPT}" <<EOF
#!/bin/bash

#NEED TO ADD CONDA ENV HERE 


# Copy inputs
mkdir -p ${WORKDIR}/reads
cp ${HIFI_READS_PATH} ${WORKDIR}/reads/hifi_reads.bam
cp ${CCS_STATS_JSON} ${WORKDIR}/reads/ccs_stats.json
EOF

# Conditionally copy parent reads
if [[ -n "${PARENT_READS_PATH}" ]]; then
  cat >> "${SBATCH_SCRIPT}" <<EOF
cp ${PARENT_READS_PATH} ${WORKDIR}/reads/parent.bam
EOF
fi

# Continue SBATCH script: JSON generation, pipeline, cleanup
cat >> "${SBATCH_SCRIPT}" <<'EOF'

# Generate inputs JSON
{
  echo "{"
  echo "  \"pipeline.hifi_reads\": \"${WORKDIR}/reads/hifi_reads.bam\",
  echo "  \"pipeline.ccs_stats\": \"${WORKDIR}/reads/ccs_stats.json\",
EOF

if [[ -n "${PARENT_READS_PATH}" ]]; then
  cat >> "${SBATCH_SCRIPT}" <<'EOF'
  echo "  \"pipeline.parent_reads\": \"${WORKDIR}/reads/parent.bam\",
EOF
fi
cat >> "${SBATCH_SCRIPT}" <<'EOF'
  echo "  \"pipeline.sample_name\": \"${SAMPLE_NAME}\"
  echo "}"
} > ${INPUT_JSON}

# Run pipeline with MiniWDL
${PIPELINE_CMD}

# Move outputs and config
mv ${WORKDIR}/*.{bam,vcf,json,inputs.json} ${RES_DIR}/ || true

# Cleanup intermediate reads
rm -rf ${WORKDIR}/reads
EOF

# Submit job
touch ${WORKDIR}/slurm_submit_marker
sbatch "${SBATCH_SCRIPT}"
echo "Submitted Slurm job for sample: ${SAMPLE_NAME}. Check logs in ${WORKDIR}."

