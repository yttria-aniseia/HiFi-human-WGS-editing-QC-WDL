#!/bin/bash
eval "$(/home/ram.ayyala/miniforge3/condabin/conda shell.bash hook)"
conda activate cromwell

path="/hpc/scratch/group.data.science/ram_temp/HiFi-human-WGS-editing-QC-WDL"
output="${path}/wdl-outputs"
mkdir -p "${output}"

# Singularity cache (correct names)
export SINGULARITY_CACHEDIR="${path}/miniwdl_cache/singularity_cache"
export SINGULARITY_TMPDIR="${path}/miniwdl_cache/tmp"
mkdir -p "$SINGULARITY_CACHEDIR" "$SINGULARITY_TMPDIR"

# Tell the conda cromwell CLI which config to use
export CROMWELL_CONFIG="${path}/cromwell.conf"

# Run the workflow
cromwell run \
  -v 1.0 \
  "${path}/workflows/family_cromwell.wdl" \
  --inputs "${path}/backends/hpc/family.hpc.inputs_v2.json" \
  --options workflow-options.json \
  --metadata-output "${output}/metadata.json" \
