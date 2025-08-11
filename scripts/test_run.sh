#!/bin/bash

eval "$(/home/ram.ayyala/miniforge3/condabin/conda shell.bash hook)"
conda activate hifi-wgs  # Your Conda env with miniwdl

path="/hpc/scratch/group.data.science/ram_temp/HiFi-human-WGS-editing-QC-WDL"
# inputs="${path}/backends/hpc/family.hpc.inputs.KOLF2.1J.json"
output="${path}/wdl-outputs"
mkdir -p "${output}"



# Set Singularity cache to high-capacity storage
export APPTAINER_CACHEDIR="/hpc/scratch/group.data.science/ram_temp/HiFi-human-WGS-editing-QC-WDL/miniwdl_cache/singularity_cache/"
export APPTAINER_TMP_DIR="/hpc/scratch/group.data.science/ram_temp/HiFi-human-WGS-editing-QC-WDL/miniwdl_cache/tmp"
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMP_DIR"

miniwdl run \
    --cfg ${path}/miniwdl.cfg \
    --dir wdl-outputs \
    ${path}/workflows/family.wdl \
    -i ${path}/backends/hpc/family.hpc.inputs_v2.json \
    --verbose 






