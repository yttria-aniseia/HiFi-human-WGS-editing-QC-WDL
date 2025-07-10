#!/bin/bash
#SBATCH --job-name=large_input_job
#SBATCH --output=logs/slurm_family_%j.out
#SBATCH --error=logs/slurm_family_%j.err
#SBATCH --nodes=1                    # Single node (multi-node rarely needed for genomics)
#SBATCH --ntasks=1                   # One task (adjust if workflow parallelizes internally)
#SBATCH --cpus-per-task=64          # Match your workflow's thread count (e.g., pbmm2/pbsv)
#SBATCH --mem=700G                   # 242GB input + 20% buffer (300G)
#SBATCH --time=48:00:00              # Adjust based on tool runtime (start with 24h)
#SBATCH --tmp=500G                   # Temporary storage for intermediate files
# Replace with your email
#SBATCH --mail-user=ram.ayyala@czbiohub.org
#SBATCH --mail-type=ALL

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






