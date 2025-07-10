#!/bin/bash

export APPTAINER_CACHEDIR="/hpc/mydata/ram.ayyala/HiFi-human-WGS-editing-QC-WDL/miniwdl_cache/singularity_cache/"
export APPTAINER_TMP_DIR=="/hpc/mydata/ram.ayyala/HiFi-human-WGS-editing-QC-WDL/miniwdl_cache/tmp"
mkdir -p "${APPTAINER_CACHEDIR}" "${APPTAINER_TMP_DIR}"


manifest="../image_manifest.txt"
singularity_cache="${APPTAINER_CACHEDIR}"

while read -r container_url; do
    # Skip empty lines and comments
    [[ -z "$container_url" || "$container_url" == \#* ]] && continue
       if [[ "${container_url}" == *@sha256:* ]]; then
        # SHA version
        repo_tag="${container_url%%@*}"
        sha="${container_url#*@sha256:}"  
        uri="docker://${repo_tag}@sha256:${sha}"
        sif_filename="$(echo "${repo_tag}" | tr '/' '-')@sha256-${sha:0:12}.sif"
    elif [[ "${container_url}" == *:* && ! "${container_url}" == *@* ]]; then  
        # Version tag format
        uri="docker://${container_url}"
        sif_filename="$(echo "${container_url}" | tr ':' '-' | tr '/' '-').sif"

    else
        echo "The following line is not a valid container uri: ${container_url}"
        continue
    fi

    #check pre-pull has already occured for this tool 
    if [[ -f "${singularity_cache}/${sif_filename}" ]]; then
        echo "Skipping ${container_url} as it already exists: ${sif_filename}"
        continue
    fi
    
    echo "Pulling ${container_url}"
    if ! singularity pull --dir "${singularity_cache}" "${uri}"; then
        echo "ERROR: Failed to pull ${container_url}"
        exit 1
    fi

done < "${manifest}"
echo "Images pulled to: ${singularity_cache}"

