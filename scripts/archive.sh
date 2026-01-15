#!/usr/bin/env bash

set -euo pipefail

DEFAULT_ARCHIVE_BASEDIR=""

DELETE_AFTER=false
INTERACTIVE=false
while getopts "di" opt; do
    case "$opt" in
        d) DELETE_AFTER=true ;;
        i) INTERACTIVE=true ;;
        *)
            echo "Usage: $0 [-d] [-i] <run_directory> [archive_base_dir]" >&2
            exit 1
            ;;
    esac
done
shift $((OPTIND-1))

if [[ $# -lt 1 ]]; then
    echo "Error: Run directory not specified" >&2
    echo "Usage: $0 [-d] [-i] <run_directory> [archive_base_dir]" >&2
    exit 1
fi

RUNNAME="${1}"
ARCHIVE_BASEDIR="${2:-$DEFAULT_ARCHIVE_BASEDIR}"
SRC="${RUNNAME}/outputs/_LAST"
DEST="${ARCHIVE_BASEDIR}/${RUNNAME}"

if [[ ! -d "${SRC}" ]]; then
    echo "Error: Source directory does not exist: ${SRC}" >&2
    exit 1
fi

echo "Source:      $(realpath "${SRC}")"
echo "Destination: ${DEST}"
echo ""

if [[ "${INTERACTIVE}" == true ]]; then
    read -p "Proceed with archive? [y/N] " -n 1 -r
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Aborted"
        exit 0
    fi
fi

mkdir -p "${DEST}"

for file in inputs.json outputs.json error.json workflow.log; do
    [[ -f "${SRC}/${file}" ]] && cp "${SRC}/${file}" "${DEST}/"
done

# Extract and copy expected_edit JSON files
if [[ -f "${SRC}/inputs.json" ]]; then
    while IFS= read -r edit_file; do
        if [[ -n "${edit_file}" && -f "${edit_file}" ]]; then
            cp "${edit_file}" "${DEST}/$(basename "${edit_file}")"
        fi
    done < <(grep -Po '"expected_edit"\s*:\s*"\K[^"]+' "${SRC}/inputs.json" 2>/dev/null || true)
fi

# Backup original and rewrite with relative paths
if [[ -f "${DEST}/inputs.json" ]]; then
    cp "${DEST}/inputs.json" "${DEST}/inputs.orig.json"
    sed 's|"expected_edit"\s*:\s*"[^"]*\/\([^/]*\.json\)"|"expected_edit": "./\1"|g' \
        "${DEST}/inputs.orig.json" > "${DEST}/inputs.json"
fi

rsync -arRPK -L --copy-unsafe-links "${SRC}/wdl" "${ARCHIVE_BASEDIR}/"
rsync -arRPK -L --copy-unsafe-links "${SRC}/out" "${ARCHIVE_BASEDIR}/"

if [[ "${DELETE_AFTER}" == true ]]; then
    if [[ "${INTERACTIVE}" == true ]]; then
        read -p "Delete source outputs? [y/N] " -n 1 -r
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            echo "Skipped deletion"
        else
            rm -rf "${RUNNAME}/outputs/"
            echo "Source outputs deleted"
        fi
    else
        rm -rf "${RUNNAME}/outputs/"
        echo "Source outputs deleted"
    fi
fi

echo "Archive complete: ${DEST}"
