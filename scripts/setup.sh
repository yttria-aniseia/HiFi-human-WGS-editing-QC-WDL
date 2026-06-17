#!/bin/bash
# setup.sh — one-time setup entry point.
#
#   1. fetches the FROZEN resource bundle from Zenodo (resumable, checksummed)
#      via scripts/fetch_resources.sh, and
#   2. builds the knock-knock Singularity image (no usable prebuilt container).
#
# Pass-through args (e.g. --dest DIR, --dbnsfp FILE) go to fetch_resources.sh.
#
# Maintainers rebuilding the bundle for a new resource version want
# scripts/build_resource_bundle.sh instead, not this script.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# --- 1. resources ------------------------------------------------------------
"${SCRIPT_DIR}/fetch_resources.sh" "$@"

# --- 2. knock-knock container ------------------------------------------------
if [[ -z "${SINGULARITY_CACHEDIR:-}" ]]; then
	echo "WARNING: SINGULARITY_CACHEDIR is not set; skipping knock-knock container build" >&2
else
	KK_SIF="${SINGULARITY_CACHEDIR}/docker___knock-knock.sif"
	if [[ ! -f "${KK_SIF}" ]]; then
		echo ">>> Building knock-knock Singularity image..."
		mkdir -p "${SINGULARITY_CACHEDIR}"
		apptainer build --fakeroot "${KK_SIF}" "${REPO_ROOT}/knock-knock.def"
		echo "knock-knock image built: ${KK_SIF}"
	else
		echo "knock-knock Singularity image already exists: ${KK_SIF}"
	fi
fi
