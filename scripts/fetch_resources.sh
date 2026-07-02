#!/bin/bash
# fetch_resources.sh — USER-FACING. Pull the frozen resource bundle from Zenodo
# (resumable, checksum-verified), extract it under one tree, and write ready-to-use
# map files with absolute paths. No upstream servers, no transforms — all the slow,
# flaky work was done once by build_resource_bundle.sh.
#
#   ./scripts/fetch_resources.sh [--dest DIR] [--dbnsfp FILE] [--keep-tar]
#
#   --dest DIR    where to extract the bundle (default: repo root)
#   --dbnsfp FILE path to your licensed dbNSFP bgzip (+ .tbi alongside); the
#                 somatic map is patched to point at it. Omit to leave a clearly
#                 marked placeholder you fill in later.
#
# If a download is interrupted, just re-run — it resumes from where it stopped.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MANIFEST="${REPO_ROOT}/resources/manifest.tsv"
DEST="${REPO_ROOT}"
DBNSFP=""
KEEP_TAR=0

while [[ $# -gt 0 ]]; do
	case "$1" in
		--dest) DEST="$2"; shift 2;;
		--dbnsfp) DBNSFP="$2"; shift 2;;
		--keep-tar) KEEP_TAR=1; shift;;
		-h|--help) sed -n '2,17p' "$0"; exit 0;;
		*) echo "Unknown arg: $1" >&2; exit 1;;
	esac
done

m() {
	local v
	v=$(grep -vE '^[[:space:]]*#' "$MANIFEST" | awk -F'\t' -v k="$1" '$1==k {print $2; exit}')
	[[ -n "$v" ]] || { echo "manifest: key '$1' not found in $MANIFEST" >&2; exit 1; }
	printf '%s' "$v"
}

VER=$(m bundle_version)
VERp=${VER//./p}
TAR_NAME=$(m bundle_tar)
RECORD=$(m bundle_zenodo_record)
WANT_SHA=$(m bundle_sha256)
TREE_NAME="editing-qc-resources-${VER}"
URL="https://zenodo.org/records/${RECORD}/files/${TAR_NAME}"

mkdir -p "${DEST}"
cd "${DEST}"

echo ">>> Fetching ${TAR_NAME} from Zenodo record ${RECORD} (resumable)"
if command -v aria2c >/dev/null 2>&1; then
	aria2c -x8 -s8 -c -o "${TAR_NAME}" "${URL}"
else
	echo "    (install aria2c for faster multi-connection downloads; using wget -c)"
	wget -c -O "${TAR_NAME}" "${URL}"
fi

if [[ "${WANT_SHA}" != "TBD" ]]; then
	echo ">>> Verifying bundle sha256"
	echo "${WANT_SHA}  ${TAR_NAME}" | sha256sum -c - || {
		echo "ERROR: checksum mismatch. Delete ${DEST}/${TAR_NAME} and re-run." >&2; exit 1; }
else
	echo "WARNING: manifest bundle_sha256 is TBD — skipping bundle verification." >&2
fi

echo ">>> Extracting"
tar -xf "${TAR_NAME}"
[[ ${KEEP_TAR} -eq 1 ]] || rm -f "${TAR_NAME}"

TREE="$(cd "${TREE_NAME}" && pwd)"   # absolute path

if [[ -f "${TREE}/SHA256SUMS" ]]; then
	echo ">>> Verifying extracted files against SHA256SUMS"
	( cd "${TREE}" && sha256sum -c --quiet SHA256SUMS )
fi

# --- write ready-to-use map files with absolute paths ------------------------
# Output keeps the historical *.template.tsv name (now holding absolute, populated
# paths) so existing input configs and docs that reference it keep working — the
# in-bundle copy under ${TREE} is the pristine <prefix> version.
echo ">>> Writing map files (absolute paths -> ${TREE})"
for kind in ref tertiary somatic; do
	src="${TREE}/GRCh38.${kind}_map.${VERp}.template.tsv"
	dst="${DEST}/GRCh38.${kind}_map.${VERp}.template.tsv"
	[[ -f "$src" ]] || { echo "missing template $src in bundle" >&2; exit 1; }
	sed "s|<prefix>|${TREE}|g" "$src" > "$dst"
	echo "    ${dst}"
done

# --- dbNSFP (license-gated, not in bundle) -----------------------------------
SOM="${DEST}/GRCh38.somatic_map.${VERp}.template.tsv"
if [[ -n "${DBNSFP}" ]]; then
	DBN_ABS="$(cd "$(dirname "${DBNSFP}")" && pwd)/$(basename "${DBNSFP}")"
	[[ -f "${DBN_ABS}" ]]       || { echo "ERROR: --dbnsfp file not found: ${DBN_ABS}" >&2; exit 1; }
	[[ -f "${DBN_ABS}.tbi" ]]   || echo "WARNING: ${DBN_ABS}.tbi not found next to dbNSFP file" >&2
	sed -i -E "s|^(dbnsfp_file)\t.*|\1\t${DBN_ABS}|; s|^(dbnsfp_file_index)\t.*|\1\t${DBN_ABS}.tbi|" "${SOM}"
	echo ">>> Patched dbNSFP path into ${SOM}"
else
	cat >&2 <<EOF

>>> dbNSFP is license-gated and is NOT in the bundle.
    Obtain dbNSFP $(m dbnsfp_ver) from $(m dbnsfp_download_page), build the
    bgzip+tabix grch38 file, then re-run with:  --dbnsfp /path/to/dbNSFP_grch38.gz
    Until then, the dbnsfp_file entries in ${SOM} are placeholders.
EOF
fi

echo ">>> Done. Map files are in ${DEST}; resources under ${TREE}"
