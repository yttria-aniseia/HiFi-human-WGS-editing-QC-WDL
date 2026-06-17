#!/bin/bash
# build_resource_bundle.sh — MAINTAINER-ONLY. Pull every redistributable resource
# from its (slow, flaky) upstream, transform it, lay it out under ONE tree, and
# tar it into the frozen bundle that gets uploaded to Zenodo. Users never run
# this — they run scripts/fetch_resources.sh against the frozen record.
#
# Run this only when bumping versions. Edit resources/manifest.tsv first.
#
#   ./scripts/build_resource_bundle.sh [--out-dir DIR] [--keep-extracted]
#
# Produces:  <out>/editing-qc-resources-<ver>/   (the tree)
#            <out>/editing-qc-resources-<ver>.tar (upload THIS to Zenodo)
#            and prints the bundle sha256 to paste into manifest.tsv.
#
# dbNSFP is license-gated and is deliberately NOT included — see manifest.tsv.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MANIFEST="${REPO_ROOT}/resources/manifest.tsv"
OUT_DIR="${REPO_ROOT}"
KEEP_EXTRACTED=0

while [[ $# -gt 0 ]]; do
	case "$1" in
		--out-dir) OUT_DIR="$2"; shift 2;;
		--keep-extracted) KEEP_EXTRACTED=1; shift;;
		-h|--help) sed -n '2,17p' "$0"; exit 0;;
		*) echo "Unknown arg: $1" >&2; exit 1;;
	esac
done

# --- tiny manifest reader -----------------------------------------------------
m() { # m KEY -> value (first match), errors if missing
	local v
	v=$(grep -vE '^[[:space:]]*#' "$MANIFEST" | awk -F'\t' -v k="$1" '$1==k {print $2; exit}')
	[[ -n "$v" ]] || { echo "manifest: key '$1' not found in $MANIFEST" >&2; exit 1; }
	printf '%s' "$v"
}

VER=$(m bundle_version)
VERp=${VER//./p}
TREE_NAME="editing-qc-resources-${VER}"
TREE="${OUT_DIR}/${TREE_NAME}"
TAR="${OUT_DIR}/${TREE_NAME}.tar"

echo ">>> Building ${TREE_NAME} into ${OUT_DIR}"
mkdir -p "${TREE}"
cd "${TREE}"

# Resumable fetch helper: aria2c (multi-conn) if present, else wget -c.
fetch() { # fetch URL [output_path]
	local url="$1" out="${2:-}"
	if command -v aria2c >/dev/null 2>&1; then
		if [[ -n "$out" ]]; then aria2c -x8 -s8 -c -o "$out" "$url"
		else aria2c -x8 -s8 -c "$url"; fi
	else
		if [[ -n "$out" ]]; then wget -c -O "$out" "$url"
		else wget -c "$url"; fi
	fi
}

# Each step below guards on its final output so the build is idempotent: if it
# aborts partway (or you re-run after fixing one upstream), completed steps are
# skipped and you resume from the failure — no re-downloading the 200 GB of inputs.

# --- 1. core reference suite (brings in the map templates + GRCh38 + slivar) --
if [[ -d "hifi-wdl-resources-${VER}/GRCh38" ]]; then
	echo ">>> hifi-wdl-resources ${VER} already extracted — skipping"
else
	echo ">>> hifi-wdl-resources ${VER}"
	fetch "$(m hifi_wdl_resources_url)" "hifi-wdl-resources-${VER}.tar"
	tar -xf "hifi-wdl-resources-${VER}.tar"
	rm -f "hifi-wdl-resources-${VER}.tar"
fi

# --- 2. hg002 assembly + chain, indexed -------------------------------------
if [[ -f "hifi-wdl-resources-${VER}/hg002/hg002v1.1.fa.gz.gzi" ]]; then
	echo ">>> hg002 assembly already indexed — skipping"
else
	echo ">>> hg002 assembly"
	mkdir -p "hifi-wdl-resources-${VER}/hg002"
	fetch "$(m hg002_chain_url)" "hifi-wdl-resources-${VER}/hg002/hg002v1.1_to_GRCh38.chain.gz"
	fetch "$(m hg002_fasta_url)"  "hifi-wdl-resources-${VER}/hg002/hg002v1.1.fa.gz"
	module load samtools 2>/dev/null || true
	samtools faidx "hifi-wdl-resources-${VER}/hg002/hg002v1.1.fa.gz"
	bgzip -r "hifi-wdl-resources-${VER}/hg002/hg002v1.1.fa.gz"
fi

# --- 3. somatic annotation reference -----------------------------------------
# The family workflow's somatic ANNOTATION path consumes only one file from the
# hifisomatic suite — the reformatted Ensembl GFF3 (somatic_map "ref_gff"). The
# rest of the hifisomatic tarball (the 3 GB GRC-exclusions FASTA, refFlat, TRF
# bed, control VCF, chr.bed) and the Severus 1000G PoN belonged to the somatic
# SV-CALLING (Severus) path, which has been removed from the workflow — so we
# extract the GFF3 only and drop ~3.2 GB of dead weight from the bundle.
SOM_GFF="hifisomatic_resources/Homo_sapiens.GRCh38.112.chr.reformatted.gff3"
if [[ -f "$SOM_GFF" ]]; then
	echo ">>> hifisomatic GFF3 already present — skipping"
else
	echo ">>> hifisomatic GFF3 (ref_gff) — extracting only the file the workflow uses"
	mkdir -p hifisomatic_resources
	fetch "$(m hifisomatic_resources_url)" hifisomatic_resources.tar.gz
	GFF_BASE="$(basename "$SOM_GFF")"
	# extract just the GFF3 wherever it lives in the archive (GNU tar --wildcards)
	tar -xzf hifisomatic_resources.tar.gz -C hifisomatic_resources --wildcards "*${GFF_BASE}"
	rm -f hifisomatic_resources.tar.gz
	# flatten to the expected path if the archive nested it in a subdir
	[[ -f "$SOM_GFF" ]] || find hifisomatic_resources -name "$GFF_BASE" -exec mv {} "$SOM_GFF" \;
	[[ -f "$SOM_GFF" ]] || { echo "ERROR: GFF3 not found in hifisomatic tarball" >&2; exit 1; }
fi

# --- 4. AnnotSV annotation cache (pinned commit + version args) --------------
if [[ -f annotsv_cache.tar.gz ]]; then
	echo ">>> annotsv_cache.tar.gz already built — skipping"
else
	echo ">>> AnnotSV annotations (commit $(m annotsv_install_commit))"
	EXO_VER=$(m annotsv_exomiser_ver)
	# Re-running INSTALL_annotations.sh re-downloads the ~12 GB Exomiser zip even
	# when the cache is already on disk, so only run it if the cache is incomplete.
	if [[ ! -d "AnnotSV_annotations/Annotations_Human" || ! -d "AnnotSV_annotations/Annotations_Exomiser/${EXO_VER}" ]]; then
		fetch "https://raw.githubusercontent.com/lgmgeo/AnnotSV/$(m annotsv_install_commit)/bin/INSTALL_annotations.sh" INSTALL_annotations.sh
		sed -i 's/rm -r AnnotSV/rm -rf AnnotSV/' INSTALL_annotations.sh
		# Upstream bug in this pinned commit: the Exomiser download block references the
		# variable ${EXOMISERVERSION} (no underscore) while the version is stored in
		# ${EXOMISER_VERSION}. The empty expansion requests ".../data/_phenotype.zip"
		# (404, a 162-byte error page) and unzip fails. Fix the typo'd token in place.
		sed -i 's/EXOMISERVERSION/EXOMISER_VERSION/g' INSTALL_annotations.sh
		chmod +x INSTALL_annotations.sh
		bash INSTALL_annotations.sh "$(m annotsv_annotation_ver)" "${EXO_VER}"
	else
		echo "    AnnotSV_annotations already complete — packaging existing cache"
	fi
	tar -czf annotsv_cache.tar.gz AnnotSV_annotations
	rm -rf AnnotSV_annotations INSTALL_annotations.sh
fi

# --- 5. VEP cache ------------------------------------------------------------
echo ">>> VEP cache (release $(m vep_ver))"
[[ -f homo_sapiens_refseq_vep_$(m vep_ver)_GRCh38.tar.gz ]] || fetch "$(m vep_cache_url)"

# --- 6. ClinVar --------------------------------------------------------------
echo ">>> ClinVar"
[[ -f clinvar.vcf.gz ]]     || fetch "$(m clinvar_vcf_url)"
[[ -f clinvar.vcf.gz.tbi ]] || fetch "$(m clinvar_tbi_url)"

# --- 7. stage + normalize map templates to a uniform <prefix>/ scheme --------
# Template sources differ:
#   ref/tertiary  ship INSIDE the hifi-wdl-resources-<ver>/ subdir (the tarball),
#                 already using "<prefix>/hifi-wdl-resources-<ver>/..." paths.
#   somatic       is a maintained, git-tracked repo file with bare paths (it is
#                 NOT in the tarball). We read it from git HEAD so a prior
#                 fetch_resources.sh run (which writes a populated copy to the
#                 repo root under the same name) can't poison the build.
# We stage all three at the tree root, add the hg002 block to ref, and rewrite
# somatic's bare paths to "<prefix>/...". fetch_resources.sh then swaps <prefix>
# for the absolute extracted-tree path. dbNSFP entries stay as a manual placeholder.
TAB=$'\t'
SUBDIR="hifi-wdl-resources-${VER}"
REF_T="GRCh38.ref_map.${VERp}.template.tsv"
SOM_T="GRCh38.somatic_map.${VERp}.template.tsv"

for t in ref tertiary; do
	f="GRCh38.${t}_map.${VERp}.template.tsv"
	if [[ ! -f "$f" ]]; then
		[[ -f "${SUBDIR}/$f" ]] || { echo "missing $f (expected in ${SUBDIR}/ from the reference tarball)" >&2; exit 1; }
		cp "${SUBDIR}/$f" "$f"
	fi
done
if [[ ! -f "$SOM_T" ]]; then
	git -C "${REPO_ROOT}" show "HEAD:${SOM_T}" > "$SOM_T" 2>/dev/null \
		|| { echo "missing $SOM_T (maintained repo file, not in tarball; not found in git HEAD)" >&2; exit 1; }
fi

# append hg002 entries (prefix-rooted) if not already present
if ! grep -q '^hg002_name' "$REF_T"; then
	cat >> "$REF_T" <<EOF
hg002_name${TAB}HG002
hg002_fasta${TAB}<prefix>/hifi-wdl-resources-${VER}/hg002/hg002v1.1.fa.gz
hg002_fasta_index${TAB}<prefix>/hifi-wdl-resources-${VER}/hg002/hg002v1.1.fa.gz.fai
hg002_fasta_bgzf_index${TAB}<prefix>/hifi-wdl-resources-${VER}/hg002/hg002v1.1.fa.gz.gzi
hg002_chain${TAB}<prefix>/hifi-wdl-resources-${VER}/hg002/hg002v1.1_to_GRCh38.chain.gz
EOF
fi

# somatic: prepend <prefix>/ to file-path values; leave scalars and dbNSFP alone
python3 - "$SOM_T" <<'PYEOF'
import sys
path_keys = {"ref_gff","vep_cache","annotsv_cache","clinvar_vcf","clinvar_vcf_index"}
f = sys.argv[1]
out = []
for line in open(f):
    if "\t" in line and not line.startswith("#"):
        k, v = line.rstrip("\n").split("\t", 1)
        if k in path_keys and not v.startswith("<prefix>/"):
            v = "<prefix>/" + v
        out.append(f"{k}\t{v}\n")
    else:
        out.append(line)
open(f, "w").writelines(out)
PYEOF

cp "$MANIFEST" manifest.tsv

# --- 8. checksums + tar -------------------------------------------------------
echo ">>> Writing SHA256SUMS"
find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS

cd "${OUT_DIR}"
echo ">>> Creating ${TAR}"
tar -cf "${TAR}" "${TREE_NAME}"
BUNDLE_SHA=$(sha256sum "${TAR}" | awk '{print $1}')

if [[ ${KEEP_EXTRACTED} -eq 0 ]]; then rm -rf "${TREE}"; fi

cat <<EOF

================================================================================
Bundle built:  ${TAR}
bundle_sha256: ${BUNDLE_SHA}

NEXT STEPS:
  1. Update resources/manifest.tsv:  bundle_sha256<TAB>${BUNDLE_SHA}
  2. Upload ${TREE_NAME}.tar to your Zenodo record, then set
     bundle_zenodo_record to the new record id.
  3. Commit manifest.tsv.
dbNSFP is NOT in this bundle (license-gated); users obtain it via fetch_resources.sh.
================================================================================
EOF
