#!/bin/bash


VER="v3.1.0"
VERp=${VER//./p}
FETCH_EXTRA=1
ZENODO_BASE="https://zenodo.org/records/17086906"

#download reference tarball for singleton
echo "Starting Reference suite download"
wget "${ZENODO_BASE}/files/hifi-wdl-resources-${VER}.tar"
echo "Reference Tarball Sucessfully downloaded"
echo "Extracting tarball"
tar -xvf hifi-wdl-resources-${VER}.tar
echo "Finished Extracting tarball"

mkdir -p hifi-wdl-resources-${VER}/hg002
wget -P hifi-wdl-resources-${VER}/hg002/ https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/changes/hg002v1.1_to_GRCh38.chain.gz
wget -P hifi-wdl-resources-${VER}/hg002/ https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.1.fasta.gz
mv hifi-wdl-resources-${VER}/hg002/hg002v1.1.fasta.gz hifi-wdl-resources-${VER}/hg002/hg002v1.1.fa.gz
module load samtools
samtools faidx hifi-wdl-resources-${VER}/hg002/hg002v1.1.fa.gz
bgzip -r hifi-wdl-resources-${VER}/hg002/hg002v1.1.fa.gz
if [[ ! -e GRCh38.ref_map.${VERp}.template.tsv.bak ]]; then
	cp GRCh38.ref_map.${VERp}.template.tsv GRCh38.ref_map.${VERp}.template.tsv.bak
	cp GRCh38.tertiary_map.${VERp}.template.tsv GRCh38.tertiary_map.${VERp}.template.tsv.bak
	cp GRCh38.somatic_map.${VERp}.template.tsv GRCh38.somatic_map.${VERp}.template.tsv.bak
else
	cp GRCh38.ref_map.${VERp}.template.tsv.bak GRCh38.ref_map.${VERp}.template.tsv
	cp GRCh38.tertiary_map.${VERp}.template.tsv.bak GRCh38.tertiary_map.${VERp}.template.tsv
	cp GRCh38.somatic_map.${VERp}.template.tsv.bak GRCh38.somatic_map.${VERp}.template.tsv
fi
sed -i "s/<prefix>\///g" GRCh38.ref_map.${VERp}.template.tsv
sed -i "s/<prefix>\///g" GRCh38.tertiary_map.${VERp}.template.tsv
sed -i "s/<prefix>\///g" GRCh38.somatic_map.${VERp}.template.tsv
TAB=$'\t'
echo "hg002_name${TAB}HG002
hg002_fasta${TAB}hifi-wdl-resources-${VER}/hg002/hg002v1.1.fa.gz
hg002_fasta_index${TAB}hifi-wdl-resources-${VER}/hg002/hg002v1.1.fa.gz.fai
hg002_fasta_bgzf_index${TAB}hifi-wdl-resources-${VER}/hg002/hg002v1.1.fa.gz.gzi
hg002_chain${TAB}hifi-wdl-resources-${VER}/hg002/hg002v1.1_to_GRCh38.chain.gz" >> GRCh38.ref_map.${VERp}.template.tsv
DBNSFP_VERSION="4.9a"


# Build knock-knock Singularity image (no pre-built container available)
if [[ -z "${SINGULARITY_CACHEDIR}" ]]; then
	echo "WARNING: SINGULARITY_CACHEDIR is not set; skipping knock-knock container build"
else
	KK_SIF="${SINGULARITY_CACHEDIR}/docker___knock-knock.sif"
	if [[ ! -f "${KK_SIF}" ]]; then
		echo "Building knock-knock Singularity image..."
		mkdir -p "${SINGULARITY_CACHEDIR}"
		apptainer build --fakeroot "${KK_SIF}" knock-knock.def
		echo "knock-knock image built: ${KK_SIF}"
	else
		echo "knock-knock Singularity image already exists: ${KK_SIF}"
	fi
fi

if [[ $FETCH_EXTRA -eq 1 ]]; then
	echo "Starting somatic suite download"
	wget https://zenodo.org/record/14847828/files/hifisomatic_resources.tar.gz
	echo "Reference Tarball Sucessfully downloaded"
	echo "Extracting tarball"
	mkdir -p hifisomatic_resources
	tar -xvzf hifisomatic_resources.tar.gz -C hifisomatic_resources
	echo "Finished Extracting tarball"

	echo "Starting 1000G panel of normal from Severus GitHub"
	wget https://github.com/KolmogorovLab/Severus/raw/refs/heads/main/pon/PoN_1000G_hg38_extended.tsv.gz
	echo "Reference Tarball Sucessfully downloaded"


	echo "Starting AnnotSV suite download"
	wget https://raw.githubusercontent.com/lgmgeo/AnnotSV/b270de3f6db45e4c4ad6b32e7fc868f2369b62c3/bin/INSTALL_annotations.sh
	echo "Install Script Sucessfully downloaded"
	echo "Starting AnnotSV cache download"
	sed -i 's/rm -r AnnotSV/rm -rf AnnotSV/' INSTALL_annotations.sh
	chmod +x INSTALL_annotations.sh
	bash INSTALL_annotations.sh "3.5" "2406"
	echo "Finished download cache"
	tar -czvf annotsv_cache.tar.gz AnnotSV_annotations

	echo "Starting VEP suite download"
	wget https://ftp.ensembl.org/pub/release-115/variation/indexed_vep_cache/homo_sapiens_refseq_vep_115_GRCh38.tar.gz
	echo "VEP cache successfully downloaded"

	echo "Starting dbNSFP download (BayesDel, REVEL, AlphaMissense, CADD)"
	# dbNSFP requires registration and license agreement; academic use is free but commercial
	# use requires a separate license. See: https://www.dbnsfp.org/download
	# v4.9a is the most recent version available for automated download without registration.
	# Newer versions (recommended) must be downloaded manually after obtaining access — update
	# DBNSFP_VERSION and the dbnsfp_file entries in your somatic map accordingly.
	if [[ "${DBNSFP_VERSION}" == "4.9a" ]]; then
		echo "WARNING: downloading dbNSFP v4.9a (automated baseline). A newer version is available at https://www.dbnsfp.org/download — download it manually and update DBNSFP_VERSION and your somatic map." >&2
	fi
	wget https://dbnsfp.s3.amazonaws.com/dbNSFP${DBNSFP_VERSION}.zip
	unzip dbNSFP${DBNSFP_VERSION}.zip "dbNSFP${DBNSFP_VERSION}_variant.chr*.gz" "dbNSFP${DBNSFP_VERSION}_variant.chr*.gz.tbi"
	# Concatenate chromosomes into single bgzipped file expected by VEP plugin
	(zcat dbNSFP${DBNSFP_VERSION}_variant.chr1.gz | head -1; \
	 for f in dbNSFP${DBNSFP_VERSION}_variant.chr*.gz; do zcat "$f" | tail -n +2; done) \
	  | bgzip -c > dbNSFP${DBNSFP_VERSION}_grch38.gz
	tabix -s 1 -b 2 -e 2 dbNSFP${DBNSFP_VERSION}_grch38.gz
	rm -f dbNSFP${DBNSFP_VERSION}_variant.chr*.gz dbNSFP${DBNSFP_VERSION}_variant.chr*.gz.tbi dbNSFP${DBNSFP_VERSION}.zip
	echo "dbNSFP successfully downloaded and indexed"

	echo "Starting ClinVar download"
	wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
	wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
	echo "ClinVar successfully downloaded"
fi
