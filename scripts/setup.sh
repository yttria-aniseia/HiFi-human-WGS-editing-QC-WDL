#!/bin/bash


#download reference tarball for singleton
echo "Starting Reference suite download"
wget "https://zenodo.org/record/15750792/files/hifi-wdl-resources-v3.0.0.tar"
echo "Reference Tarball Sucessfully downloaded"
echo "Extracting tarball"
tar -xvf hifi-wdl-resources-v3.0.0.tar
echo "Finished Extracting tarball"

mkdir -p hifi-wdl-resources-v3.0.0/hg002
wget -P hifi-wdl-resources-v3.0.0/hg002/ https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/changes/hg002v1.1_to_GRCh38.chain.gz
wget -P hifi-wdl-resources-v3.0.0/hg002/ https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.1.fasta.gz
mv hifi-wdl-resources-v3.0.0/hg002/hg002v1.1.fasta.gz hifi-wdl-resources-v3.0.0/hg002/hg002v1.1.fa.gz
module load samtools
samtools faidx hifi-wdl-resources-v3.0.0/hg002/hg002v1.1.fa.gz
bgzip -r hifi-wdl-resources-v3.0.0/hg002/hg002v1.1.fa.gz
if [[ ! -e GRCh38.ref_map.v3p0p0.template.tsv.bak ]]; then
	cp GRCh38.ref_map.v3p0p0.template.tsv GRCh38.ref_map.v3p0p0.template.tsv.bak
	cp GRCh38.tertiary_map.v3p0p0.template.tsv GRCh38.tertiary_map.v3p0p0.template.tsv.bak
	cp GRCh38.somatic_map.v3p0p0.template.tsv GRCh38.somatic_map.v3p0p0.template.tsv.bak
else
	cp GRCh38.ref_map.v3p0p0.template.tsv.bak GRCh38.ref_map.v3p0p0.template.tsv
	cp GRCh38.tertiary_map.v3p0p0.template.tsv.bak GRCh38.tertiary_map.v3p0p0.template.tsv
	cp GRCh38.somatic_map.v3p0p0.template.tsv.bak GRCh38.somatic_map.v3p0p0.template.tsv
fi
sed -i "s/<prefix>\///g" GRCh38.ref_map.v3p0p0.template.tsv
sed -i "s/<prefix>\///g" GRCh38.tertiary_map.v3p0p0.template.tsv
sed -i "s/<prefix>\///g" GRCh38.somatic_map.v3p0p0.template.tsv
echo "hg002_name	HG002
hg002_fasta	hifi-wdl-resources-v3.0.0/hg002/hg002v1.1.fa.gz
hg002_fasta_index	hifi-wdl-resources-v3.0.0/hg002/hg002v1.1.fa.gz.fai
hg002_fasta_bgzf_index	hifi-wdl-resources-v3.0.0/hg002/hg002v1.1.fa.gz.gzi
hg002_chain	hifi-wdl-resources-v3.0.0/hg002/hg002v1.1_to_GRCh38.chain.gz" >> GRCh38.ref_map.v3p0p0.template.tsv

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
wget https://github.com/lgmgeo/AnnotSV/raw/master/bin/INSTALL_annotations.sh
echo "Install Script Sucessfully downloaded"
echo "Starting AnnotSV cache download"
sed -i 's/rm -r AnnotSV/rm -rf AnnotSV/' INSTALL_annotations.sh
chmod +x INSTALL_annotations.sh
bash INSTALL_annotations.sh
echo "Finished download cache"
tar -czvf annotsv_cache.tar.gz AnnotSV_annotations

echo "Starting VEP suite download"
wget https://ftp.ensembl.org/pub/release-112/variation/indexed_vep_cache/homo_sapiens_refseq_vep_112_GRCh38.tar.gz
echo "Reference Tarball Sucessfully downloaded"



#prepull images
bash scripts/create_image_manifest.sh
bash scripts/pre-pull.sh
