version 1.1

import "../wdl-common/wdl/structs.wdl"

task create_sample_reference_consensus {
  meta {
    description: "Create sample-specific reference consensus for edit region using bcftools consensus with sample variants"
  }

  parameter_meta {
    ref_fasta: "Reference genome FASTA"
    ref_fasta_index: "Reference genome FASTA index"
    small_variant_vcf: "Sample small variant VCF (phased)"
    small_variant_vcf_index: "Index for small variant VCF"
    sv_vcf: "Sample structural variant VCF (phased)"
    sv_vcf_index: "Index for structural variant VCF"
    crispr_edit_json: "JSON file describing CRISPR edit with target region"
  }

  input {
    File ref_fasta
    File ref_fasta_index
    File small_variant_vcf
    File small_variant_vcf_index
    File sv_vcf
    File sv_vcf_index
    File crispr_edit_json
    String sample_id
    RuntimeAttributes runtime_attributes
  }

  String out_prefix = "~{sample_id}_sample_ref"
  Int padding = 5000  # Add padding around edit region for better alignment context

  command <<<
    set -euxo pipefail

    python3 <<'EOF'
import json

# Read CRISPR edit JSON to get target region
with open("~{crispr_edit_json}", 'r') as f:
    crispr_data = json.load(f)

# Get first edit's target region
edit = crispr_data['edits'][0]
chr = edit['target']['chr']
start = max(1, edit['target']['start'] - ~{padding})  # Add padding but don't go below 1
end = edit['target']['end'] + ~{padding}
symbol = edit['target']['symbol']

# Write region info for bash
with open('region.txt', 'w') as f:
    f.write(f"{chr}:{start}-{end}\n")
with open('symbol.txt', 'w') as f:
    f.write(f"{symbol}\n")
with open('chr.txt', 'w') as f:
    f.write(f"{chr}\n")
EOF

    REGION=$(cat region.txt)
    SYMBOL=$(cat symbol.txt)
    CHR=$(cat chr.txt)

    echo "Creating sample-specific reference for region: ${REGION}"

    # Extract region from reference genome
    samtools faidx ~{ref_fasta} "${REGION}" > region_ref.fasta

    # Index the small variant VCF if needed and extract region
    tabix -p vcf ~{small_variant_vcf} || true
    bcftools view -r "${REGION}" ~{small_variant_vcf} -O z -o small_variants_region.vcf.gz
    tabix -p vcf small_variants_region.vcf.gz

    # Index the SV VCF if needed and extract region
    tabix -p vcf ~{sv_vcf} || true
    bcftools view -r "${REGION}" ~{sv_vcf} -O z -o sv_region.vcf.gz
    tabix -p vcf sv_region.vcf.gz

    # Merge small variants and SVs
    bcftools concat -a -D small_variants_region.vcf.gz sv_region.vcf.gz -O z -o merged_variants.vcf.gz
    tabix -p vcf merged_variants.vcf.gz

    # Filter out symbolic alleles (except <DEL>, <*>, <NON_REF> which are supported)
    # This removes alleles like <INV>, <DUP>, <CNV>, <INS>, etc.
    bcftools view -e 'ALT~"<INV>" || ALT~"<DUP>" || ALT~"<CNV>" || ALT~"<INS>"' \
      merged_variants.vcf.gz -O z -o merged_variants_filtered.vcf.gz
    tabix -p vcf merged_variants_filtered.vcf.gz

    echo "Filtered VCF stats:"
    echo "  Total variants: $(bcftools view -H merged_variants.vcf.gz | wc -l)"
    echo "  After filtering: $(bcftools view -H merged_variants_filtered.vcf.gz | wc -l)"

    # Apply variants to reference using bcftools consensus
    # This creates a sample-specific reference that includes the variants
    cat region_ref.fasta | bcftools consensus merged_variants_filtered.vcf.gz > sample_ref_consensus.fasta

    # Rename header to include sample info
    python3 <<'EOF'
symbol = open('symbol.txt').read().strip()
region = open('region.txt').read().strip()

with open('sample_ref_consensus.fasta', 'r') as f:
    lines = f.readlines()

# Replace header
with open('~{out_prefix}_consensus.fasta', 'w') as out:
    if lines:
        out.write(f">{symbol}_sample_ref_{region}\n")
        for line in lines[1:]:
            out.write(line)
EOF

    # Index the consensus
    samtools faidx ~{out_prefix}_consensus.fasta

    echo "Sample-specific reference consensus created for ${SYMBOL}"
  >>>

  output {
    File sample_ref_consensus_fasta = "~{out_prefix}_consensus.fasta"
    File sample_ref_consensus_fasta_index = "~{out_prefix}_consensus.fasta.fai"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87"
    cpu: 2
    memory: "8 GB"
    disk: "20 GB"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
  }
}

task extract_region_reads {
  meta {
    description: "Extract aligned reads from edit region preserving haplotype tags"
  }

  parameter_meta {
    haplotagged_bam: "Haplotagged BAM file from phasing"
    haplotagged_bam_index: "Index for haplotagged BAM"
    crispr_edit_json: "JSON file describing CRISPR edit with target region"
  }

  input {
    File haplotagged_bam
    File haplotagged_bam_index
    File crispr_edit_json
    String sample_id
    Int threads = 4
    RuntimeAttributes runtime_attributes
  }

  String out_prefix = "~{sample_id}_region_reads"
  Int padding = 5000  # Must match padding in create_sample_reference_consensus

  command <<<
    set -euxo pipefail

    python3 <<'EOF'
import json

# Read CRISPR edit JSON to get target region
with open("~{crispr_edit_json}", 'r') as f:
    crispr_data = json.load(f)

# Get first edit's target region
edit = crispr_data['edits'][0]
chr = edit['target']['chr']
start = max(1, edit['target']['start'] - ~{padding})
end = edit['target']['end'] + ~{padding}
symbol = edit['target']['symbol']

# Write region info for bash
with open('region.txt', 'w') as f:
    f.write(f"{chr}:{start}-{end}\n")
with open('symbol.txt', 'w') as f:
    f.write(f"{symbol}\n")
EOF

    REGION=$(cat region.txt)
    SYMBOL=$(cat symbol.txt)

    echo "Extracting reads from region: ${REGION}"

    # Extract reads from the edit region (with haplotype tags preserved)
    samtools view -h -@ ~{threads} ~{haplotagged_bam} "${REGION}" | \
      samtools fastq -T HP,PS - | \
      gzip > ~{out_prefix}.fastq.gz

    # Count extracted reads
    READ_COUNT=$(zcat ~{out_prefix}.fastq.gz | grep -c "^@" || echo 0)
    echo "Extracted ${READ_COUNT} reads from ${REGION}"
  >>>

  output {
    File region_reads_fastq = "~{out_prefix}.fastq.gz"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87"
    cpu: threads
    memory: "8 GB"
    disk: "50 GB"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
  }
}

task minimap2_align_to_sample_reference {
  meta {
    description: "Align extracted reads to sample-specific reference using minimap2"
  }

  parameter_meta {
    region_reads_fastq: "FASTQ reads extracted from edit region"
    sample_ref_consensus_fasta: "Sample-specific reference consensus FASTA"
  }

  input {
    File region_reads_fastq
    File sample_ref_consensus_fasta
    String sample_id
    Int threads = 4
    RuntimeAttributes runtime_attributes
  }

  String out_prefix = "~{sample_id}_realigned"

  command <<<
    set -euxo pipefail

    # Align reads to sample-specific reference using minimap2
    minimap2 -ax map-hifi -t ~{threads} -Y \
      ~{sample_ref_consensus_fasta} \
      ~{region_reads_fastq} > ~{out_prefix}.sam

    echo "Alignment complete"
  >>>

  output {
    File realigned_sam = "~{out_prefix}.sam"
  }

  runtime {
    docker: "quay.io/biocontainers/minimap2@sha256:fdc9ef8bfbd31bab59a61b2e90dd226647deed971556175a4cd004f0bcdc7608"
    cpu: threads
    memory: "16 GB"
    disk: "50 GB"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
  }
}

task sam_to_sorted_tagged_bam {
  meta {
    description: "Convert SAM to sorted BAM, generate stats, and transfer haplotype tags from original BAM"
  }

  parameter_meta {
    sam_file: "SAM alignment file"
    original_bam: "Original haplotagged BAM (source of HP/PS tags)"
    original_bam_index: "Index for original BAM"
    crispr_edit_json: "JSON file with region info"
  }

  input {
    File sam_file
    File original_bam
    File original_bam_index
    File crispr_edit_json
    String sample_id
    Int threads = 4
    RuntimeAttributes runtime_attributes
  }

  String out_prefix = "~{sample_id}_realigned"
  String tagged_prefix = "~{sample_id}_realigned_tagged"
  Int padding = 5000

  command <<<
    set -euxo pipefail

    # Step 1: Convert to BAM, sort, and index
    samtools view -b -@ ~{threads} ~{sam_file} | \
      samtools sort -@ ~{threads} -o ~{out_prefix}.bam

    samtools index -@ ~{threads} ~{out_prefix}.bam

    # Count aligned reads and generate stats
    ALIGNED_COUNT=$(samtools view -c ~{out_prefix}.bam)
    echo "Sorted ${ALIGNED_COUNT} reads"

    samtools flagstat ~{out_prefix}.bam > ~{out_prefix}_flagstat.txt

    # Step 2: Transfer haplotype tags from original BAM
    python3 <<'EOF'
import json
import pysam

# Read CRISPR edit JSON to get target region
with open("~{crispr_edit_json}", 'r') as f:
    crispr_data = json.load(f)

edit = crispr_data['edits'][0]
chr = edit['target']['chr']
start = max(1, edit['target']['start'] - ~{padding})
end = edit['target']['end'] + ~{padding}

# Open BAMs
original = pysam.AlignmentFile("~{original_bam}", "rb")
realigned = pysam.AlignmentFile("~{out_prefix}.bam", "rb")

# Create output BAM with same header as realigned
output = pysam.AlignmentFile("~{tagged_prefix}.bam", "wb", template=realigned)

# Build dict of read_name -> (HP, PS) from original BAM in the region
tags_by_read = {}
for read in original.fetch(chr, start, end):
    hp_tag = read.get_tag("HP") if read.has_tag("HP") else None
    ps_tag = read.get_tag("PS") if read.has_tag("PS") else None
    if hp_tag is not None or ps_tag is not None:
        tags_by_read[read.query_name] = (hp_tag, ps_tag)

print(f"Extracted tags from {len(tags_by_read)} reads in original BAM")

# Transfer tags to realigned reads
transferred = 0
for read in realigned.fetch():
    if read.query_name in tags_by_read:
        hp_tag, ps_tag = tags_by_read[read.query_name]
        if hp_tag is not None:
            read.set_tag("HP", hp_tag, value_type='i')
        if ps_tag is not None:
            read.set_tag("PS", ps_tag, value_type='i')
        transferred += 1
    output.write(read)

print(f"Transferred tags to {transferred} reads in realigned BAM")

original.close()
realigned.close()
output.close()
EOF

    # Index the tagged BAM
    samtools index ~{tagged_prefix}.bam

    echo "Combined SAM conversion, sorting, and tag transfer complete"
  >>>

  output {
    File sorted_bam = "~{out_prefix}.bam"
    File sorted_bam_index = "~{out_prefix}.bam.bai"
    File alignment_stats = "~{out_prefix}_flagstat.txt"
    File tagged_bam = "~{tagged_prefix}.bam"
    File tagged_bam_index = "~{tagged_prefix}.bam.bai"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87"
    cpu: threads
    memory: "8 GB"
    disk: "50 GB"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
  }
}

task create_haplotype_consensus_from_realigned {
  meta {
    description: "Create haplotype consensus sequences from reads realigned to sample-specific reference"
  }

  parameter_meta {
    realigned_bam: "BAM file with reads realigned to sample-specific reference"
    realigned_bam_index: "Index for realigned BAM"
    sample_ref_consensus_fasta: "Sample-specific reference consensus FASTA (used as coordinate system)"
    crispr_edit_json: "JSON file describing CRISPR edit with target region"
  }

  input {
    File realigned_bam
    File realigned_bam_index
    File sample_ref_consensus_fasta
    File crispr_edit_json
    String sample_id
    RuntimeAttributes runtime_attributes
  }

  String out_prefix = "~{sample_id}_haplotype_consensus"

  command <<<
    set -euxo pipefail

    python3 <<'EOF'
import json

# Read CRISPR edit JSON to get target info
with open("~{crispr_edit_json}", 'r') as f:
    crispr_data = json.load(f)

# Get first edit's metadata
edit = crispr_data['edits'][0]
symbol = edit['target']['symbol']

# Write symbol for bash
with open('symbol.txt', 'w') as f:
    f.write(f"{symbol}\n")
EOF

    SYMBOL=$(cat symbol.txt)

    echo "Creating haplotype consensus sequences for ${SYMBOL}"

    # Get the sequence name from the sample reference (should be first/only sequence)
    REF_NAME=$(head -n 1 ~{sample_ref_consensus_fasta} | cut -d' ' -f1 | sed 's/>//')
    echo "Reference sequence name: ${REF_NAME}"

    # Split BAM by haplotype tag (HP)
    samtools split -d HP --output-fmt BAM --write-index -f "hap%#.%." -u "unphased.bam" ~{realigned_bam}

    # Check if haplotype files were created and generate consensus for each
    if [ -f "hap0.bam" ]; then
      echo "Generating consensus for haplotype 1..."
      samtools consensus -a -m simple --min-MQ 10 --min-BQ 10 -f fasta -X hifi hap0.bam > hap1_consensus.fasta || echo ">hap1_no_consensus" > hap1_consensus.fasta

      # Get depth info
      samtools depth hap0.bam | awk '{sum+=$3; count++} END {if(count>0) print "Haplotype 1 mean depth:", sum/count; else print "Haplotype 1 mean depth: 0"}' > hap1_depth.txt
      cat hap1_depth.txt
    else
      echo "No haplotype 1 reads found"
      echo ">hap1_no_consensus" > hap1_consensus.fasta
      echo "Haplotype 1 mean depth: 0" > hap1_depth.txt
    fi

    if [ -f "hap1.bam" ]; then
      echo "Generating consensus for haplotype 2..."
      samtools consensus -a -m simple --min-MQ 10 --min-BQ 10 -f fasta -X hifi hap1.bam > hap2_consensus.fasta || echo ">hap2_no_consensus" > hap2_consensus.fasta

      # Get depth info
      samtools depth hap1.bam | awk '{sum+=$3; count++} END {if(count>0) print "Haplotype 2 mean depth:", sum/count; else print "Haplotype 2 mean depth: 0"}' > hap2_depth.txt
      cat hap2_depth.txt
    else
      echo "No haplotype 2 reads found"
      echo ">hap2_no_consensus" > hap2_consensus.fasta
      echo "Haplotype 2 mean depth: 0" > hap2_depth.txt
    fi

    # Combine consensus sequences with informative headers
    python3 <<'EOF'
import sys

symbol = open('symbol.txt').read().strip()

# Combine haplotype consensus files
with open('~{out_prefix}_combined.fasta', 'w') as out:
    # Haplotype 1
    with open('hap1_consensus.fasta', 'r') as f:
        lines = f.readlines()
        if len(lines) > 1:
            seq = ''.join(line.strip() for line in lines[1:])
            out.write(f">{symbol}_hap1_realigned\n{seq}\n")
        else:
            out.write(f">{symbol}_hap1_no_consensus\nN\n")

    # Haplotype 2
    with open('hap2_consensus.fasta', 'r') as f:
        lines = f.readlines()
        if len(lines) > 1:
            seq = ''.join(line.strip() for line in lines[1:])
            out.write(f">{symbol}_hap2_realigned\n{seq}\n")
        else:
            out.write(f">{symbol}_hap2_no_consensus\nN\n")

print("Combined haplotype consensus sequences created")
EOF

    # Create summary file
    echo "Haplotype Consensus Summary for ${SYMBOL}" > ~{out_prefix}_summary.txt
    echo "==========================================" >> ~{out_prefix}_summary.txt
    echo "" >> ~{out_prefix}_summary.txt
    cat hap1_depth.txt >> ~{out_prefix}_summary.txt
    cat hap2_depth.txt >> ~{out_prefix}_summary.txt
    echo "" >> ~{out_prefix}_summary.txt
    echo "Consensus sequences:" >> ~{out_prefix}_summary.txt
    grep -c "^>" ~{out_prefix}_combined.fasta >> ~{out_prefix}_summary.txt || echo "0" >> ~{out_prefix}_summary.txt
  >>>

  output {
    File haplotype_consensus_fasta = "~{out_prefix}_combined.fasta"
    File consensus_summary = "~{out_prefix}_summary.txt"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87"
    cpu: 2
    memory: "8 GB"
    disk: "20 GB"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
  }
}

task annotate_consensus_genbank {
  meta {
    description: "Create GenBank files annotating edit parts on haplotype consensus sequences"
  }

  parameter_meta {
    parts_alignment_paf: "PAF file from minimap2 alignment of edit parts to consensus sequences"
    consensus_fasta: "FASTA file with haplotype consensus sequences"
    crispr_edit_json: "JSON file describing CRISPR edit structure"
    coverage_threshold: "Minimum coverage threshold for faithful alignments"
  }

  input {
    File parts_alignment_paf
    File consensus_fasta
    File crispr_edit_json
    Float coverage_threshold = 0.8
    String sample_id
    RuntimeAttributes runtime_attributes
  }

  String out_prefix = "~{sample_id}_edit_annotations"

  command <<<
    set -euxo pipefail

    cat > create_genbank_annotations.py <<'EOF'
#!/usr/bin/env python3
"""
create_genbank_annotations.py - Annotate edit parts on haplotype consensus sequences

Creates GenBank files for each haplotype consensus sequence with proper source and gene
features indicating genomic origin, plus edit part annotations.
"""

import argparse
import json
import re
from collections import defaultdict
from datetime import datetime
from typing import Dict, List, Tuple
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

# Type aliases
Hit = Tuple[int, int, int, int, str, int]  # q_start, q_end, t_start, t_end, strand, n_match

def parse_paf_line(line: str) -> Dict:
    """Parse a single PAF format line"""
    f = line.rstrip("\n").split("\t")
    if len(f) < 12:
        raise ValueError("PAF line has fewer than 12 mandatory columns")
    return dict(
        query_name=f[0], query_len=int(f[1]),
        q_start=int(f[2]), q_end=int(f[3]), strand=f[4],
        target_name=f[5], target_len=int(f[6]),
        t_start=int(f[7]), t_end=int(f[8]),
        n_match=int(f[9]), n_align=int(f[10]), mapq=int(f[11]),
    )

def load_faithful_alignments(paf_path: str, cov_thr: float) -> Tuple[Dict[str, Dict[str, List[Hit]]], Dict[str, int]]:
    """
    Load only faithful alignments (high coverage, low gaps) from PAF.
    Returns: ({target_consensus_name: {query_part_name: [Hit, ...]}}, {query_name: length})
    """
    faithful = defaultdict(lambda: defaultdict(list))
    qlen = {}

    with open(paf_path) as fh:
        for ln in fh:
            if not ln.strip():
                continue
            r = parse_paf_line(ln)
            qlen.setdefault(r["query_name"], r["query_len"])

            # Calculate coverage and gap metrics
            cov_hit = r["n_match"] / r["query_len"]
            gap_ok = abs(r["n_match"] - r["n_align"]) <= 0.05 * r["n_match"]

            # Only keep faithful alignments
            if cov_hit >= cov_thr and gap_ok:
                hit: Hit = (
                    r["q_start"], r["q_end"],
                    r["t_start"], r["t_end"],
                    r["strand"], r["n_match"]
                )
                faithful[r["target_name"]][r["query_name"]].append(hit)

    return faithful, qlen

def load_edit_metadata(json_path: str) -> Dict:
    """Load edit metadata from CRISPR edit JSON"""
    with open(json_path) as f:
        data = json.load(f)

    # Get first edit's metadata
    edit = data['edits'][0]
    return {
        'edit_id': edit['id'],
        'chromosome': edit['target']['chr'],
        'start': edit['target']['start'],
        'end': edit['target']['end'],
        'symbol': edit['target']['symbol']
    }

def create_source_and_gene_features(seq_len: int, metadata: Dict, haplotype: str) -> List[SeqFeature]:
    """
    Create source and gene features for the genomic context of the consensus sequence.

    Args:
        seq_len: Length of the consensus sequence
        metadata: Edit metadata with genomic coordinates
        haplotype: 'hap1' or 'hap2'

    Returns:
        List containing source and gene SeqFeature objects
    """
    features = []

    # Source feature covering entire sequence
    source_feature = SeqFeature(
        location=FeatureLocation(0, seq_len, strand=+1),
        type="source",
        qualifiers={
            "organism": ["Homo sapiens"],
            "mol_type": ["genomic DNA"],
            "db_xref": ["taxon:9606"],
            "chromosome": [metadata['chromosome'].replace('chr', '')],
            "map": [f"{metadata['start']}-{metadata['end']}"],
            "note": [f"Consensus sequence from ~{sample_id} {haplotype}"]
        }
    )
    features.append(source_feature)

    # Gene feature covering entire sequence
    gene_feature = SeqFeature(
        location=FeatureLocation(0, seq_len, strand=+1),
        type="gene",
        qualifiers={
            "label": [metadata['symbol']],
            "gene": [metadata['symbol']],
            "note": [f"color: #a6acb3"]
        }
    )
    features.append(gene_feature)

    return features

def create_edit_part_features(alignments: Dict[str, List[Hit]], query_lens: Dict[str, int]) -> List[SeqFeature]:
    """
    Create SeqFeature objects for each aligned edit part.

    Args:
        alignments: {query_part_name: [Hit, ...]}
        query_lens: {query_part_name: length}

    Returns:
        List of SeqFeature objects for edit parts
    """
    features = []

    for query_name, hits in alignments.items():
        # Parse part type from query name (e.g., "SYMBOL_exonX_left_HA")
        color = "#ff0000"
        if '_left_HA' in query_name:
            query_name = "5' Homology Arm"
            part_type = 'left_HA'
            color = "#00ffff"
        elif '_right_HA' in query_name:
            query_name = "3' Homology Arm"
            part_type = 'right_HA'
            color = "#00ffff"
        elif '_payload' in query_name:
            query_name = 'Payload'
            part_type = 'payload'
            color = "#75c6a9"
        elif '_HDR_template' in query_name:
            query_name = query_name.replace("_", " ")
            part_type = 'HDR_template'
            color = "#c0c0c0"
        elif '_WT_template' in query_name:
            query_name = query_name.replace("_", " ")
            part_type = 'WT_template'
            color = "#808080"
        else:
            query_name = query_name.replace("_", " ")
            part_type = 'unknown'

        for hit in hits:
            q_start, q_end, t_start, t_end, strand, n_match = hit

            # Calculate coverage for this hit
            coverage = n_match / query_lens[query_name] if query_name in query_lens else 0.0

            # Create feature location
            if strand == '+':
                location = FeatureLocation(t_start, t_end, strand=+1)
            else:
                location = FeatureLocation(t_start, t_end, strand=-1)

            # Create SeqFeature
            feature = SeqFeature(
                location=location,
                type="misc_feature",
                qualifiers={
                    "label": [query_name],
                    "note": [f"color: {color}"],
                    "part_type": [part_type],
                    "coverage": [f"{coverage:.3f}"],
                    "alignment_score": [str(n_match)]
                }
            )

            features.append(feature)

    return features

def annotate_consensus_sequences(faithful_alignments: Dict[str, Dict[str, List[Hit]]],
                                 fasta_path: str,
                                 query_lens: Dict[str, int],
                                 edit_metadata: Dict,
                                 output_prefix: str) -> Tuple[int, int]:
    """
    Annotate haplotype consensus sequences with source, gene, and edit part features.

    Returns: (num_sequences_annotated, num_edit_features_added)
    """
    # Load consensus sequences
    consensus_seqs = {rec.id: rec for rec in SeqIO.parse(fasta_path, "fasta")}

    annotated_records = []
    total_edit_features = 0

    for seq_id, alignments in faithful_alignments.items():
        if seq_id not in consensus_seqs:
            print(f"Warning: Sequence {seq_id} has alignments but no sequence in FASTA")
            continue

        # Get consensus sequence
        original_rec = consensus_seqs[seq_id]

        # Determine haplotype from sequence ID
        if '_hap1' in seq_id:
            haplotype = 'hap1'
        elif '_hap2' in seq_id:
            haplotype = 'hap2'
        else:
            haplotype = 'unknown'

        # Create source and gene features
        source_gene_features = create_source_and_gene_features(
            len(original_rec.seq), edit_metadata, haplotype
        )

        # Create edit part features
        edit_features = create_edit_part_features(alignments, query_lens)

        # Combine all features
        all_features = source_gene_features + edit_features

        # Create annotated SeqRecord
        annotated_rec = SeqRecord(
            seq=original_rec.seq,
            id=original_rec.id,
            name=original_rec.id[:16],  # GenBank name max 16 chars
            description=f"~{sample_id} {edit_metadata['symbol']} {haplotype} consensus {edit_metadata['chromosome']}:{edit_metadata['start']}-{edit_metadata['end']}",
            features=all_features,
            annotations={
                "source": "~{sample_id}",
                "organism": "Homo sapiens",
                "molecule_type": "DNA",
                "topology": "linear",
                "date": datetime.now().strftime("%d-%b-%Y").upper()
            }
        )

        annotated_records.append(annotated_rec)
        total_edit_features += len(edit_features)

    # Write GenBank file
    output_gbk = f"{output_prefix}.gbk"
    SeqIO.write(annotated_records, output_gbk, "genbank")

    return len(annotated_records), total_edit_features

def main():
    parser = argparse.ArgumentParser(
        description="Create GenBank annotations for CRISPR edit parts on haplotype consensus sequences"
    )
    parser.add_argument("paf", help="PAF alignment file")
    parser.add_argument("fasta", help="FASTA file with consensus sequences")
    parser.add_argument("json", help="CRISPR edit JSON file")
    parser.add_argument("-c", "--coverage_threshold", type=float, default=0.8,
                       help="Minimum coverage threshold for faithful alignments")
    parser.add_argument("-o", "--output_prefix", default="edit_annotations",
                       help="Output file prefix")

    args = parser.parse_args()

    print("Loading faithful alignments from PAF...")
    faithful_alignments, query_lens = load_faithful_alignments(args.paf, args.coverage_threshold)
    print(f"Found {len(faithful_alignments)} consensus sequences with faithful alignments")

    print("Loading edit metadata...")
    edit_metadata = load_edit_metadata(args.json)
    print(f"Edit: {edit_metadata['edit_id']} at {edit_metadata['chromosome']}:{edit_metadata['start']}-{edit_metadata['end']}")

    print("Annotating consensus sequences and creating GenBank files...")
    num_seqs, num_features = annotate_consensus_sequences(
        faithful_alignments,
        args.fasta,
        query_lens,
        edit_metadata,
        args.output_prefix
    )

    # Write summary
    summary_file = f"{args.output_prefix}_summary.txt"
    with open(summary_file, 'w') as f:
        f.write(f"CRISPR Edit Consensus Annotation Summary\n")
        f.write(f"=========================================\n")
        f.write(f"Edit: {edit_metadata['edit_id']}\n")
        f.write(f"Gene: {edit_metadata['symbol']}\n")
        f.write(f"Location: {edit_metadata['chromosome']}:{edit_metadata['start']}-{edit_metadata['end']}\n")
        f.write(f"Consensus sequences annotated: {num_seqs}\n")
        f.write(f"Edit part features added: {num_features}\n")
        if num_seqs > 0:
            f.write(f"Average features per sequence: {num_features/num_seqs:.2f}\n")
        f.write(f"Coverage threshold: {args.coverage_threshold}\n")

    print(f"✓ Created GenBank file with {num_seqs} annotated consensus sequences ({num_features} edit part features)")
    print(f"✓ Summary written to {summary_file}")

if __name__ == "__main__":
    main()
EOF

    chmod +x create_genbank_annotations.py

    # Run the annotation script
    python3 create_genbank_annotations.py \
      ~{parts_alignment_paf} \
      ~{consensus_fasta} \
      ~{crispr_edit_json} \
      -c ~{coverage_threshold} \
      -o ~{out_prefix}
  >>>

  output {
    File genbank_annotations = "~{out_prefix}.gbk"
    File annotation_summary = "~{out_prefix}_summary.txt"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87"
    cpu: 2
    memory: "8 GB"
    disk: "20 GB"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
  }
}
