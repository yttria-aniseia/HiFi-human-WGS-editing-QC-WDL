version 1.1

import "../wdl-common/wdl/structs.wdl"
import "crispr_edit_tasks.wdl" as CrisprEditTasks
import "crispr_consensus_tasks.wdl" as CrisprConsensusTasks

workflow crispr_edit_consensus {
  meta {
    description: "Generate haplotype consensus sequences and annotate CRISPR edit parts in GenBank format using sample-specific reference"
  }

  parameter_meta {
    sample_id: "Sample identifier"
    haplotagged_bam: "Haplotagged BAM file from phasing"
    haplotagged_bam_index: "Index for haplotagged BAM"
    ref_fasta: "Reference genome FASTA"
    ref_fasta_index: "Reference genome FASTA index"
    small_variant_vcf: "Sample small variant VCF (phased)"
    small_variant_vcf_index: "Index for small variant VCF"
    sv_vcf: "Sample structural variant VCF (phased)"
    sv_vcf_index: "Index for structural variant VCF"
    crispr_edit_json: "JSON file describing expected CRISPR edit with target region"
    coverage_threshold: "Minimum coverage threshold for faithful alignments (0-1)"
  }

  input {
    String sample_id
    File haplotagged_bam
    File haplotagged_bam_index
    File ref_fasta
    File ref_fasta_index
    File small_variant_vcf
    File small_variant_vcf_index
    File sv_vcf
    File sv_vcf_index
    File crispr_edit_json
    Float coverage_threshold = 0.8
    RuntimeAttributes runtime_attributes
  }

  # Step 1: Create sample-specific reference consensus using bcftools consensus
  call CrisprConsensusTasks.create_sample_reference_consensus {
    input:
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      small_variant_vcf = small_variant_vcf,
      small_variant_vcf_index = small_variant_vcf_index,
      sv_vcf = sv_vcf,
      sv_vcf_index = sv_vcf_index,
      crispr_edit_json = crispr_edit_json,
      sample_id = sample_id,
      runtime_attributes = runtime_attributes
  }

  # Step 2a: Extract reads from edit region
  call CrisprConsensusTasks.extract_region_reads {
    input:
      haplotagged_bam = haplotagged_bam,
      haplotagged_bam_index = haplotagged_bam_index,
      crispr_edit_json = crispr_edit_json,
      sample_id = sample_id,
      runtime_attributes = runtime_attributes
  }

  # Step 2b: Align reads to sample-specific reference using minimap2
  call CrisprConsensusTasks.minimap2_align_to_sample_reference {
    input:
      region_reads_fastq = extract_region_reads.region_reads_fastq,
      sample_ref_consensus_fasta = create_sample_reference_consensus.sample_ref_consensus_fasta,
      sample_id = sample_id,
      runtime_attributes = runtime_attributes
  }

  # Step 2c: Convert SAM to sorted BAM and transfer haplotype tags
  call CrisprConsensusTasks.sam_to_sorted_tagged_bam {
    input:
      sam_file = minimap2_align_to_sample_reference.realigned_sam,
      original_bam = haplotagged_bam,
      original_bam_index = haplotagged_bam_index,
      crispr_edit_json = crispr_edit_json,
      sample_id = sample_id,
      runtime_attributes = runtime_attributes
  }

  # Step 3: Create haplotype consensus sequences from realigned reads
  call CrisprConsensusTasks.create_haplotype_consensus_from_realigned {
    input:
      realigned_bam = sam_to_sorted_tagged_bam.tagged_bam,
      realigned_bam_index = sam_to_sorted_tagged_bam.tagged_bam_index,
      sample_ref_consensus_fasta = create_sample_reference_consensus.sample_ref_consensus_fasta,
      crispr_edit_json = crispr_edit_json,
      sample_id = sample_id,
      runtime_attributes = runtime_attributes
  }

  # Step 4: Create FASTA with edit parts for alignment (including payload_components for detailed annotation)
  call CrisprEditTasks.create_parts_query_fasta {
    input:
      crispr_edit_json = crispr_edit_json,
      include_payload_components = true,
      runtime_attributes = runtime_attributes
  }

  # Step 5: Align edit parts to haplotype consensus sequences
  call CrisprEditTasks.minimap2_remap_parts {
    input:
      filtered_reads_fasta = create_haplotype_consensus_from_realigned.haplotype_consensus_fasta,
      parts_query_fasta = create_parts_query_fasta.parts_query_fasta,
      sample_id = sample_id,
      runtime_attributes = runtime_attributes
  }

  # Step 6: Create annotated GenBank files
  call CrisprConsensusTasks.annotate_consensus_genbank {
    input:
      parts_alignment_paf = minimap2_remap_parts.parts_alignment_paf,
      consensus_fasta = create_haplotype_consensus_from_realigned.haplotype_consensus_fasta,
      crispr_edit_json = crispr_edit_json,
      coverage_threshold = coverage_threshold,
      sample_id = sample_id,
      runtime_attributes = runtime_attributes
  }

  output {
    File sample_ref_consensus_fasta = create_sample_reference_consensus.sample_ref_consensus_fasta
    File realigned_bam = sam_to_sorted_tagged_bam.tagged_bam
    File realigned_bam_index = sam_to_sorted_tagged_bam.tagged_bam_index
    File alignment_stats = sam_to_sorted_tagged_bam.alignment_stats
    File haplotype_consensus_fasta = create_haplotype_consensus_from_realigned.haplotype_consensus_fasta
    File consensus_summary = create_haplotype_consensus_from_realigned.consensus_summary
    File parts_alignment_paf = minimap2_remap_parts.parts_alignment_paf
    File genbank_annotations = annotate_consensus_genbank.genbank_annotations
    File annotation_summary = annotate_consensus_genbank.annotation_summary
  }
}
