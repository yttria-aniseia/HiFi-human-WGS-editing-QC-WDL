version 1.1

import "../wdl-common/wdl/structs.wdl"
import "crispr_edit_tasks.wdl" as CrisprEditTasks

workflow crispr_edit_qc {
  meta {
    description: "CRISPR edit QC workflow to identify and analyze edit genotypes from HiFi reads"
  }

  parameter_meta {
    sample_id: "Sample identifier"
    fasta_reads: "FASTA converted reads from upstream workflow"
    crispr_edit_json: "JSON file describing expected CRISPR edit"
    min_identity: "Minimum alignment identity for initial filtering (0-1)"
    coverage_threshold: "Minimum coverage threshold for faithful alignments (0-1)"
    threads: "Number of threads for minimap2"
  }

  input {
    String sample_id
    File fasta_reads
    File crispr_edit_json
    Float min_identity = 0.80
    Float coverage_threshold = 0.8
    Int threads = 4
    RuntimeAttributes runtime_attributes
  }

  call CrisprEditTasks.create_parts_query_fasta {
    input:
    crispr_edit_json = crispr_edit_json,
    runtime_attributes = runtime_attributes
  }

  call CrisprEditTasks.minimap2_align_reads {
    input:
    fasta_reads = fasta_reads,
    parts_query_fasta = create_parts_query_fasta.parts_query_fasta,
    min_identity = min_identity,
    threads = threads,
    sample_id = sample_id,
    runtime_attributes = runtime_attributes
  }

  call CrisprEditTasks.samtools_filter_reads {
    input:
    fasta_reads = fasta_reads,
    matched_read_names = minimap2_align_reads.matched_read_names,
    sample_id = sample_id,
    runtime_attributes = runtime_attributes
  }

  call CrisprEditTasks.minimap2_remap_parts {
    input:
    filtered_reads_fasta = samtools_filter_reads.filtered_reads_fasta,
    parts_query_fasta = create_parts_query_fasta.parts_query_fasta,
    sample_id = sample_id,
    runtime_attributes = runtime_attributes
  }

  call CrisprEditTasks.parse_edit_parts {
    input:
    parts_alignment_paf = minimap2_remap_parts.parts_alignment_paf,
    coverage_threshold = coverage_threshold,
    sample_id = sample_id,
    runtime_attributes = runtime_attributes
  }


  output {
    File? parts_query_fasta = create_parts_query_fasta.parts_query_fasta
    File? filtered_reads_fasta = samtools_filter_reads.filtered_reads_fasta
    File? parts_alignment_paf = minimap2_remap_parts.parts_alignment_paf
    File? summary_table = parse_edit_parts.summary_table
    File? faithful_table = parse_edit_parts.faithful_table
    File? wide_table = parse_edit_parts.wide_table
    File? wide_faithful_table = parse_edit_parts.wide_faithful_table
  }
}