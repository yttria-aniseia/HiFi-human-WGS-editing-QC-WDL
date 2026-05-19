version 1.1

import "../wdl-common/wdl/structs.wdl"
import "crispr_edit_tasks.wdl" as CrisprEditTasks
import "crispr_consensus_tasks.wdl" as CrisprConsensusTasks

workflow crispr_edit_consensus {
  meta {
    description: "Annotate MSA-based group consensus sequences with CRISPR edit parts in GenBank format"
  }

  parameter_meta {
    sample_id: "Sample identifier"
    consensus_fasta: "Per-group consensus FASTA from build_group_consensus (crispr_edit_qc output)"
    crispr_edit_json: "JSON file describing expected CRISPR edit"
    coverage_threshold: "Minimum coverage threshold for faithful alignments (0-1)"
  }

  input {
    String sample_id
    File consensus_fasta
    File crispr_edit_json
    Float coverage_threshold = 0.8
    RuntimeAttributes runtime_attributes
  }

  # Build parts query FASTA (including payload components for detailed annotation)
  call CrisprEditTasks.create_parts_query_fasta {
    input:
      crispr_edit_json          = crispr_edit_json,
      include_payload_components = true,
      runtime_attributes        = runtime_attributes
  }

  # Align edit parts to group consensus sequences
  call CrisprEditTasks.blastn_remap_parts as blastn_remap_consensus_parts {
    input:
      filtered_reads_fasta = consensus_fasta,
      parts_query_fasta    = create_parts_query_fasta.parts_query_fasta,
      sample_id            = sample_id,
      runtime_attributes   = runtime_attributes
  }

  # Annotate consensus sequences and emit GenBank
  call CrisprConsensusTasks.annotate_consensus_genbank {
    input:
      parts_alignment_tsv = blastn_remap_consensus_parts.parts_alignment_tsv,
      consensus_fasta     = consensus_fasta,
      crispr_edit_json    = crispr_edit_json,
      coverage_threshold  = coverage_threshold,
      sample_id           = sample_id,
      runtime_attributes  = runtime_attributes
  }

  output {
    File parts_alignment_tsv = blastn_remap_consensus_parts.parts_alignment_tsv
    File genbank_annotations = annotate_consensus_genbank.genbank_annotations
    File annotation_summary  = annotate_consensus_genbank.annotation_summary
  }
}
