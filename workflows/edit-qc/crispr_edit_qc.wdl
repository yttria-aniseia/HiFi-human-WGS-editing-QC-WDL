version 1.1

import "../wdl-common/wdl/structs.wdl"
import "crispr_edit_tasks.wdl" as CrisprEditTasks

workflow crispr_edit_qc {
  meta {
    description: "CRISPR edit QC workflow to identify and analyze edit genotypes from HiFi reads"
  }

  parameter_meta {
    sample_id: "Sample identifier"
    fasta_reads: "FASTA converted reads from upstream workflow (whole genome; mutually exclusive with pre_filtered_reads_fasta)"
    pre_filtered_reads_fasta: "Pre-filtered FASTA already restricted to the edit region (skips minimap2 filtering when provided)"
    haplotagged_bam: "Optional whole-genome aligned BAM. When provided, reads with high-confidence primary alignments >offtarget_window_bp from the edit site are excluded before edit-part search; unmapped and low-MAPQ reads are retained."
    haplotagged_bam_index: "Index for haplotagged_bam"
    offtarget_window_bp: "Distance (bp) from edit target interval within which on-target primary alignments are kept"
    offtarget_mapq_threshold: "Primary alignments with MAPQ <= this value are kept as low-confidence (not excluded as off-target)"
    crispr_edit_json: "JSON file describing expected CRISPR edit"
    min_identity: "Minimum alignment identity for initial filtering (0-1)"
    coverage_threshold: "Minimum coverage threshold for faithful alignments (0-1)"
    threads: "Number of threads for minimap2 (initial read filtering step)"
    clip_padding: "Flanking bases to retain around part-hit span when clipping reads for MSA"
    gap_merge_threshold: "Max gap difference (bp) between two reads to merge into the same structural group"
    wt_gap_threshold:    "Max gap (bp) between left_HA and right_HA to call WT_candidate (covers cut-site sequence)"
    min_flank_bp:        "Minimum flanking bases beyond part hits required to call HDR or WT (vs inconclusive)"
    min_part_hit_bp:     "Minimum query-span (bp) for a blastn hit to be counted; filters repetitive short matches"
  }

  input {
    String sample_id
    File? fasta_reads
    File? pre_filtered_reads_fasta
    File? haplotagged_bam
    File? haplotagged_bam_index
    File crispr_edit_json
    Float min_identity = 0.80
    Float coverage_threshold = 0.8
    Int offtarget_window_bp = 500000
    Int offtarget_mapq_threshold = 20
    Int threads = 8
    Int clip_padding = 100
    Int gap_merge_threshold = 5
    Int wt_gap_threshold    = 20
    Int min_flank_bp = 50
    Int min_part_hit_bp = 100
    RuntimeAttributes runtime_attributes
  }

  # ── Build parts query FASTA from edit JSON ────────────────────────────────
  call CrisprEditTasks.create_parts_query_fasta {
    input:
    crispr_edit_json = crispr_edit_json,
    runtime_attributes = runtime_attributes
  }

  # ── Optional off-target read exclusion using whole-genome BAM ─────────────
  if (!defined(pre_filtered_reads_fasta) && defined(haplotagged_bam)) {
    call CrisprEditTasks.filter_offtarget_reads {
      input:
      haplotagged_bam       = select_first([haplotagged_bam]),
      haplotagged_bam_index = select_first([haplotagged_bam_index]),
      crispr_edit_json      = crispr_edit_json,
      sample_id             = sample_id,
      window_bp             = offtarget_window_bp,
      mapq_threshold        = offtarget_mapq_threshold,
      runtime_attributes    = runtime_attributes
    }
  }

  # ── Initial read filter: minimap2 (skipped when pre_filtered_reads_fasta is supplied) ─
  if (!defined(pre_filtered_reads_fasta)) {
    File minimap2_input_fasta = select_first([filter_offtarget_reads.filtered_reads_fasta,
                                              fasta_reads])
    call CrisprEditTasks.minimap2_align_reads {
      input:
      fasta_reads = minimap2_input_fasta,
      parts_query_fasta = create_parts_query_fasta.parts_query_fasta,
      min_identity = min_identity,
      threads = threads,
      sample_id = sample_id,
      runtime_attributes = runtime_attributes
    }

    call CrisprEditTasks.extract_reads {
      input:
      fasta_reads = minimap2_input_fasta,
      matched_read_names = minimap2_align_reads.matched_read_names,
      sample_id = sample_id,
      runtime_attributes = runtime_attributes
    }
  }

  File actual_filtered_fasta = select_first([pre_filtered_reads_fasta,
                                             extract_reads.filtered_reads_fasta])

  # ── Pass-1 blastn: map parts to filtered reads ────────────────────────────
  call CrisprEditTasks.blastn_remap_parts as blastn_pass1 {
    input:
    filtered_reads_fasta = actual_filtered_fasta,
    parts_query_fasta = create_parts_query_fasta.parts_query_fasta,
    sample_id = sample_id + "_pass1",
    runtime_attributes = runtime_attributes
  }

  # ── Clip reads to union of part-hit spans ─────────────────────────────────
  call CrisprEditTasks.clip_reads_to_parts {
    input:
    parts_alignment_tsv = blastn_pass1.parts_alignment_tsv,
    filtered_reads_fasta = actual_filtered_fasta,
    sample_id = sample_id,
    clip_padding = clip_padding,
    min_part_hit_bp = min_part_hit_bp,
    runtime_attributes = runtime_attributes
  }

  # ── Structural grouping: signature + prototype clustering ─────────────────
  call CrisprEditTasks.structural_grouping {
    input:
    parts_alignment_tsv = blastn_pass1.parts_alignment_tsv,
    sample_id = sample_id,
    gap_merge_threshold = gap_merge_threshold,
    min_part_hit_bp = min_part_hit_bp,
    runtime_attributes = runtime_attributes
  }

  # ── Per-group MSA and error-corrected consensus (three steps) ────────────
  call CrisprEditTasks.prepare_group_fastas {
    input:
    clipped_reads_fasta = clip_reads_to_parts.clipped_reads_fasta,
    groups_tsv = structural_grouping.groups_tsv,
    sample_id = sample_id,
    runtime_attributes = runtime_attributes
  }

  call CrisprEditTasks.run_group_msa {
    input:
    group_inputs_tarball = prepare_group_fastas.group_inputs_tarball,
    sample_id = sample_id,
    runtime_attributes = runtime_attributes
  }

  call CrisprEditTasks.build_group_consensus {
    input:
    msa_tarball = run_group_msa.msa_tarball,
    singletons_fasta = prepare_group_fastas.singletons_fasta,
    groups_tsv = structural_grouping.groups_tsv,
    sample_id = sample_id,
    runtime_attributes = runtime_attributes
  }

  # ── Wide per-part coverage table (backwards-compatible) ──────────────────
  call CrisprEditTasks.make_wide_coverage_table {
    input:
    parts_alignment_tsv = blastn_pass1.parts_alignment_tsv,
    coverage_threshold = coverage_threshold,
    min_part_hit_bp = min_part_hit_bp,
    sample_id = sample_id,
    runtime_attributes = runtime_attributes
  }

  # ── Final categorization ──────────────────────────────────────────────────
  call CrisprEditTasks.categorize_reads {
    input:
    filtered_reads_fasta = actual_filtered_fasta,
    groups_tsv = structural_grouping.groups_tsv,
    subhaps_tsv = build_group_consensus.subhaps_tsv,
    pass2_alignment_tsv = blastn_pass1.parts_alignment_tsv,
    sample_id = sample_id,
    gap_merge_threshold = gap_merge_threshold,
    wt_gap_threshold    = wt_gap_threshold,
    min_flank_bp = min_flank_bp,
    min_part_hit_bp = min_part_hit_bp,
    runtime_attributes = runtime_attributes
  }

  output {
    # Intermediate / diagnostic
    File? parts_query_fasta      = create_parts_query_fasta.parts_query_fasta
    File? filtered_reads_fasta   = extract_reads.filtered_reads_fasta
    File? filtered_reads_fastq   = extract_reads.filtered_reads_fastq
    File? offtarget_filter_counts_tsv = filter_offtarget_reads.counts_tsv
    File? parts_alignment_tsv    = blastn_pass1.parts_alignment_tsv
    File? clipped_reads_fasta    = clip_reads_to_parts.clipped_reads_fasta

    # Per-part coverage table (wide faithful format, backwards-compatible)
    File? wide_faithful_table    = make_wide_coverage_table.wide_faithful_table

    # MSA outputs
    File? msa_tarball            = build_group_consensus.all_msa_tarball
    File? consensus_fasta        = build_group_consensus.consensus_fasta

    # Primary result
    File? categories_tsv         = categorize_reads.categories_tsv
  }
}
