version 1.1
import "humanwgs_structs.wdl"
import "wdl-common/wdl/workflows/backend_configuration/backend_configuration.wdl" as BackendConfiguration
import "process_trgt_catalog/process_trgt_catalog.wdl" as ProcessTrgtCatalog
import "upstream/upstream.wdl" as Upstream
import "joint/joint.wdl" as Joint
import "downstream/downstream.wdl" as Downstream
import "wdl-common/wdl/tasks/bcftools.wdl" as Bcftools
import "wdl-common/wdl/tasks/trgt.wdl" as Trgt
import "tertiary/tertiary.wdl" as TertiaryAnalysis
import "wdl-common/wdl/tasks/utilities.wdl" as Utilities
import "edit-qc/bcftools_aux.wdl" as Bcftools_aux
import "edit-qc/bcftools_norm.wdl" as Bcftools_norm
import "edit-qc/truvari_parent_filter.wdl" as TruvariParentFilter
import "somatic_ports/somatic_annotation.wdl" as Somatic_annotation
import "somatic_ports/somatic_calling.wdl" as Somatic_calling

workflow humanwgs_family {
  meta {
    description: "PacBio HiFi human whole genome sequencing pipeline, with joint calling for related samples."
  }

  parameter_meta {
    family: {
      name: "Family struct describing samples, relationships, and unaligned BAM paths"
    }
    phenotypes: {
      name: "Comma-delimited list of HPO codes for phenotypes"
    }
    ref_map_file: {
      name: "TSV containing reference genome file paths; must match backend"
    }
    tertiary_map_file: {
      name: "TSV containing tertiary analysis file paths and thresholds; must match backend"
    }
    somatic_map_file: {
      name: "TSV containing somatic analysis file paths; must match backend"
    }
    max_reads_per_alignment_chunk: {
      name: "Maximum reads per alignment chunk"
    }
    pharmcat_min_coverage: {
      name: "Minimum coverage for PharmCAT"
    }
    glnexus_mem_gb: {
      name: "Override GLnexus memory request (GB)"
    }
    gpu: {
      name: "Use GPU when possible"
    }
    backend: {
      name: "Backend where the workflow will be executed",
      choices: ["GCP", "Azure", "AWS-HealthOmics", "HPC"]
    }
    zones: {
      name: "Zones where compute will take place; required if backend is set to 'GCP'"
    }
    cpuPlatform: {
      help: "Optional minimum CPU platform to use for tasks on GCP"
    }
    gpuType: {
      name: "GPU type to use; required if gpu is set to `true` for cloud backends; must match backend"
    }
    container_registry: {
      name: "Container registry where workflow images are hosted. If left blank, PacBio's public Quay.io registry will be used. Must be set if backend is set to 'AWS-HealthOmics'"
    }
    preemptible: {
      name: "Where possible, run tasks preemptibly"
    }
    debug_version: {
      name: "Debug version for testing purposes"
    }
  }

  input {
    Family family

    String phenotypes = "HP:0000001"

    File ref_map_file
    File? tertiary_map_file
    File somatic_map_file

    Int max_reads_per_alignment_chunk = 500000
    Int pharmcat_min_coverage = 10
    Int? glnexus_mem_gb

    Boolean gpu = false

    # Backend configuration
    String backend
    String? zones
    String? cpuPlatform
    String? gpuType
    String? container_registry

    Boolean preemptible = true

    String? debug_version
  }

  call BackendConfiguration.backend_configuration {
    input:
      backend            = backend,
      zones              = zones,
      cpuPlatform        = cpuPlatform,
      gpuType            = gpuType,
      container_registry = container_registry
  }

  RuntimeAttributes default_runtime_attributes = if preemptible then backend_configuration.spot_runtime_attributes else backend_configuration.on_demand_runtime_attributes

  Map [String, String] ref_map = read_map(ref_map_file)

  call ProcessTrgtCatalog.process_trgt_catalog {
    input:
      trgt_catalog               = ref_map["trgt_tandem_repeat_bed"],  # !FileCoercion
      ref_fasta                  = ref_map["fasta"],                   # !FileCoercion
      ref_index                  = ref_map["fasta_index"],             # !FileCoercion
      default_runtime_attributes = default_runtime_attributes
  }

  Boolean single_sample = length(family.samples) == 1

  Map[String, String] pedigree_sex = {
    "MALE": "1",
    "FEMALE": "2",
    "": "."
  }

  ####################################
  # 1a) UPSTREAM: per‐sample calls  #
  ####################################
  scatter (sample in family.samples) {
    String sample_id = sample.sample_id
    Boolean is_trio_kid = defined(sample.father_id) && defined(sample.mother_id)  # !UnusedDeclaration
    Boolean is_duo_kid = defined(sample.father_id) != defined(sample.mother_id)   # !UnusedDeclaration

    call Upstream.upstream {
      input:
        sample_id                     = sample.sample_id,
        sex                           = sample.sex,
        hifi_reads                    = sample.hifi_reads,
        fail_reads                    = sample.fail_reads,
        ref_map_file                  = ref_map_file,
        trgt_catalog                  = process_trgt_catalog.full_catalog,
        fail_reads_bed                = process_trgt_catalog.include_fail_reads_bed,
        fail_reads_bait_fasta         = process_trgt_catalog.fail_reads_bait_fasta,
        fail_reads_bait_index         = process_trgt_catalog.fail_reads_bait_index,
        max_reads_per_alignment_chunk = max_reads_per_alignment_chunk,
        single_sample                 = single_sample,
        gpu                           = gpu,
        default_runtime_attributes    = default_runtime_attributes
    }

    # write sample metadata similar to pedigree format
    # family_id, sample_id, father_id, mother_id, sex, affected
    Array[String] sample_metadata = [
      family.family_id,
      sample.sample_id,
      select_first([sample.father_id, "."]),
      select_first([sample.mother_id, "."]),
      pedigree_sex[upstream.inferred_sex],
      if sample.affected then "2" else "1"
    ]
  }

    ####################################
    # 1b)         JOINT:                #
    ####################################
    # after your upstream scatter and single_sample Boolean…
  if (!single_sample) {
    call Joint.joint {
      input:
        family_id                  = family.family_id,
        sample_ids                 = sample_id,
        gvcfs                      = upstream.small_variant_gvcf,
        gvcf_indices               = upstream.small_variant_gvcf_index,
        discover_tars              = upstream.discover_tar,
        aligned_bams               = upstream.out_bam,
      aligned_bam_indices        = upstream.out_bam_index,
        ref_map_file               = ref_map_file,
        glnexus_mem_gb             = glnexus_mem_gb,
        default_runtime_attributes = default_runtime_attributes
    }
  }

  scatter (sample_index in range(length(family.samples))) {
    call Bcftools_aux.bcftools_merge_assembly_align as merge_sv_vcfs_align_assembly {
      input:
      vcfs = select_all([upstream.large_sv_filtered_vcf[sample_index], select_all(select_first([joint.split_joint_structural_variant_vcfs,upstream.sv_vcf]))[sample_index]]),
      out_prefix  = "~{sample_id[sample_index]}.merged_structural_variants",
      runtime_attributes = default_runtime_attributes
    }
  }

    ####################################
    # 1c) DOWNSTREAM: per‐sample calls  #
    ####################################
  scatter (sample_index in range(length(family.samples))) {
    call Downstream.downstream {
      input:
        sample_id                  = sample_id[sample_index],
        small_variant_vcf          = select_first([joint.split_joint_small_variant_vcfs, upstream.small_variant_vcf])[sample_index],
        small_variant_vcf_index    = select_first([joint.split_joint_small_variant_vcf_indices, upstream.small_variant_vcf_index])[sample_index],
        sv_vcf                     = merge_sv_vcfs_align_assembly.merged_vcf[sample_index],
        sv_vcf_index               = merge_sv_vcfs_align_assembly.merged_vcf_index[sample_index],
        trgt_vcf                   = upstream.trgt_vcf[sample_index],
        trgt_vcf_index             = upstream.trgt_vcf_index[sample_index],
        trgt_catalog               = process_trgt_catalog.full_catalog,
        aligned_bam                = upstream.out_bam[sample_index],
        aligned_bam_index          = upstream.out_bam_index[sample_index],
        pharmcat_min_coverage      = pharmcat_min_coverage,
        ref_map_file               = ref_map_file,
        default_runtime_attributes = default_runtime_attributes
    }
  }

    #############################################################
    #       MERGING SVs and small variants ACROSS SAMPLES       #
    #############################################################



    if (!single_sample) {
        call Bcftools.bcftools_merge as merge_small_variant_vcfs {
          input:
            vcfs               = downstream.phased_small_variant_vcf,
            vcf_indices        = downstream.phased_small_variant_vcf_index,
            out_prefix         = "~{family.family_id}.joint.~{ref_map['name']}.small_variants.phased",
            runtime_attributes = default_runtime_attributes
        }

        call Bcftools.bcftools_merge as merge_sv_vcfs {
          input:
            vcfs               = downstream.phased_sv_vcf,
            vcf_indices        = downstream.phased_sv_vcf_index,
            out_prefix         = "~{family.family_id}.joint.~{ref_map['name']}.structural_variants.phased",
            runtime_attributes = default_runtime_attributes
        }

        call Trgt.trgt_merge {
          input:
            vcfs               = downstream.phased_trgt_vcf,
            vcf_indices        = downstream.phased_trgt_vcf_index,
            ref_fasta          = ref_map["fasta"],                              # !FileCoercion
            ref_index          = ref_map["fasta_index"],                        # !FileCoercion
            out_prefix         = "~{family.family_id}.merged.~{ref_map['name']}.trgt",
            runtime_attributes = default_runtime_attributes
        }
  }


  #############################################################
  # 3)                  Tertiary Analysis                     #
  #############################################################

  if (defined(tertiary_map_file)) {
    call TertiaryAnalysis.tertiary_analysis {
      input:
      sample_metadata            = sample_metadata,
      phenotypes                 = phenotypes,
      is_trio_kid                = is_trio_kid,
      is_duo_kid                 = is_duo_kid,
      small_variant_vcf          = select_first([merge_small_variant_vcfs.merged_vcf, downstream.phased_small_variant_vcf[0]]),
      small_variant_vcf_index    = select_first([merge_small_variant_vcfs.merged_vcf_index, downstream.phased_small_variant_vcf_index[0]]),
      sv_vcf                     = select_first([merge_sv_vcfs.merged_vcf, downstream.phased_sv_vcf[0]]),
      sv_vcf_index               = select_first([merge_sv_vcfs.merged_vcf_index, downstream.phased_sv_vcf_index[0]]),
      ref_map_file               = ref_map_file,
      tertiary_map_file          = select_first([tertiary_map_file]),
      default_runtime_attributes = default_runtime_attributes
    }
  }

  #############################################################
  #                  SEVERUS & ANNOTATION                     #
  #############################################################

  if (defined(tertiary_map_file) && defined(somatic_map_file)) {
    Map [String, String] somatic_map = read_map(somatic_map_file)

    scatter (sample in family.samples) {
      Array[File] hifi_reads = sample.hifi_reads
    }

    ####################################
    # 3a)      SEVERUS & ANNOTATION    #
    ####################################
    Map[String, File] bam_files_by_id = as_map(zip(sample_id, downstream.merged_haplotagged_bam))
    Map[String, File] bam_index_by_id = as_map(zip(sample_id, downstream.merged_haplotagged_bam_index))
    Map[String, Int] sidx_by_id = as_map(zip(sample_id, range(length(sample_id))))

    call Bcftools_norm.bcftools_norm_split_multiallelic as normalize_small_variants {
      input:
      input_vcf           = select_first([merge_small_variant_vcfs.merged_vcf, downstream.phased_small_variant_vcf]),
      ref_fasta           = ref_map["fasta"],
      ref_fasta_index     = ref_map["fasta_index"],
      out_prefix          = "~{family.family_id}.small_variants.normalized",
      runtime_attributes  = default_runtime_attributes
    }

    call Bcftools_norm.bcftools_norm_split_multiallelic as normalize_sv {
      input:
      input_vcf           = select_first([merge_sv_vcfs.merged_vcf, downstream.phased_sv_vcf]),
      input_vcf_index     = select_first([merge_sv_vcfs.merged_vcf_index, downstream.phased_sv_vcf_index]),
      ref_fasta           = ref_map["fasta"],
      ref_fasta_index     = ref_map["fasta_index"],
      out_prefix          = "~{family.family_id}.sv.normalized",
      runtime_attributes  = default_runtime_attributes
    }

    call Somatic_annotation.vep_annotate as annotateGermline {
      input:
      input_vcf           = normalize_small_variants.normalized_vcf,
      vep_cache           = somatic_map["vep_cache"],                           # !FileCoercion
      ref_fasta           = ref_map["fasta"],
      ref_fasta_index     = ref_map["fasta_index"],
      threads             = 16
    }

    call TruvariParentFilter.truvari_collapse as truvari_collapse_sv {
      input:
      merged_vcf         = normalize_sv.normalized_vcf,
      merged_vcf_index   = normalize_sv.normalized_vcf_index,
      out_prefix         = "~{family.family_id}.sv.parent_filter",
      runtime_attributes = default_runtime_attributes
    }

    call TruvariParentFilter.bcftools_split_for_truvari as bcftools_split_sv {
      input:
      truvari_merge_vcf  = truvari_collapse_sv.truvari_merge_vcf,
      out_prefix         = "~{family.family_id}.sv.parent_filter",
      runtime_attributes = default_runtime_attributes
    }

    call TruvariParentFilter.truvari_consistency as truvari_consistency_sv {
      input:
      split_vcfs         = bcftools_split_sv.split_vcfs,
      split_vcfs_index   = bcftools_split_sv.split_vcfs_index,
      sample_order_file  = bcftools_split_sv.sample_order_file,
      out_prefix         = "~{family.family_id}.sv.parent_filter",
      runtime_attributes = default_runtime_attributes
    }

    call TruvariParentFilter.truvari_collapse as truvari_collapse_small_variants {
      input:
      merged_vcf         = annotateGermline.annotated_vcf,
      merged_vcf_index   = annotateGermline.annotated_vcf_index,
      out_prefix         = "~{family.family_id}.small_variants.parent_filter",
      runtime_attributes = default_runtime_attributes
    }

    call TruvariParentFilter.bcftools_split_for_truvari as bcftools_split_small_variants {
      input:
      truvari_merge_vcf  = truvari_collapse_small_variants.truvari_merge_vcf,
      out_prefix         = "~{family.family_id}.small_variants.parent_filter",
      runtime_attributes = default_runtime_attributes
    }

    call TruvariParentFilter.truvari_consistency as truvari_consistency_small_variants {
      input:
      split_vcfs         = bcftools_split_small_variants.split_vcfs,
      split_vcfs_index   = bcftools_split_small_variants.split_vcfs_index,
      sample_order_file  = bcftools_split_small_variants.sample_order_file,
      out_prefix         = "~{family.family_id}.small_variants.parent_filter",
      runtime_attributes = default_runtime_attributes
    }

    scatter (sample_index in range(length(family.samples))) {

      # Sex matching, if sample is male, we use male parent, if female we use female parent.
      # if no sex given for sample, use father sample by default.
      # Generally Prefer father_id if present, else use mother_id if present, else null unless sample is Female
      String? parent_id =
      if (family.samples[sample_index].sex == "MALE") then
      if (defined(family.samples[sample_index].father_id)) then select_first([family.samples[sample_index].father_id])
      else if (defined(family.samples[sample_index].mother_id)) then select_first([family.samples[sample_index].mother_id])
      else None
      else if (family.samples[sample_index].sex == "FEMALE") then
      if (defined(family.samples[sample_index].mother_id)) then select_first([family.samples[sample_index].mother_id])
      else if (defined(family.samples[sample_index].father_id)) then select_first([family.samples[sample_index].father_id])
      else None
      else
      if (defined(family.samples[sample_index].father_id)) then select_first([family.samples[sample_index].father_id])
      else if (defined(family.samples[sample_index].mother_id)) then select_first([family.samples[sample_index].mother_id])
      else None

      File? parental_bam_find       = if (defined(parent_id)) then bam_files_by_id[select_first([parent_id])] else None
      File? parental_bam_index_find = if (defined(parent_id)) then bam_index_by_id[select_first([parent_id])] else None
      # call Somatic_calling.Severus_sv as phased_severus {
        #   input:
        #     wt_bam             = downstream.merged_haplotagged_bam[sample_index],
        #     wt_bam_index       = downstream.merged_haplotagged_bam_index[sample_index],
        #     parental_bam       = parental_bam_find,
        #     parental_bam_index = parental_bam_index_find,
        #     trf_bed            = somatic_map["trf_bed"],                    # !FileCoercion
        #     phased_vcf         = downstream.phased_small_variant_vcf[sample_index],
        #     threads            = 32,
        #     min_supp_reads     = 3,
        #     PON_tsv            = somatic_map["severus_pon_tsv"]             # !FileCoercion
      # }
      # call Somatic_calling.tabix_vcf as tabix_vcf {
        #  input:
        #    vcf        = select_first([phased_severus.output_vcf]),
        #    contig_bed = somatic_map["ref_bed"],                                   # !FileCoercion
        #    threads    = 2
      # }
      # call Somatic_calling.svpack_filter_annotated as svpack_filter_annotated {
        #  input:
        #    sv_vcf                 = tabix_vcf.output_vcf,
        #    population_vcfs        = [somatic_map["control_vcf"]],                 # !FileCoercion
        #    population_vcf_indices = [somatic_map["control_vcf_index"]],           # !FileCoercion
        #    gff                    = somatic_map["ref_gff"]                        # !FileCoercion
      # }
      # call Somatic_calling.recover_mate_bnd as recover_mate_bnd {
        #  input:
        #    sv_vcf_original    = select_first([phased_severus.output_vcf]),
        #    sv_svpack_filtered = svpack_filter_annotated.output_vcf
      # }
      # call Somatic_annotation.annotsv as annotateSeverusSVfiltered {
        #   input:
        #     sv_vcf        = select_first([recover_mate_bnd.output_vcf]),
        #     sv_vcf_index  = select_first([recover_mate_bnd.output_vcf_index]),
        #     annotsv_cache = somatic_map["annotsv_cache"],                           # !FileCoercion
        #     threads       = 2
      # }


      # parent filtering
      Int parent_index = if defined(parent_id) then sidx_by_id[select_first([parent_id])] else 999

      call TruvariParentFilter.filter_parent_variants as parent_filter_small_variants {
        input:
        sample_vcf              = bcftools_split_small_variants.split_vcfs[sample_index],
        sample_vcf_index        = bcftools_split_small_variants.split_vcfs_index[sample_index],
        truvari_consistency_tsv = truvari_consistency_small_variants.truvari_consistency_tsv,
        sample_index            = sample_index,
        parent_index            = parent_index,
        runtime_attributes      = default_runtime_attributes
      }

      call TruvariParentFilter.filter_parent_variants as parent_filter_sv {
        input:
        sample_vcf              = bcftools_split_sv.split_vcfs[sample_index],
        sample_vcf_index        = bcftools_split_sv.split_vcfs_index[sample_index],
        truvari_consistency_tsv = truvari_consistency_sv.truvari_consistency_tsv,
        sample_index            = sample_index,
        parent_index            = parent_index,
        runtime_attributes      = default_runtime_attributes
      }

      # AnnotSV after split because we need its tsv conversion
      call Somatic_annotation.annotsv as annotate_parent_filter_sv {
        input:
        sv_vcf        = parent_filter_sv.filtered_vcf,
        sv_vcf_index  = parent_filter_sv.filtered_vcf_index,
        annotsv_cache = somatic_map["annotsv_cache"],                             # !FileCoercion
        threads       = 2
      }

      call Somatic_annotation.prioritize_small_variants as prioritizeSomatic {
        input:
        vep_annotated_vcf   = parent_filter_small_variants.filtered_vcf,
        threads             = 2,
        sample              = family.samples[sample_index].sample_id
      }
      call Somatic_annotation.prioritize_sv_intogen as prioritize_sv {
        input:
        annotSV_tsv   = annotate_parent_filter_sv.annotated_tsv,
        threads       = 2
      }
      # call Somatic_annotation.prioritize_sv_intogen as prioritize_Severus {
        #   input:
        #     annotSV_tsv   = annotateSeverusSVfiltered.annotated_tsv,
        #     threads       = 2
      # }
    }
  }
  
    Array[Array[String]] stats = [
      flatten([['sample_id'], sample_id]),
      flatten([['read_count'], downstream.stat_read_count]),
      flatten([['read_length_mean'], downstream.stat_read_length_mean]),
      flatten([['read_length_median'], downstream.stat_read_length_median]),
      flatten([['read_length_n50'], downstream.stat_read_length_n50]),
      flatten([['read_quality_mean'], downstream.stat_read_quality_mean]),
      flatten([['read_quality_median'], downstream.stat_read_quality_median]),
      flatten([['mapped_read_count'], downstream.stat_mapped_read_count]),
      flatten([['mapped_read_percent'], downstream.stat_mapped_read_percent]),
      flatten([['gap_compressed_identity_mean'], downstream.stat_gap_compressed_identity_mean]),
      flatten([['gap_compressed_identity_median'], downstream.stat_gap_compressed_identity_median]),
      flatten([['depth_mean'], upstream.stat_depth_mean]),
      flatten([['inferred_sex'], upstream.inferred_sex]),
      flatten([['stat_phased_basepairs'], downstream.stat_phased_basepairs]),
      flatten([['phase_block_ng50'], downstream.stat_phase_block_ng50]),
      flatten([['cpg_combined_count'], downstream.stat_combined_cpg_count]),
      flatten([['cpg_hap1_count'], downstream.stat_hap1_cpg_count]),
      flatten([['cpg_hap2_count'], downstream.stat_hap2_cpg_count]),
      flatten([['methbat_methylated_count'], downstream.stat_methbat_methylated_count]),
      flatten([['methbat_unmethylated_count'], downstream.stat_methbat_unmethylated_count]),
      flatten([['methbat_asm_count'], downstream.stat_methbat_asm_count]),
      flatten([['SNV_count'], downstream.stat_SNV_count]),
      flatten([['TSTV_ratio'], downstream.stat_TSTV_ratio]),
      flatten([['HETHOM_ratio'], downstream.stat_HETHOM_ratio]),
      flatten([['INDEL_count'], downstream.stat_INDEL_count]),
      flatten([['sv_DUP_count'], downstream.stat_sv_DUP_count]),
      flatten([['sv_DEL_count'], downstream.stat_sv_DEL_count]),
      flatten([['sv_INS_count'], downstream.stat_sv_INS_count]),
      flatten([['sv_INV_count'], downstream.stat_sv_INV_count]),
      flatten([['sv_SWAP_count'], downstream.stat_sv_SWAP_count]),
      flatten([['sv_BND_count'], downstream.stat_sv_BND_count]),
      flatten([['trgt_genotyped_count'], upstream.stat_trgt_genotyped_count]),
      flatten([['trgt_uncalled_count'], upstream.stat_trgt_uncalled_count])
  ]

  call Utilities.consolidate_stats {
    input:
      id                 = family.family_id,
      stats              = stats,
      msg_array          = flatten([process_trgt_catalog.msg, flatten(upstream.msg)]),
      runtime_attributes = default_runtime_attributes
  }
  
  output {
    # to maintain order of samples
    Array[String] sample_ids = sample_id
    File  stats_file         = consolidate_stats.output_tsv
    File  msg_file           = consolidate_stats.messages

    ## WT OUTPUTS
    # bam stats
    Array[File]   bam_statistics                      = downstream.bam_statistics
    Array[File]   read_length_plot                    = downstream.read_length_plot
    Array[File?]  read_quality_plot                   = downstream.read_quality_plot
    Array[File]   mapq_distribution_plot              = downstream.mapq_distribution_plot
    Array[File]   mg_distribution_plot                = downstream.mg_distribution_plot
    Array[String] stat_read_count                     = downstream.stat_read_count
    Array[String] stat_read_length_mean               = downstream.stat_read_length_mean
    Array[String] stat_read_length_median             = downstream.stat_read_length_median
    Array[String] stat_read_length_n50                = downstream.stat_read_length_n50
    Array[String] stat_read_quality_mean              = downstream.stat_read_quality_mean
    Array[String] stat_read_quality_median            = downstream.stat_read_quality_median
    Array[String] stat_mapped_read_count              = downstream.stat_mapped_read_count
    Array[String] stat_mapped_read_percent            = downstream.stat_mapped_read_percent
    Array[String] stat_gap_compressed_identity_mean   = downstream.stat_gap_compressed_identity_mean
    Array[String] stat_gap_compressed_identity_median = downstream.stat_gap_compressed_identity_median

    # merged, haplotagged alignments
    Array[File]   merged_haplotagged_bam       = downstream.merged_haplotagged_bam
    Array[File]   merged_haplotagged_bam_index = downstream.merged_haplotagged_bam_index

    # mosdepth outputs
    Array[File]   mosdepth_summary                 = upstream.mosdepth_summary
    Array[File]   mosdepth_region_bed              = upstream.mosdepth_region_bed
    Array[File]   mosdepth_region_bed_index        = upstream.mosdepth_region_bed_index
    Array[File]   mosdepth_depth_distribution_plot = upstream.mosdepth_depth_distribution_plot
    Array[String] stat_depth_mean                  = upstream.stat_depth_mean
    Array[String] inferred_sex                     = upstream.inferred_sex

    # mosdepth hg002 outputs
    Array[File]   mosdepth_hg002_summary                 = upstream.mosdepth_hg002_summary
    Array[File]   mosdepth_hg002_region_bed              = upstream.mosdepth_hg002_region_bed
    Array[File]   mosdepth_hg002_region_bed_index        = upstream.mosdepth_hg002_region_bed_index
    Array[File]   mosdepth_hg002_depth_distribution_plot = upstream.mosdepth_hg002_depth_distribution_plot
    Array[String] stat_mean_depth_hg002                  = upstream.stat_mean_depth_hg002
    Array[String] inferred_sex_hg002                     = upstream.inferred_sex_hg002

    # phasing stats
    Array[File]   phase_stats           = downstream.phase_stats
    Array[File]   phase_blocks          = downstream.phase_blocks
    Array[File]   phase_haplotags       = downstream.phase_haplotags
    Array[String] stat_phased_basepairs = downstream.stat_phased_basepairs
    Array[String] stat_phase_block_ng50 = downstream.stat_phase_block_ng50

    # methylation outputs and profile
    Array[File?]  cpg_combined_bed                 = downstream.cpg_combined_bed
    Array[File?]  cpg_combined_bed_index           = downstream.cpg_combined_bed_index
    Array[File?]  cpg_hap1_bed                     = downstream.cpg_hap1_bed
    Array[File?]  cpg_hap1_bed_index               = downstream.cpg_hap1_bed_index
    Array[File?]  cpg_hap2_bed                     = downstream.cpg_hap2_bed
    Array[File?]  cpg_hap2_bed_index               = downstream.cpg_hap2_bed_index
    Array[File?]  cpg_combined_bw                  = downstream.cpg_combined_bw
    Array[File?]  cpg_hap1_bw                      = downstream.cpg_hap1_bw
    Array[File?]  cpg_hap2_bw                      = downstream.cpg_hap2_bw
    Array[String] stat_cpg_hap1_count              = downstream.stat_hap1_cpg_count
    Array[String] stat_cpg_hap2_count              = downstream.stat_hap2_cpg_count
    Array[String] stat_cpg_combined_count          = downstream.stat_combined_cpg_count
    Array[File?]  methbat_profile                  = downstream.methbat_profile
    Array[String] stat_methbat_methylated_count   = downstream.stat_methbat_methylated_count
    Array[String] stat_methbat_unmethylated_count = downstream.stat_methbat_unmethylated_count
    Array[String] stat_methbat_asm_count          = downstream.stat_methbat_asm_count

    # sv outputs
    Array[File] phased_sv_vcf                 = downstream.phased_sv_vcf
    Array[File] phased_sv_vcf_index           = downstream.phased_sv_vcf_index
    File sv_supporting_reads                  = select_first([joint.sv_supporting_reads, upstream.sv_supporting_reads[0]])
    Array[File] sv_copynum_bedgraph           = select_first([joint.sv_copynum_bedgraph, select_all(upstream.sv_copynum_bedgraph)])
    Array[File] sv_depth_bw                   = select_first([joint.sv_depth_bw, select_all(upstream.sv_depth_bw)])
    Array[File] sv_gc_bias_corrected_depth_bw = select_first([joint.sv_gc_bias_corrected_depth_bw, select_all(upstream.sv_gc_bias_corrected_depth_bw)])
    Array[File] sv_maf_bw                     = select_first([joint.sv_maf_bw, select_all(upstream.sv_maf_bw)])
    Array[File] sv_copynum_summary            = select_first([joint.sv_copynum_summary, select_all(upstream.sv_copynum_summary)])

    # sv stats
    Array[String] stat_sv_DUP_count  = downstream.stat_sv_DUP_count
    Array[String] stat_sv_DEL_count  = downstream.stat_sv_DEL_count
    Array[String] stat_sv_INS_count  = downstream.stat_sv_INS_count
    Array[String] stat_sv_INV_count  = downstream.stat_sv_INV_count
    Array[String] stat_sv_SWAP_count = downstream.stat_sv_SWAP_count
    Array[String] stat_sv_BND_count  = downstream.stat_sv_BND_count

    # small variant outputs
    Array[File] phased_small_variant_vcf       = downstream.phased_small_variant_vcf
    Array[File] phased_small_variant_vcf_index = downstream.phased_small_variant_vcf_index
    Array[File] small_variant_gvcf             = upstream.small_variant_gvcf
    Array[File] small_variant_gvcf_index       = upstream.small_variant_gvcf_index

    # small variant stats
    Array[File]   small_variant_stats             = downstream.small_variant_stats
    Array[File]   bcftools_roh_out                = downstream.bcftools_roh_out
    Array[File]   bcftools_roh_bed                = downstream.bcftools_roh_bed
    Array[String] stat_small_variant_SNV_count    = downstream.stat_SNV_count
    Array[String] stat_small_variant_INDEL_count  = downstream.stat_INDEL_count
    Array[String] stat_small_variant_TSTV_ratio   = downstream.stat_TSTV_ratio
    Array[String] stat_small_variant_HETHOM_ratio = downstream.stat_HETHOM_ratio
    Array[File]   snv_distribution_plot           = downstream.snv_distribution_plot
    Array[File]   indel_distribution_plot         = downstream.indel_distribution_plot

    # trgt outputs
    Array[File]   phased_trgt_vcf           = downstream.phased_trgt_vcf
    Array[File]   phased_trgt_vcf_index     = downstream.phased_trgt_vcf_index
    Array[File]   trgt_spanning_reads       = upstream.trgt_spanning_reads
    Array[File]   trgt_spanning_reads_index = upstream.trgt_spanning_reads_index
    Array[File]   trgt_coverage_dropouts    = downstream.trgt_coverage_dropouts
    Array[String] stat_trgt_genotyped_count = upstream.stat_trgt_genotyped_count
    Array[String] stat_trgt_uncalled_count  = upstream.stat_trgt_uncalled_count

    # paraphase outputs
    Array[File?] paraphase_output_json         = upstream.paraphase_output_json
    Array[File?] paraphase_realigned_bam       = upstream.paraphase_realigned_bam
    Array[File?] paraphase_realigned_bam_index = upstream.paraphase_realigned_bam_index
    Array[File?] paraphase_vcfs                = upstream.paraphase_vcfs

    # per sample mitorsaw outputs
    Array[File] mitorsaw_vcf       = upstream.mitorsaw_vcf
    Array[File] mitorsaw_vcf_index = upstream.mitorsaw_vcf_index
    Array[File] mitorsaw_hap_stats = upstream.mitorsaw_hap_stats

    # PGx outputs
    Array[File]  pbstarphase_json        = downstream.pbstarphase_json
    Array[File?] pharmcat_match_json     = downstream.pharmcat_match_json
    Array[File?] pharmcat_phenotype_json = downstream.pharmcat_phenotype_json
    Array[File?] pharmcat_report_html    = downstream.pharmcat_report_html
    Array[File?] pharmcat_report_json    = downstream.pharmcat_report_json

    # Assembly outputs
    # hifiasm upstream outputs
    Array[File] asm_hap1 = upstream.asm_1
    Array[File] asm_hap2 = upstream.asm_2

    # pav outputs
    Array[File?] pav_vcf = upstream.pav_vcf
    Array[File?] pav_vcf_index = upstream.pav_vcf_index

    # liftover outputs
    Array[File?] asm_lifted_vcf = upstream.lifted_vcf
    Array[File?] asm_reject_vcf = upstream.reject_vcf

    # LARGE SV FILTER OUTPUTS
    Array[File?] pav_large_sv_filtered_vcf = upstream.large_sv_filtered_vcf
    Array[File?] pav_large_sv_filtered_vcf_index = upstream.large_sv_filtered_vcf_index
    # confusing with joint merge
    #Array[File] merged_assembly_aligned_sv_vcf = select_first([merge_sv_vcfs_align_assembly.merged_vcf])
    #Array[File] merged_assembly_aligned_sv_vcf_index = select_first([merge_sv_vcfs_align_assembly.merged_vcf_index])

    # joint call outputs
    File? joint_small_variants_vcf       = merge_small_variant_vcfs.merged_vcf
    File? joint_small_variants_vcf_index = merge_small_variant_vcfs.merged_vcf_index
    File? joint_sv_vcf                   = merge_sv_vcfs.merged_vcf
    File? joint_sv_vcf_index             = merge_sv_vcfs.merged_vcf_index
    File? joint_trgt_vcf                 = trgt_merge.merged_vcf
    File? joint_trgt_vcf_index           = trgt_merge.merged_vcf_index

    # Parent filtering outputs
    # sample consistency table
    File? truvari_small_variant_consistency_tsv        = truvari_consistency_small_variants.truvari_consistency_tsv
    File? truvari_small_variant_consistency_report_txt = truvari_consistency_small_variants.truvari_consistency_report_txt
    File? truvari_sv_consistency_tsv                   = truvari_consistency_sv.truvari_consistency_tsv
    File? truvari_sv_consistency_report_txt            = truvari_consistency_sv.truvari_consistency_report_txt

    # truvari-collapsed variants
    Array[File]? truvari_small_variant_vcf            = bcftools_split_small_variants.split_vcfs
    Array[File]? truvari_small_variant_vcf_index      = bcftools_split_small_variants.split_vcfs_index
    Array[File]? truvari_sv_vcf                       = bcftools_split_sv.split_vcfs
    Array[File]? truvari_sv_vcf_index                 = bcftools_split_sv.split_vcfs_index

    Array[File]? parent_filtered_small_variant_vcf = parent_filter_small_variants.filtered_vcf
    Array[File]? parent_filtered_sv_vcf            = parent_filter_sv.filtered_vcf

    # Somatic SV calling
    # Array[File] Severus_somatic_vcf                           = select_all(select_first([phased_severus.output_vcf]))
    # Array[File] Severus_all_vcf                               = select_all(select_first([phased_severus.output_all_vcf]))
    # Array[File] Severus_breakpoint_cluster                    = select_all(select_first([phased_severus.output_breakpoint_clusters]))
    # Array[File] Severus_breakpoint_cluster_all                = select_all(select_first([phased_severus.output_breakpoint_clusters_all]))
    # Array[File] Severus_cluster_plots                         = select_all(select_first([phased_severus.output_somatic_sv_plots]))

    # Array[File] Severus_tabix_vcf                             = select_first([tabix_vcf.output_vcf])
    # Array[File] Severus_filtered_vcf                          = select_first([recover_mate_bnd.output_vcf])
    # Array[File] Severus_filtered_vcf_index                    = select_first([recover_mate_bnd.output_vcf_index])

    # annotation analysis outputs
    File? merged_vep_annotated_vcf                            = annotateGermline.annotated_vcf

    # SV annotations only come in tsv form.
    Array[File]? sv_annotated_tsv                             = annotate_parent_filter_sv.annotated_tsv
    Array[File]? sv_annotated_ranked_tsv                      = select_first([prioritize_sv.annotSV_ranked_tsv])
    Array[File]? sv_annotated_concerning_tsv                  = select_first([prioritize_sv.annotSV_concerning_tsv])
    Array[File]? sv_annotated_ccg_tsv                         = select_first([prioritize_sv.annotSV_ccg_tsv])
    Array[File]? sv_annotated_ccg_ranked_tsv                  = select_first([prioritize_sv.annotSV_ccg_ranked_tsv])

    Array[File]? small_variant_annotated_tsv                  = prioritizeSomatic.vep_annotated_tsv
    Array[File]? small_variant_annotated_ranked_tsv           = prioritizeSomatic.vep_annotated_ranked_tsv
    Array[File]? small_variant_annotated_filtered_tsv         = prioritizeSomatic.vep_annotated_filtered_tsv

    # Array[File] Severus_annotated_tsv                         = select_first([annotateSeverusSVfiltered.annotated_tsv])
    # Array[File] Severus_annotated_ranked_tsv                  = select_first([prioritize_Severus.annotSV_ranked_tsv])
    # Array[File] Severus_annotated_ccg_tsv                     = select_first([prioritize_Severus.annotSV_ccg_tsv])
    # Array[File] Severus_annotated_ccg_ranked_tsv              = select_first([prioritize_Severus.annotSV_ccg_ranked_tsv])
    # tertiary analysis outputs
    File? tertiary_small_variant_filtered_vcf           = tertiary_analysis.small_variant_filtered_vcf
    File? tertiary_small_variant_filtered_vcf_index     = tertiary_analysis.small_variant_filtered_vcf_index
    File? tertiary_small_variant_filtered_tsv           = tertiary_analysis.small_variant_filtered_tsv
    File? tertiary_small_variant_compound_het_vcf       = tertiary_analysis.small_variant_compound_het_vcf
    File? tertiary_small_variant_compound_het_vcf_index = tertiary_analysis.small_variant_compound_het_vcf_index
    File? tertiary_small_variant_compound_het_tsv       = tertiary_analysis.small_variant_compound_het_tsv
    File? tertiary_sv_filtered_vcf                      = tertiary_analysis.sv_filtered_vcf
    File? tertiary_sv_filtered_vcf_index                = tertiary_analysis.sv_filtered_vcf_index
    File? tertiary_sv_filtered_tsv                      = tertiary_analysis.sv_filtered_tsv

    # qc messages
    Array[String] msg = flatten(
      [
        process_trgt_catalog.msg,
        flatten(upstream.msg)
      ]
    )

    # workflow metadata
    String workflow_name    = "humanwgs_family"
    String workflow_version = "v3.1.0" + if defined(debug_version) then "~{"-" + debug_version}" else ""
  }
}
