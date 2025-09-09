version 1.0

import "../wdl-common/wdl/structs.wdl"
import "../wdl-common/wdl/workflows/pbmm2/pbmm2.wdl" as Pbmm2
import "../wdl-common/wdl/tasks/sawfish.wdl" as Sawfish
import "../wdl-common/wdl/workflows/deepvariant/deepvariant.wdl" as DeepVariant
import "../wdl-common/wdl/tasks/samtools.wdl" as Samtools
import "../wdl-common/wdl/tasks/mosdepth.wdl" as Mosdepth
import "../wdl-common/wdl/tasks/trgt.wdl" as Trgt
import "../wdl-common/wdl/tasks/paraphase.wdl" as Paraphase
import "../edit-qc/samtools_aux.wdl" as Samtools_aux
import "../edit-qc/mosdepth_himem.wdl" as Mosdepth_himem
import "../assembly/assembly.wdl" as Assembly
import "../wdl-common/wdl/tasks/mitorsaw.wdl" as Mitorsaw

workflow upstream {
  meta {
    description: "Given a set of HiFi reads for a human sample, run steps upstream of phasing."
  }

  parameter_meta {
    sample_id: {
      name: "Sample ID"
    }
    sex: {
      name: "Sample sex",
      choices: ["MALE", "FEMALE"]
    }
    hifi_reads: {
      name: "HiFi reads (BAMs)"
    }
    ref_map_file: {
      name: "TSV containing reference genome information"
    }
    max_reads_per_alignment_chunk: {
      name: "Maximum reads per alignment chunk"
    }
    single_sample: {
      name: "Single sample workflow"
    }
    gpu: {
      name: "Use GPU for DeepVariant"
}
    default_runtime_attributes: {
      name: "Runtime attribute structure"
    }
  }

  input {
    String sample_id
    String? sex
    Array[File] hifi_reads

    File ref_map_file

    Int max_reads_per_alignment_chunk

    Boolean single_sample = false

    Boolean gpu

    RuntimeAttributes default_runtime_attributes
  }

  Map[String, String] ref_map = read_map(ref_map_file)

  scatter (hifi_read_bam in hifi_reads) {
    call Pbmm2.pbmm2 as pbmm2 {
      input:
        sample_id                  = sample_id,
        bam                        = hifi_read_bam,
        max_reads_per_chunk        = max_reads_per_alignment_chunk,
        ref_fasta                  = ref_map["fasta"],              # !FileCoercion
        ref_index                  = ref_map["fasta_index"],        # !FileCoercion
        ref_name                   = ref_map["name"],
        default_runtime_attributes = default_runtime_attributes
    }
    call Pbmm2.pbmm2_align_wgs as pbmm2_align_hg002 {
      input:
        sample_id          = sample_id,
        bam                = hifi_read_bam,
        ref_fasta          = ref_map["hg002_fasta"],       # !FileCoercion
        ref_index          = ref_map["hg002_fasta_index"], # !FileCoercion
        ref_name           = ref_map["hg002_name"],
        runtime_attributes = default_runtime_attributes
    }
  }

  # merge aligned bams if there are multiple
  if (length(flatten(pbmm2.aligned_bams)) > 1) {
    call Samtools.samtools_merge {
      input:
        bams               = flatten(pbmm2.aligned_bams),
        out_prefix         = "~{sample_id}.~{ref_map['name']}",
        runtime_attributes = default_runtime_attributes
    }
    call Samtools.samtools_merge as samtools_merge_hg002 {
      input:
        bams               = pbmm2_align_hg002.aligned_bam,
        out_prefix         = "~{sample_id}.~{ref_map['hg002_name']}",
        runtime_attributes = default_runtime_attributes
    }
  }

  # select the merged bam if it exists, otherwise select the first (only) aligned bam
  File aligned_bam_data  = select_first([samtools_merge.merged_bam, flatten(pbmm2.aligned_bams)[0]])
  File aligned_bam_index = select_first([samtools_merge.merged_bam_index, flatten(pbmm2.aligned_bam_indices)[0]])

  File aligned_bam_data_hg002  = select_first([samtools_merge_hg002.merged_bam, pbmm2_align_hg002.aligned_bam[0]])
  File aligned_bam_index_hg002 = select_first([samtools_merge_hg002.merged_bam_index, pbmm2_align_hg002.aligned_bam_index[0]])


  call Mosdepth_himem.mosdepth as mosdepth_hg002 {
    input:
      sample_id          = sample_id,
      ref_name           = ref_map["hg002_name"],
      aligned_bam        = aligned_bam_data_hg002,
      aligned_bam_index  = aligned_bam_index_hg002,
      infer_sex          = true,
      runtime_attributes = default_runtime_attributes
  }

  call Mosdepth.mosdepth as mosdepth {
    input:
      sample_id          = sample_id,
      ref_name           = ref_map["name"],
      aligned_bam        = aligned_bam_data,
      aligned_bam_index  = aligned_bam_index,
      infer_sex          = true,
      runtime_attributes = default_runtime_attributes
  }

  String qc_sex = 
    if (defined(sex) && (mosdepth.inferred_sex != sex)) 
    then "~{sample_id}: Reported sex ~{sex} does not match inferred sex ~{mosdepth.inferred_sex}."
    else ""

  call DeepVariant.deepvariant {
    input:
      sample_id                  = sample_id,
      aligned_bams               = [aligned_bam_data],
      aligned_bam_indices        = [aligned_bam_index],
      ref_fasta                  = ref_map["fasta"],       # !FileCoercion
      ref_index                  = ref_map["fasta_index"], # !FileCoercion
      ref_name                   = ref_map["name"],
      gpu                        = gpu,
      default_runtime_attributes = default_runtime_attributes
  }

  call Sawfish.sawfish_discover {
    input:
      sample_id               = sample_id,
      sex                     = mosdepth.inferred_sex,
      aligned_bam             = aligned_bam_data,
      aligned_bam_index       = aligned_bam_index,
      ref_fasta               = ref_map["fasta"],                       # !FileCoercion
      ref_index               = ref_map["fasta_index"],                 # !FileCoercion
      exclude_bed             = ref_map["sawfish_exclude_bed"],         # !FileCoercion
      exclude_bed_index       = ref_map["sawfish_exclude_bed_index"],   # !FileCoercion
      expected_male_bed       = ref_map["sawfish_expected_bed_male"],   # !FileCoercion
      expected_female_bed     = ref_map["sawfish_expected_bed_female"], # !FileCoercion
      small_variant_vcf       = deepvariant.vcf,
      small_variant_vcf_index = deepvariant.vcf_index,
      out_prefix              = "~{sample_id}",
      runtime_attributes      = default_runtime_attributes
  }

  call Trgt.trgt {
    input:
      sample_id          = sample_id,
      sex                = mosdepth.inferred_sex,
      aligned_bam        = aligned_bam_data,
      aligned_bam_index  = aligned_bam_index,
      ref_fasta          = ref_map["fasta"],                  # !FileCoercion
      ref_index          = ref_map["fasta_index"],            # !FileCoercion
      trgt_bed           = ref_map["trgt_tandem_repeat_bed"], # !FileCoercion
      out_prefix         = "~{sample_id}.~{ref_map['name']}",
      runtime_attributes = default_runtime_attributes
  }

  call Paraphase.paraphase {
    input:
      aligned_bam        = aligned_bam_data,
      aligned_bam_index  = aligned_bam_index,
      ref_fasta          = ref_map["fasta"],         # !FileCoercion
      ref_index          = ref_map["fasta_index"],   # !FileCoercion
      sample_id          = sample_id,
      runtime_attributes = default_runtime_attributes
  }

  call Mitorsaw.mitorsaw {
    input:
      aligned_bam        = aligned_bam_data,
      aligned_bam_index  = aligned_bam_index,
      ref_fasta          = ref_map["fasta"],                  # !FileCoercion
      ref_index          = ref_map["fasta_index"],            # !FileCoercion
      out_prefix         = "~{sample_id}.~{ref_map['name']}",
      runtime_attributes = default_runtime_attributes
  }

  if (single_sample) {
    call Sawfish.sawfish_call {
      input:
        sample_ids          = [sample_id],
        discover_tars       = [sawfish_discover.discover_tar],
        aligned_bams        = [aligned_bam_data],
        aligned_bam_indices = [aligned_bam_index],
        ref_fasta           = ref_map["fasta"],                                      # !FileCoercion
        ref_index           = ref_map["fasta_index"],                                # !FileCoercion
        out_prefix          = "~{sample_id}.~{ref_map['name']}.structural_variants",
        runtime_attributes  = default_runtime_attributes
    }

    File copynum_bedgraph_output           = sawfish_call.copynum_bedgraph[0]
    File depth_bw_output                   = sawfish_call.depth_bw[0]
    File gc_bias_corrected_depth_bw_output = sawfish_call.gc_bias_corrected_depth_bw[0]
    File maf_bw_output                     = sawfish_call.maf_bw[0]
  }
  if (length(hifi_reads) > 1) {
    call Samtools_aux.samtools_cat as samtools_cat {
      input:
        bams               = hifi_reads,
        out_prefix         = "~{sample_id}.~{ref_map['name']}",
        runtime_attributes = default_runtime_attributes
    }
  }
  File assembly_input_bam = select_first([samtools_cat.merged_bam, hifi_reads[0]])
  call Assembly.assembly {
      input:
        hifi_read_bam   = assembly_input_bam,
        sample_id   = "~{sample_id}",
        ref_map     = ref_map,
        default_runtime_attributes = default_runtime_attributes
  }


  output {
    # alignments
    File out_bam       = aligned_bam_data
    File out_bam_index = aligned_bam_index

    #hg002 alingments
    File out_bam_hg002       = aligned_bam_data_hg002
    File out_bam_hg002_index = aligned_bam_index_hg002

    # mosdepth outputs
    File   mosdepth_summary                 = mosdepth.summary
    File   mosdepth_region_bed              = mosdepth.region_bed
    File   mosdepth_region_bed_index        = mosdepth.region_bed_index
    File   mosdepth_depth_distribution_plot = mosdepth.depth_distribution_plot
    String inferred_sex                     = mosdepth.inferred_sex
    String stat_mean_depth                  = mosdepth.stat_mean_depth

    # hg002 mosdepth_hg002 outputs
    File   mosdepth_hg002_summary                 = mosdepth_hg002.summary
    File   mosdepth_hg002_region_bed              = mosdepth_hg002.region_bed
    File   mosdepth_hg002_region_bed_index        = mosdepth_hg002.region_bed_index
    File   mosdepth_hg002_depth_distribution_plot = mosdepth_hg002.depth_distribution_plot
    String inferred_sex_hg002                     = mosdepth_hg002.inferred_sex
    String stat_mean_depth_hg002                  = mosdepth_hg002.stat_mean_depth

    # per sample sv signatures
    File discover_tar = sawfish_discover.discover_tar

    # sawfish outputs for single sample
    File? sv_vcf                        = sawfish_call.vcf
    File? sv_vcf_index                  = sawfish_call.vcf_index
    File? sv_supporting_reads           = sawfish_call.supporting_reads
    File? sv_copynum_bedgraph           = copynum_bedgraph_output
    File? sv_depth_bw                   = depth_bw_output
    File? sv_gc_bias_corrected_depth_bw = gc_bias_corrected_depth_bw_output
    File? sv_maf_bw                     = maf_bw_output

    # small variant outputs
    File small_variant_vcf        = deepvariant.vcf
    File small_variant_vcf_index  = deepvariant.vcf_index
    File small_variant_gvcf       = deepvariant.gvcf
    File small_variant_gvcf_index = deepvariant.gvcf_index

    # trgt outputs
    File   trgt_vcf                  = trgt.vcf
    File   trgt_vcf_index            = trgt.vcf_index
    File   trgt_spanning_reads       = trgt.bam
    File   trgt_spanning_reads_index = trgt.bam_index
    String stat_trgt_genotyped_count = trgt.stat_genotyped_count
    String stat_trgt_uncalled_count  = trgt.stat_uncalled_count

    # paraphase outputs
    File? paraphase_output_json         = paraphase.out_json
    File? paraphase_realigned_bam       = paraphase.bam
    File? paraphase_realigned_bam_index = paraphase.bam_index
    File? paraphase_vcfs                = paraphase.vcfs_tar

    # assembly outputs
    #catted bam
    File? cat_bam = samtools_cat.merged_bam
    File? cat_bam_index = samtools_cat.merged_bam_index

    #bam to fasta outputs
    File fasta_output = assembly.fasta_output

    # hifiasm assembly outputs
    File asm_1 = assembly.asm_1
    File asm_2 = assembly.asm_2

    # quast report output
    File quast_transposed_report = assembly.transposed_report
    File quast_icarus_html = assembly.icarus_html
    File quast_report_html = assembly.report_html
    File quast_report_pdf = assembly.report_pdf
    File quast_report_txt = assembly.report_txt

    # pav outputs
    File? pav_vcf = assembly.pav_vcf
    File? pav_vcf_index = assembly.pav_vcf_index

    #liftover outputs
    File? lifted_vcf = assembly.lifted_vcf
    File? reject_vcf = assembly.reject_vcf

    #LARGE SV FILTER OUTPUTS
    File? large_sv_filtered_vcf = assembly.large_sv_filtered_vcf
    File? large_sv_filtered_vcf_index = assembly.large_sv_filtered_vcf_index

    # per sample mitorsaw outputs
    File mitorsaw_vcf       = mitorsaw.vcf
    File mitorsaw_vcf_index = mitorsaw.vcf_index
    File mitorsaw_hap_stats = mitorsaw.hap_stats

    # qc messages
    Array[String] msg = flatten(
      [
        flatten(pbmm2.msg),
        [qc_sex],
        trgt.msg,
        sawfish_discover.msg
      ]
    )
  }
}
