version 1.1
import "../wdl-common/wdl/structs.wdl"
import "assembly_tasks.wdl" as assembly_tasks

workflow assembly {
  meta {
    description: "Given a set of HiFi reads for a human sample, run steps assembly pipeline."
  }

  parameter_meta {
    hifi_reads: {
      name: "HiFi reads (BAMs)"
    }
    sample_id: {
      name: "Sample ID information"
    }
    ref_map_file: {
      name: "TSV containing reference genome information"
    }
  }

  input {
    Array[File] hifi_read_bams

    String sample_id

    Map[String, String] ref_map

    # Pre-computed hifiasm GFA outputs; when provided, hifiasm is skipped
    File? precomputed_asm_hap1
    File? precomputed_asm_hap2

    RuntimeAttributes default_runtime_attributes
  }

  call assembly_tasks.samtools_bam_to_fasta as bam_to_fasta {
    input:
      input_bams         = hifi_read_bams,
      threads            = 32,
      runtime_attributes = default_runtime_attributes
  }

  # Skip hifiasm if pre-computed GFAs are provided
  if (!defined(precomputed_asm_hap1)) {
    call assembly_tasks.hifiasm_assembly as hifiasm_assembly {
      input:
        input_fasta        = bam_to_fasta.input_fasta,
        threads            = 32,
        runtime_attributes = default_runtime_attributes
    }
  }

  File asm_hap1_gfa = select_first([precomputed_asm_hap1, hifiasm_assembly.input_1_asm])
  File asm_hap2_gfa = select_first([precomputed_asm_hap2, hifiasm_assembly.input_2_asm])

  call assembly_tasks.gfa_to_fa as gfa_to_fa {
    input:
      input_hap1_gfa     = asm_hap1_gfa,
      input_hap2_gfa     = asm_hap2_gfa,
      runtime_attributes = default_runtime_attributes
  }
  # call assembly_tasks.quast as quast {
  #   input:
  #     input_fa_1 = gfa_to_fa.fasta_hap1,
  #     input_fa_2 = gfa_to_fa.fasta_hap2,
  #     input_fa = bam_to_fasta.input_fasta,
  #     ref = ref_map["hg002_fasta"],
  #     runtime_attributes = default_runtime_attributes
  # }
  call assembly_tasks.pav as pav {
    input:
      ref            = ref_map["hg002_fasta"],
      ref_index      = ref_map["hg002_fasta_index"],
      ref_bgzf_index = ref_map["hg002_fasta_bgzf_index"],
      input_fa_1     = gfa_to_fa.fasta_hap1,
      input_fa_2     = gfa_to_fa.fasta_hap2,
      sample_name    = sample_id,
      runtime_attributes = default_runtime_attributes
  }
  if (defined(pav.pav_vcf)) {
    call assembly_tasks.liftover as liftover {
      input:
      vcf                = select_first([pav.pav_vcf]),
      chain              = ref_map["hg002_chain"],
      ref                = ref_map["fasta"],
      runtime_attributes = default_runtime_attributes
    }
    call assembly_tasks.pav_sv_filter_size as pav_sv_filter_size {
      input:
        vcf                = liftover.lifted_vcf,
        runtime_attributes = default_runtime_attributes
    }
  }

  output {
    # bam to fasta outputs
    File fasta_output = bam_to_fasta.input_fasta

    # hifiasm assembly outputs
    File asm_1 = asm_hap1_gfa
    File asm_2 = asm_hap2_gfa
    # File transposed_report = quast.transposed_report
    # File icarus_html       = quast.icarus_html
    # File report_html       = quast.report_html
    # File report_pdf        = quast.report_pdf
    # File report_txt        = quast.report_txt

    # pav outputs
    File? pav_vcf = pav.pav_vcf
    File? pav_vcf_index = pav.pav_vcf_index

    #liftover outputs
    File? lifted_vcf = liftover.lifted_vcf
    File? reject_vcf = liftover.reject_vcf

    #LARGE SV FILTER OUTPUTS
    File? large_sv_filtered_vcf = pav_sv_filter_size.large_sv_filtered_vcf
    File? large_sv_filtered_vcf_index = pav_sv_filter_size.large_sv_filtered_vcf_index
  }
}
