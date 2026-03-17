version 1.1

struct Sample {
  String sample_id

  String? sex
  Boolean affected

  Array[File] hifi_reads
  Array[File]? fail_reads

  String? father_id
  String? mother_id
  File? expected_edit

  # Pre-computed assembly outputs; when provided, hifiasm is skipped
  File? precomputed_asm_hap1  # hifiasm hap1 GFA
  File? precomputed_asm_hap2  # hifiasm hap2 GFA
}

struct Family {
  String family_id
  Array[Sample] samples
}
