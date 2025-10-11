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
}

struct Family {
  String family_id
  Array[Sample] samples
}
