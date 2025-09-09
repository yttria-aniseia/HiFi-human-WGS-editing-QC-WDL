version 1.0

import "../wdl-common/wdl/structs.wdl"

task pbsv_discover {
  meta {
    description: "Discover structural variant signatures with `pbsv discover`."
  }

  parameter_meta {
    aligned_bam: {
      name: "Aligned BAM"
    }
    aligned_bam_index: {
      name: "Aligned BAM index"
    }
    trf_bed: {
      name: "Tandem repeat BED used to normalize representation of variation within tandem repeats"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    svsig: {
      name: "Structural variant signature file"
    }
  }

  input {
    File aligned_bam
    File aligned_bam_index

    File? trf_bed

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 2
  Int mem_gb    = 16
  Int disk_size = ceil((size(aligned_bam, "GB") + size(trf_bed, "GB")) * 2 + 20)

  String out_prefix = basename(aligned_bam, ".bam")

  command <<<
    set -euo pipefail

    pbsv --version

    pbsv discover \
      --log-level INFO \
      --hifi \
      ~{"--tandem-repeats " + trf_bed} \
      ~{aligned_bam} \
      ~{out_prefix}.svsig.gz
  >>>

  output {
    File svsig = "~{out_prefix}.svsig.gz"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pbsv@sha256:3a8529853c1e214809dcdaacac0079de70d0c037b41b43bb8ba7c3fc5f783e26"
    cpu: threads
    memory: mem_gb + " GB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
  }
}
