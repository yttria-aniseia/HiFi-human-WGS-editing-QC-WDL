version 1.0

import "../wdl-common/wdl/structs.wdl"

task samtools_cat {
  meta {
    description: "Cat multiple unsorted unaligned BAM files into a single BAM file."
  }

  parameter_meta {
    bams: {
      name: "BAMs"
    }
    out_prefix: {
      name: "Output BAM name"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    merged_bam: {
      name: "Merged BAM"
    }
    merged_bam_index: {
      name: "Merged BAM index"
    }
  }

  input {
    Array[File] bams

    String out_prefix

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 8
  Int mem_gb    = 4
  Int disk_size = ceil(size(bams, "GB") * 2 + 20)

  command <<<
    set -euo pipefail

    samtools --version

    samtools cat \
      ~{if threads > 1 then "--threads " + (threads - 1) else ""} \
      -o ~{out_prefix}.bam \
      ~{sep=' ' bams}

    samtools index ~{out_prefix}.bam
  >>>

  output {
    File merged_bam       = "~{out_prefix}.bam"
    File merged_bam_index = "~{out_prefix}.bam.bai"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87"
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
