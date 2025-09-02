version 1.0

import "../wdl-common/wdl/structs.wdl"

task bcftools_merge {
  meta {
    description: "Merge multiple sample VCFs into a single joint VCF."
  }

  parameter_meta {
    vcfs: {
      name: "VCFs"
    }
    vcf_indices: {
      name: "VCF indices"
    }
    out_prefix: {
      name: "Output VCF name"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    merged_vcf: {
      name: "Merged VCF"
    }
    merged_vcf_index: {
      name: "Merged VCF index"
    }
  }

  input {
    Array[File] vcfs
    Array[File] vcf_indices

    String out_prefix

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 2
  Int mem_gb    = 4
  Int disk_size = ceil(size(vcfs, "GB") * 2 + 20)

  command <<<
    set -euo pipefail

    bcftools --version

    bcftools merge \
      ~{if threads > 1 then "--threads " + (threads - 1) else ""} \
      --output-type z \
      --output ~{out_prefix}.vcf.gz \
      ~{sep=" " vcfs}
    bcftools index --tbi ~{out_prefix}.vcf.gz
  >>>

  output {
    File merged_vcf       = "~{out_prefix}.vcf.gz"
    File merged_vcf_index = "~{out_prefix}.vcf.gz.tbi"
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

task bcftools_merge_assembly_align {
  meta {
    description: "Merge multiple assembly and aligned VCFs into a single joint VCF."
  }

  parameter_meta {
    vcfs: {
      name: "VCFs"
    }
    out_prefix: {
      name: "Output VCF name"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    merged_vcf: {
      name: "Merged VCF"
    }  
    merged_vcf_index: {
      name: "Merged VCF index"
    }
  }

  input {
    Array[File] vcfs

    String out_prefix

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 2
  Int mem_gb    = 4
  Int disk_size = ceil(size(vcfs, "GB") * 2 + 20)

  command <<<
    set -euo pipefail
    for vcf in ~{sep=' ' vcfs}; do
        cp "$vcf" ./
        # Regenerate index to ensure compatibility
        tabix -p vcf "$(basename "$vcf")"
    done

    bcftools --version
    bcftools concat \
      --allow-overlaps \
      --remove-duplicates \
      ~{if threads > 1 then "--threads " + (threads - 1) else ""} \
      --output-type z \
      --output ~{out_prefix}.vcf.gz \
      *.vcf.gz

    #bcftools merge \
    #  ~{if threads > 1 then "--threads " + (threads - 1) else ""} \
    #  --output-type z \
    #  --output ~{out_prefix}.vcf.gz \
    #  *.vcf.gz
    bcftools index --tbi ~{out_prefix}.vcf.gz
  >>>

  output {
    File merged_vcf       = "~{out_prefix}.vcf.gz"
    File merged_vcf_index = "~{out_prefix}.vcf.gz.tbi"
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