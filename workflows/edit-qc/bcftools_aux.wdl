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

task count_vcf_variants {
  meta {
    description: "Count the number of variant records in a VCF for a specific sample, excluding homozygous-ref and missing genotypes."
  }

  parameter_meta {
    vcf: {
      name: "Input VCF"
    }
    sample_id: {
      name: "Sample to count variants for (if omitted, counts all records)"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    variant_count: {
      name: "Number of variant records"
    }
  }

  input {
    File vcf
    String? sample_id
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(vcf, "GB") + 5)

  command <<<
    set -euo pipefail
    if [ -n "~{default="" sample_id}" ]; then
      # Count records where this sample has a non-ref, non-missing genotype
      bcftools view -H -s "~{sample_id}" -i 'GT!="mis" && GT!="ref"' ~{vcf} | wc -l
    else
      bcftools view -H ~{vcf} | wc -l
    fi
  >>>

  output {
    String variant_count = read_string(stdout())
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87"
    cpu: 2
    memory: "8 GB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
  }
}

task count_tsv_rows {
  meta {
    description: "Count data rows in a TSV file (excluding header)."
  }

  parameter_meta {
    tsv: {
      name: "Input TSV"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    row_count: {
      name: "Number of data rows"
    }
  }

  input {
    File tsv
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(tsv, "GB") + 5)

  command <<<
    set -euo pipefail
    tail -n +2 ~{tsv} | wc -l
  >>>

  output {
    String row_count = read_string(stdout())
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87"
    cpu: 2
    memory: "8 GB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
  }
}

task bcftools_qual_filter {
  meta {
    description: "Filter variants to FILTER=PASS. Small variants and sawfish SVs additionally require GQ>=20 and, if AD is present, VAF>0.1. PAV assembly SVs (identified by presence of INFO/COV_PROP) skip the GQ filter and instead require all haplotype COV_PROP values >0.1."
  }

  parameter_meta {
    vcf: { name: "Input VCF" }
    out_prefix: { name: "Output file prefix" }
    runtime_attributes: { name: "Runtime attribute structure" }
    filtered_vcf: { name: "Filtered VCF" }
    filtered_vcf_index: { name: "Filtered VCF index" }
  }

  input {
    File vcf
    String out_prefix
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(vcf, "GB") * 2 + 10)

  command <<<
    set -euo pipefail

    bcftools --version

    # PAV assembly SVs carry INFO/COV_PROP typed as String in their VCF header.
    # Reheader to Float (or add the tag as Float if absent) so bcftools can
    # evaluate MIN(INFO/COV_PROP) numerically.  Records without the tag remain
    # absent (compared as ".") in either case.
    if bcftools view -h ~{vcf} | grep -q '##INFO=<ID=COV_PROP'; then
      bcftools view -h ~{vcf} | \
        sed 's|##INFO=<ID=COV_PROP,Number=\.,Type=String|##INFO=<ID=COV_PROP,Number=.,Type=Float|' \
        > fixed_header.txt
    else
      bcftools view -h ~{vcf} | grep -v '^#CHROM' > fixed_header.txt
      echo '##INFO=<ID=COV_PROP,Number=.,Type=Float,Description="Assembly haplotype coverage proportion (PAV)">' \
        >> fixed_header.txt
      bcftools view -h ~{vcf} | grep '^#CHROM' >> fixed_header.txt
    fi
    bcftools reheader -h fixed_header.txt -o reheaded.vcf.gz ~{vcf}
    bcftools index --tbi reheaded.vcf.gz

    # Unified filter:
    #   Non-PAV (INFO/COV_PROP absent): require GQ>=20; if AD present, VAF>0.1.
    #   PAV (INFO/COV_PROP present):    require MIN(COV_PROP)>0.1; no GQ check.
    bcftools view \
      --output-type z \
      -i '(FILTER="PASS" || FILTER=".") && (
        (INFO/COV_PROP="." && GQ>=20 && (AD[0:0]="." || (AD[0:0]+AD[0:1])=0 || AD[0:1]/(AD[0:0]+AD[0:1]) > 0.1)) ||
        (INFO/COV_PROP!="." && MIN(INFO/COV_PROP) > 0.1)
      )' \
      reheaded.vcf.gz \
      -o ~{out_prefix}.vcf.gz
    bcftools index --tbi ~{out_prefix}.vcf.gz
  >>>

  output {
    File filtered_vcf       = "~{out_prefix}.vcf.gz"
    File filtered_vcf_index = "~{out_prefix}.vcf.gz.tbi"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87"
    cpu: 4
    memory: "8 GB"
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
    sample_id: {
      name: "Sample ID to normalize VCF sample names to before merging"
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
    String sample_id

    String out_prefix

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 2
  Int mem_gb    = 4
  Int disk_size = ceil(size(vcfs, "GB") * 2 + 20)

  command <<<
    set -euo pipefail
    echo "~{sample_id}" > sample_name.txt
    for vcf in ~{sep=' ' vcfs}; do
        base="$(basename $vcf)"
        bcftools reheader -s sample_name.txt "$vcf" -o "$base"
        bcftools index --tbi "$base"
    done

    bcftools --version
    bcftools concat \
      --allow-overlaps \
      --remove-duplicates \
      ~{if threads > 1 then "--threads " + (threads - 1) else ""} \
      *.vcf.gz | \
    bcftools filter -e 'SVCLAIM="D"' -O z -Wtbi -o "~{out_prefix}.vcf.gz"
    # ^ filter out sawfish depth-evidenced CNV records

    #bcftools merge \
    #  ~{if threads > 1 then "--threads " + (threads - 1) else ""} \
    #  --output-type z \
    #  --output ~{out_prefix}.vcf.gz \
    #  *.vcf.gz
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
