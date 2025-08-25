version 1.1

import "../wdl-common/wdl/structs.wdl"

task truvari_collapse_consistency {
  meta {
    description: "Use truvari collapse and consistency to identify parent variants for filtering"
  }

  parameter_meta {
    merged_vcf: {
      name: "Merged family VCF with all samples"
    }
    merged_vcf_index: {
      name: "Merged family VCF index"
    }
    sample_ids: {
      name: "Array of sample IDs in order"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    truvari_collapsed_vcf: {
      name: "Truvari collapsed VCF"
    }
    truvari_consistency_tsv: {
      name: "Truvari consistency results TSV"
    }
  }

  input {
    File merged_vcf
    File merged_vcf_index
    Array[String] sample_ids
    String out_prefix = "family_merged"
    RuntimeAttributes runtime_attributes
  }

  Int threads = 4
  Int mem_gb = 8
  Int disk_size = ceil(size(merged_vcf, "GB") * 4 + 50)

  command <<<
    set -euxo pipefail

    # Check truvari version
    truvari version

    # Step 1: truvari collapse to merge variants
    truvari collapse \
      -i ~{merged_vcf} \
      -o ~{out_prefix}_truvari_merge.vcf \
      -c ~{out_prefix}_truvari_collapsed.vcf.gz

    #bgzip ~{out_prefix}_truvari_collapsed.vcf
    #tabix -p vcf ~{out_prefix}_truvari_collapsed.vcf.gz

    # Step 2: bcftools split to create per-sample VCFs, excluding uncalled and reference genotypes
    mkdir -p truvari_merge_split
    bcftools +split -W ~{out_prefix}_truvari_merge.vcf \
      -Oz \
      -o truvari_merge_split \
      -e 'GT="mis" || GT="ref"'

    # Step 3: need to pass exact input sample order to consistency
    bcftools query -l ~{out_prefix}_truvari_merge.vcf > sample_order.txt

    SPLIT_VCFS=()
    while IFS= read -r sample; do
      vcf_file="truvari_merge_split/${sample}.vcf.gz"
      if [[ -f "$vcf_file" ]]; then
        SPLIT_VCFS+=("$vcf_file")
      fi
    done < sample_order.txt

    # Step 4: truvari consistency to identify consistent variants across samples
    truvari consistency \
      -d \
      -o truvari_consistency.tsv \
      "${SPLIT_VCFS[@]}"

    # sort by position
    (head -n1 "truvari_consistency.tsv" && tail -n+2 "truvari_consistency.tsv" | sort -k1,1 -k2,2g) > "~{out_prefix}_truvari_consistency.tsv"
  >>>

  output {
    File truvari_collapsed_vcf = "~{out_prefix}_truvari_collapsed.vcf.gz"
    File truvari_collapsed_vcf_index = "~{out_prefix}_truvari_collapsed.vcf.gz.tbi"
    Array[File] split_vcfs = suffix(".vcf.gz", prefix("truvari_merge_split/", read_lines("sample_order.txt")))
    Array[File] split_vcfs_index = suffix(".vcf.gz.csi", prefix("truvari_merge_split/", read_lines("sample_order.txt")))
    File truvari_consistency_tsv = "~{out_prefix}_truvari_consistency.tsv"
    File sample_order_file = "sample_order.txt"
  }

  runtime {
    docker: "quay.io/biocontainers/truvari@sha256:022be7a2dfe5decf1aff8d007a914f56ed68d0bc224888e1764820922bc11935"
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

task filter_parent_variants {
  meta {
    description: "Filter variants present in parent from sample using truvari consistency results"
  }

  parameter_meta {
    sample_vcf: {
      name: "Sample VCF to filter"
    }
    sample_vcf_index: {
      name: "Sample VCF index"
    }
    truvari_consistency_tsv: {
      name: "Truvari consistency results"
    }
    sample_id: {
      name: "Sample ID"
    }
    parent_id: {
      name: "Parent sample ID"
    }
    sample_index: {
      name: "Sample index in consistency table"
    }
    parent_index: {
      name: "Parent index in consistency table"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    filtered_vcf: {
      name: "Parent-filtered VCF"
    }
    filtered_vcf_index: {
      name: "Parent-filtered VCF index"
    }
  }

  input {
    File sample_vcf
    File sample_vcf_index
    File truvari_consistency_tsv
    String sample_id
    String parent_id
    Int sample_index
    Int parent_index
    RuntimeAttributes runtime_attributes
  }

  Int threads = 2
  Int mem_gb = 8
  Int disk_size = ceil(size(sample_vcf, "GB") * 2 + size(truvari_consistency_tsv, "GB") + 20)
  String filtered_vcf_file = sub(basename(sample_vcf), "\\.vcf.gz$", "") + ".parent_filtered.vcf.gz"
  String filtered_vcf_index_file = "~{filtered_vcf_file}.tbi"
  String filtered_consistency_tsv = sub(basename(truvari_consistency_tsv), "\\.tsv$", "") + ".parent_filtered.tsv"

  command <<<
    set -euxo pipefail

    # Filter variants based on truvari consistency using variant IDs
    # truvari_consistency.tsv format: "CHROM    POS    ID    REF    ALT    FLAG    COUNT"
    if [[ -n "~{parent_index}" ]]; then
      # filter truvari_consistency.tsv to `(FLAG & sample_idx) && !(FLAG & parent_idx)
      sample_mask=$(( 1 << ~{sample_index} ))
      parent_mask=$(( 1 << ~{parent_index} ))
      flagcol=6
      awk -vflag="$flagcol" -vsample="$sample_mask" -vparent="$parent_mask" \
        'NR == 1 || and($flag,sample) && !and($flag,parent)' \
        ~{truvari_consistency_tsv} > ~{filtered_consistency_tsv}
    
      # assume position-sorted truvari_consistency.tsv
      bgzip ~{filtered_consistency_tsv}
      tabix -S1 -s1 -b2 -e2 ~{filtered_consistency_tsv}.gz

      cat <<- EOF
##INFO=<ID=FLAG,Number=1,Type=Integer,Description="Bitmask of samples from family with the variant">
##INFO=<ID=COUNT,Number=1,Type=Integer,Description="Count of samples from family with the variant">
EOF > truvari_consistency.hdr

      bcftools annotate \
        -a ~{filtered_consistency_tsv}.gz \
        -c 'CHROM,POS,~ID,REF,ALT,=FLAG,=COUNT' \
        -h truvari_consistency.hdr \
        -i 'FLAG!=0' \
        ~{sample_vcf} > ~{filtered_vcf_file}
    else
      # No parent specified, copy original VCF
      cp ~{sample_vcf} ~{filtered_vcf_file}
    fi
    
    # Index the filtered VCF
    tabix -p vcf ~{filtered_vcf_file}
  >>>

  output {
    File filtered_vcf = "~{filtered_vcf_file}"
    File filtered_vcf_index = "~{filtered_vcf_index_file}"
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