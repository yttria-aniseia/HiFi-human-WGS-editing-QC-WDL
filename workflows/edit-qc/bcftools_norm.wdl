version 1.1

import "../wdl-common/wdl/structs.wdl"

task bcftools_norm_split_multiallelic {
    input {
        File input_vcf
        File? input_vcf_index
        File ref_fasta
        File ref_fasta_index
        String out_prefix = "normalized"
        Int threads = 4
        Int mem_gb = 8
        RuntimeAttributes runtime_attributes
    }

    Float file_size = ceil(size(input_vcf, "GB") + size(ref_fasta, "GB") + size(ref_fasta_index, "GB"))
    String output_vcf = "~{out_prefix}.norm.vcf.gz"
    
    command <<<
        set -euxo pipefail
        
        bcftools --version

        bcftools sort \
            --max-mem "~{mem_gb}G" \
            -O z \
            ~{input_vcf} \
            | \
        bcftools norm \
            -m - \
            -f ~{ref_fasta} \
            --threads ~{threads} \
            -O z \
            | \
        bcftools sort \
            --max-mem "~{mem_gb}G" \
            -Wtbi \
            -O z \
            -o ~{output_vcf} \
    >>>
    
    output {
        File normalized_vcf = output_vcf
        File normalized_vcf_index = "~{output_vcf}.tbi"
    }
    
    runtime {
        docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87"
        cpu: threads
        memory: "~{mem_gb} GB"
        disk: file_size + " GB"
        maxRetries: 2
        preemptible: 1
    }
}