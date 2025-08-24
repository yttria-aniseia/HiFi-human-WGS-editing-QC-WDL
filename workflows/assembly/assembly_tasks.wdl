version 1.1

import "../wdl-common/wdl/structs.wdl"
task samtools_bam_to_fasta {
    input {
        File input_bam
        Int threads
        RuntimeAttributes runtime_attributes
    }

    Int mem_gb    = ceil(threads * 4)
    Int disk_size = ceil(size(input_bam, "GB") * 3 + 70)
    command <<<
        set -euxo pipefail
        samtools index -@ ~{threads} ~{input_bam}
        name=$(basename "~{input_bam}" .bam)
        samtools fasta -@ ~{threads} ~{input_bam} > ${name}.fasta
    >>>

    output {
        File input_fasta = sub(basename(input_bam), "\\.bam$","") + ".fasta"
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


task hifiasm_assembly {
    input {
        File input_fasta
        Int threads
        RuntimeAttributes runtime_attributes
    }
    Int mem_gb    = ceil(threads * 6)
    Int disk_size = ceil(size(input_fasta, "GB") * 3 + 70)

    command <<<
        set -euxo pipefail
        name=$(basename "~{input_fasta}" .fasta)
        hifiasm -o ${name}.asm -t ~{threads} ~{input_fasta}
    >>>

    output {
        File input_1_asm = sub(basename(input_fasta), "\\.fasta$", "") + ".asm.bp.hap1.p_ctg.gfa"
        File input_2_asm = sub(basename(input_fasta), "\\.fasta$", "") + ".asm.bp.hap2.p_ctg.gfa"
    }
    runtime {
    docker: "quay.io/biocontainers/hifiasm@sha256:5dc4c88cabceb56445f44e785dae252e13fb8131e9bde54028bfb4102a3f424d"
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

task gfa_to_fa {
    input {
        File input_hap1_gfa
        File input_hap2_gfa
        RuntimeAttributes runtime_attributes
    }
    Int threads = 32
    Int mem_gb    = ceil(threads * 4)
    Int disk_size = ceil(size(input_hap1_gfa, "GB") * 3 + size(input_hap2_gfa, "GB") + 70)

    command <<<
        set -euxo pipefail
        name=$(basename "~{input_hap1_gfa}" .asm.bp.hap1.p_ctg.gfa)
        gfatools gfa2fa ~{input_hap1_gfa} > ${name}.asm.bp.hap1.p_ctg.fa
        gfatools gfa2fa ~{input_hap2_gfa} > ${name}.asm.bp.hap2.p_ctg.fa
        bgzip ${name}.asm.bp.hap1.p_ctg.fa
        bgzip ${name}.asm.bp.hap2.p_ctg.fa
    >>>

    output {
        File fasta_hap1 = sub(basename(input_hap1_gfa), "\\.gfa$", "") + ".fa.gz"
        File fasta_hap2 = sub(basename(input_hap2_gfa), "\\.gfa$", "") + ".fa.gz"
    }
    runtime {
    docker: "quay.io/pacbio/gfatools@sha256:da024e74381236932a432f4879bed853be4decb62e0b522e8a3206b5a5c1e527"
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

task quast {
    input {
        File input_fa_1       # haplotype 1 assembly
        File input_fa_2       # haplotype 2 assembly
        File input_fa         # full fasta before assembly
        File ref
        RuntimeAttributes runtime_attributes
    }
    Int threads = 32
    Int mem_gb    = ceil(threads * 4)
    Int disk_size = ceil(size(input_fa_1, "GB") * 3 + size(input_fa_2, "GB") + 70)
    String name=sub(basename(input_fa_1), "\\.asm.bp.hap1.p_ctg.fa.gz$", "") + "_quast"

    command <<<
        set -euxo pipefail
        name=$(basename "~{input_fa_1}" .asm.bp.hap1.p_ctg.fa.gz)
        mkdir ${name}_quast
        quast ~{input_fa_1} ~{input_fa_2} --pacbio ~{input_fa} -r ~{ref} --threads ~{threads} --conserved-genes-finding --no-sv --no-snps -o ${name}_quast
    >>>

    output {
       #File quast_output = sub(basename(input_fa_1), "\\.asm.bp.hap1.p_ctg.fa.gz$", "") + "_quast"
       # Array[File] quast_output = glob(name+"/**/*")
       #Directory quast_output = name + "_quast"
       File transposed_report = name + "/transposed_report.tsv"
       File icarus_html = name + "/icarus.html"
       File report_html = name + "/report.html"
       File report_pdf = name + "/report.pdf"
       File report_txt = name + "/report.txt"
    }
    runtime {
    docker: "quay.io/biocontainers/quast@sha256:a8d5771c90e4fb6daf32885c56ccece7523344c4e7cd02aa1601a75dff885c22"
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

task pav {
  input {
    File ref  # The HG002 reference FASTA
    File ref_index
    File ref_bgzf_index
    File input_fa_1       # haplotype 1 assembly
    File input_fa_2       # haplotype 2 assembly
    String sample_name
    String pav_sif_path = "/hpc/scratch/group.data.science/ram_temp/HiFi-human-WGS-editing-QC-WDL/miniwdl_cache/singularity_cache/pav_latest.sif/pav_latest.sif"  # Pre-downloaded container
    Boolean no_temp_cleanup = true
    RuntimeAttributes runtime_attributes
  }

  Int threads = 32
  Int mem_gb = ceil(threads * 4)
  Int disk_size = ceil(size(input_fa_1, "GB") * 2 + size(input_fa_2, "GB") * 2 + size(ref, "GB") + 20)
  command <<<
    set -euo pipefail

    # Create PAV output directory
    OUTPUT_DIR="pav_outputs_~{sample_name}"
    mkdir -p "${OUTPUT_DIR}"
    mkdir -p "${OUTPUT_DIR}/ref"
    mkdir -p "${OUTPUT_DIR}/assemblies"

    # Copy HG002 reference
    cp "~{ref}" "${OUTPUT_DIR}/ref/"
    cp "~{ref_index}" "${OUTPUT_DIR}/ref/"
    cp "~{ref_bgzf_index}" "${OUTPUT_DIR}/ref/"
    REF_BASENAME=$(basename "~{ref}")

    # Copy and prepare assembly files
    HAP1_BASENAME=$(basename "~{input_fa_1}")
    HAP2_BASENAME=$(basename "~{input_fa_2}")

    cp "~{input_fa_1}" "${OUTPUT_DIR}/assemblies/"
    cp "~{input_fa_2}" "${OUTPUT_DIR}/assemblies/"


    # Create config.json for HG002
    cat <<EOF > "${OUTPUT_DIR}/config.json"
    {
      "reference": "ref/${REF_BASENAME}"
    }
    EOF

    # Create assemblies.tsv
    # Clean sample name (replace dots with underscores)
    CLEAN_NAME=$(echo "~{sample_name}" | sed 's/\./_/g')

    # Extract base names for assemblies.tsv (removing .asm.bp.hap1.p_ctg pattern)
    INPUT_1=$(echo "${HAP1_BASENAME}" | sed 's/\.asm\.bp\.hap1\.p_ctg\.fa\.gz$//')
    INPUT_2=$(echo "${HAP2_BASENAME}" | sed 's/\.asm\.bp\.hap2\.p_ctg\.fa\.gz$//')

    cat <<EOF > "${OUTPUT_DIR}/assemblies.tsv"
    NAME	HAP1	HAP2
    ${CLEAN_NAME}	assemblies/${HAP1_BASENAME}	assemblies/${HAP2_BASENAME}
    EOF

    # Change to PAV directory and run
    cd "${OUTPUT_DIR}"

    # --nt parameter prevents the .fai read error common in pav 2.x.x, allegedly will be fixed in pav 3
    snakemake -s /opt/pav/Snakefile --nt ${CLEAN_NAME}.vcf.gz
  >>>

  output {
    File? pav_vcf = "pav_outputs_~{sub(sample_name, "\\.", "_")}/~{sub(sample_name, "\\.", "_")}.vcf.gz"
    File? pav_vcf_index = "pav_outputs_~{sub(sample_name, "\\.", "_")}/~{sub(sample_name, "\\.", "_")}.vcf.gz.tbi"
  }

  runtime {
    docker: "becklab/pav"
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


task liftover {
    input {
        File vcf
        File chain
        File ref
        RuntimeAttributes runtime_attributes
    }
    Int threads   = 2
    Int mem_gb    = ceil(threads * 2)
    Int disk_size = ceil(size(vcf, "GB") * 3 + size(ref, "GB") + 70)
    command <<<
        set -euxo pipefail

        name=$(basename "~{vcf}" .vcf.gz)
        OUTPUT="${name}_liftover_hg002_hg38.vcf"
        CrossMap vcf ~{chain} ~{vcf} ~{ref} $OUTPUT
    >>>

    output {
        File lifted_vcf = sub(basename(vcf), "\\.vcf.gz$", "") + "_liftover_hg002_hg38.vcf"
        File reject_vcf = sub(basename(vcf), "\\.vcf.gz$", "") + "_liftover_hg002_hg38.vcf.unmap"

    }
    runtime {
    docker: "quay.io/biocontainers/crossmap@sha256:cc8bbb5eaa8028ebf0f44a2499d1342c7e94e1989d4c3d7e0cc456e36b12053a"
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

task pav_sv_filter_size {
    input {
        File vcf
        RuntimeAttributes runtime_attributes
    }
    Int threads = 6
    Int mem_gb = 4
    Int disk_size = 500
    command <<<
        set -euxo pipefail

        name=$(basename "~{vcf}" .vcf)
        OUTPUT="${name}.large_sv_filtered.vcf"
        bcftools view -i 'INFO/SVTYPE !="." && (INFO/SVLEN >= 1000 || INFO/SVLEN <= -1000)' ~{vcf} | bcftools sort -Oz -o ${OUTPUT}.gz
        tabix -p vcf ${OUTPUT}.gz
    >>>

    output {
        File large_sv_filtered_vcf = sub(basename(vcf), "\\.vcf$", "") + ".large_sv_filtered.vcf.gz"
        File large_sv_filtered_vcf_index = sub(basename(vcf), "\\.vcf$", "") + ".large_sv_filtered.vcf.gz.tbi"
    }
    runtime {
    docker: "quay.io/biocontainers/bcftools:1.17--h3cc50cf_1"
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
