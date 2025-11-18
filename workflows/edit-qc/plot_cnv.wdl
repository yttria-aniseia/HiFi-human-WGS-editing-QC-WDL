version 1.1

import "../wdl-common/wdl/structs.wdl"

task cnvpytor_plot {
  meta {
    description: "Generate CNVpytor Manhattan plots from aligned BAM files"
  }

  parameter_meta {
    sample_id: "Sample identifier"
    aligned_bam: "Aligned BAM file from pbmm2"
    aligned_bam_index: "BAM index file"
    bin_size: "Bin size for CNVpytor analysis"
    runtime_attributes: "Runtime attributes for backend configuration"
  }

  input {
    String sample_id
    File aligned_bam
    File aligned_bam_index
    Int bin_size = 100000
    RuntimeAttributes runtime_attributes
  }

  Int threads = 4
  Int mem_gb = 16
  Int disk_size = ceil(size(aligned_bam, "GiB") * 3 + 50)

  command <<<
    set -euo pipefail

    # Define output names
    PYTOR_FILE="~{sample_id}.cnvpytor.pytor"
    CALLS_FILE="~{sample_id}.cnvpytor.calls.~{bin_size}.tsv"
    MANHATTAN_FILE="~{sample_id}.cnvpytor.manhattan.~{bin_size}.png"

    # Step 1: Import read depth from BAM
    echo "Step 1: Importing read depth from BAM..."
    cnvpytor -root "${PYTOR_FILE}" -rd "~{aligned_bam}" -chrom chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY

    # Step 2: Generate histograms
    echo "Step 2: Generating histograms..."
    cnvpytor -root "${PYTOR_FILE}" -his ~{bin_size}

    # Step 3: Partition
    echo "Step 3: Partitioning..."
    cnvpytor -root "${PYTOR_FILE}" -partition ~{bin_size}

    # Step 4: Call CNVs
    echo "Step 4: Calling CNVs..."
    cnvpytor -root "${PYTOR_FILE}" -call ~{bin_size} > "${CALLS_FILE}"

    # Step 5: Generate Manhattan plot
    echo "Step 5: Generating Manhattan plot..."
    cnvpytor -root "${PYTOR_FILE}" -plot manhattan ~{bin_size} -o "${MANHATTAN_FILE}"

    ls
    echo "CNVpytor analysis complete for sample ~{sample_id}"
  >>>

  output {
    File pytor_file = "~{sample_id}.cnvpytor.pytor"
    File calls_tsv = "~{sample_id}.cnvpytor.calls.~{bin_size}.tsv"
    File manhattan_plot = "~{sample_id}.cnvpytor.manhattan.~{bin_size}.global.0000.png"
  }

  runtime {
    docker: "quay.io/biocontainers/cnvpytor@sha256:42e480d51b4c640ebb46ede4ded05486c9974c0d0bde3a3cef9884ad56674383"
    cpu: threads
    memory: mem_gb + " GiB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
  }
}
