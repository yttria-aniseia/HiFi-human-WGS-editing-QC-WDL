version 1.1

import "../wdl-common/wdl/structs.wdl"

# Task to generate guidescan genome index
task guidescan_index {
  meta {
    description: "Generate guidescan index for the reference genome"
  }

  parameter_meta {
    ref_fasta: "Reference genome FASTA file"
    index_name: "Name for the guidescan index (e.g., 'hg38')"
  }

  input {
    File ref_fasta
    String index_name = "hg38"
    Int mem_gb = 16
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(3 * size(ref_fasta, "GB")) + 20

  command <<<
    set -euo pipefail

    # Generate guidescan index
    guidescan index --index ~{index_name} ~{ref_fasta}
  >>>

  output {
    Directory guidescan_index = index_name
  }

  runtime {
    docker: "quay.io/biocontainers/guidescan@sha256:9d93021243780b1faff47f1df4c1ae495177ff65ccc8273f0ec590caad5c82f0"
    cpu: 4
    memory: "~{mem_gb} GiB"
    disk: disk_size + " GB"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
  }
}

# Combined task to create kmer CSV, enumerate off-targets, and convert to BED
task guidescan_search {
  meta {
    description: "Search for CRISPR off-target sites: create kmer CSV, enumerate, and convert to BED"
  }

  parameter_meta {
    sample_id: "Sample identifier"
    crispr_edit_json: "JSON file describing expected CRISPR edit"
    guidescan_index: "Guidescan genome index directory"
    mismatches: "Maximum number of mismatches to allow"
    threads: "Number of threads for guidescan enumerate"
    region_extend: "Base pairs to extend region in both directions (plus gRNA length in strand direction)"
  }

  input {
    String sample_id
    File crispr_edit_json
    Directory guidescan_index
    Int mismatches = 4
    Int threads = 24
    Int mem_gb = 32
    Int region_extend = 20
    RuntimeAttributes runtime_attributes
  }

  command <<<
    set -euo pipefail

    # Step 1: Create kmer CSV from expected_edit JSON
    python3 <<'EOF1'
import json

with open("~{crispr_edit_json}", 'r') as f:
    edit_data = json.load(f)

# Write header
with open("~{sample_id}.kmers.csv", 'w') as outfile:
    outfile.write("id,sequence,pam,chromosome,position,sense\n")

    # Process each edit
    for edit in edit_data.get('edits', []):
        edit_id = edit.get('id', '~{sample_id}')
        target = edit.get('target', {})
        edit_info = edit.get('edit', {})

        grna = edit_info.get('guide', '')
        strand = edit_info.get('strand', '+')
        chrom = target.get('chr', '')

        # For position, use start of target region
        position = target.get('start', 0)

        # Extract PAM from the guide sequence (last 3 bases for SpCas9)
        if strand == '-':
            sense = '-'
            # For reverse strand, assume NGG PAM
            pam = 'NGG'
        else:
            sense = '+'
            # For forward strand, assume NGG PAM
            pam = 'NGG'

        outfile.write(f"{edit_id},{grna},{pam},{chrom},{position},{sense}\n")
EOF1

    # Step 2: Run guidescan enumerate
    guidescan enumerate \
      --alt-pam NGA,NAG \
      --mismatches ~{mismatches} \
      --kmers-file ~{sample_id}.kmers.csv \
      --output ~{sample_id}-offtarget.csv \
      --threads ~{threads} \
      --format csv \
      ~{guidescan_index}

    # Step 3: Convert CSV to BED format, excluding on-target matches
    python3 <<'EOF2'
import csv
import json

# Load expected edit locations to filter out on-target matches
with open("~{crispr_edit_json}", 'r') as f:
    edit_data = json.load(f)

# Create set of on-target regions to exclude
on_target_regions = []
for edit in edit_data.get('edits', []):
    target = edit.get('target', {})
    chrom = target.get('chr', '')
    start = target.get('start', 0)
    end = target.get('end', 0)
    on_target_regions.append((chrom, start, end))

# Read CSV and convert to BED
with open("~{sample_id}-offtarget.csv", 'r') as infile, \
     open("~{sample_id}-offtarget.bed", 'w') as outfile:

    reader = csv.DictReader(infile)

    for row in reader:
        chrom = row.get('match_chrm', '')
        pos = int(row.get('match_position', 0))
        strand = row.get('match_strand', '+')
        sequence = row.get('match_sequence', '')
        distance = row.get('match_distance', '0')

        # Skip on-target matches (distance 0 and overlapping expected edit)
        if distance == '0':
            is_on_target = False
            for target_chr, target_start, target_end in on_target_regions:
                if chrom == target_chr and target_start <= pos <= target_end:
                    is_on_target = True
                    break
            if is_on_target:
                continue

        # Calculate BED start and end based on strand
        # Extend by region_extend bp in both directions, plus gRNA length (20bp) in strand direction
        region_extend = ~{region_extend}
        grna_length = 20
        if strand == '-':
            # For - strand: extend region_extend + gRNA length upstream (toward lower coordinates)
            start = pos - region_extend - grna_length
            end = pos + region_extend
        else:
            # For + strand: extend region_extend + gRNA length downstream (toward higher coordinates)
            start = pos - region_extend
            end = pos + region_extend + grna_length

        # Ensure start is not negative
        start = max(0, start)

        # Write BED line (tab-delimited)
        outfile.write(f"{chrom}\t{start}\t{end}\t{sequence}\t{distance}\t{strand}\n")
EOF2
  >>>

  output {
    File kmer_csv = "~{sample_id}.kmers.csv"
    File offtarget_csv = "~{sample_id}-offtarget.csv"
    File offtarget_bed = "~{sample_id}-offtarget.bed"
  }

  runtime {
    docker: "quay.io/biocontainers/guidescan@sha256:9d93021243780b1faff47f1df4c1ae495177ff65ccc8273f0ec590caad5c82f0"
    cpu: threads
    memory: "~{mem_gb} GiB"
    disk: "50 GB"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
  }
}
