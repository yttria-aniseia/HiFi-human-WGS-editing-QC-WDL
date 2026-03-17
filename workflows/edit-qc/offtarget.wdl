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
    Int mem_gb = 32
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(3 * size(ref_fasta, "GB")) + 20

  command <<<
    set -euo pipefail

    # Generate guidescan index
    # This creates three files: <index_name>.forward, <index_name>.gs, <index_name>.reverse
    guidescan index --index ~{index_name} ~{ref_fasta}
  >>>

  output {
    File index_forward = "~{index_name}.forward"
    File index_gs = "~{index_name}.gs"
    File index_reverse = "~{index_name}.reverse"
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

# Task to create kmer CSV from CRISPR edit JSON
task create_kmer_csv {
  meta {
    description: "Create kmer CSV file from CRISPR edit JSON for guidescan"
  }

  parameter_meta {
    sample_id: "Sample identifier"
    crispr_edit_json: "JSON file describing expected CRISPR edit"
  }

  input {
    String sample_id
    File crispr_edit_json
    RuntimeAttributes runtime_attributes
  }

  command <<<
    set -euo pipefail

    python3 <<'EOF'
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
EOF
  >>>

  output {
    File kmer_csv = "~{sample_id}.kmers.csv"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87"
    cpu: 2
    memory: "4 GiB"
    disk: "10 GB"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
  }
}

# Task to enumerate off-target sites using guidescan
task guidescan_enumerate {
  meta {
    description: "Enumerate CRISPR off-target sites using guidescan"
  }

  parameter_meta {
    sample_id: "Sample identifier"
    kmer_csv: "CSV file containing gRNA sequences and positions"
    index_forward: "Guidescan index forward file"
    index_gs: "Guidescan index gs file"
    index_reverse: "Guidescan index reverse file"
    index_prefix: "Prefix for the guidescan index files"
    mismatches: "Maximum number of mismatches to allow"
    threads: "Number of threads for guidescan enumerate"
  }

  input {
    String sample_id
    File kmer_csv
    File index_forward
    File index_gs
    File index_reverse
    String index_prefix = "index"
    Int mismatches = 4
    Int threads = 24
    Int mem_gb = 32
    RuntimeAttributes runtime_attributes
  }

  command <<<
    set -euo pipefail

    # Symlink the index files with the correct prefix
    ln -s ~{index_forward} ~{index_prefix}.forward
    ln -s ~{index_gs} ~{index_prefix}.gs
    ln -s ~{index_reverse} ~{index_prefix}.reverse

    # Run guidescan enumerate
    guidescan enumerate \
      --alt-pam NGA,NAG \
      --mismatches ~{mismatches} \
      --kmers-file ~{kmer_csv} \
      --output ~{sample_id}-offtarget.csv \
      --threads ~{threads} \
      --format csv \
      ~{index_prefix}
  >>>

  output {
    File offtarget_csv = "~{sample_id}-offtarget.csv"
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

# Task to convert guidescan CSV output to BED format
task convert_offtarget_to_bed {
  meta {
    description: "Convert guidescan off-target CSV to BED format, excluding on-target matches"
  }

  parameter_meta {
    sample_id: "Sample identifier"
    offtarget_csv: "CSV file from guidescan enumerate"
    crispr_edit_json: "JSON file describing expected CRISPR edit"
  }

  input {
    String sample_id
    File offtarget_csv
    File crispr_edit_json
    RuntimeAttributes runtime_attributes
  }

  command <<<
    set -euo pipefail

    python3 <<'EOF'
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

# Read CSV, convert to BED, and find best match for expected cut site
best_match = None
with open("~{offtarget_csv}", 'r') as infile, \
     open("~{sample_id}-offtarget.bed", 'w') as outfile:

    reader = csv.DictReader(infile)

    for row in reader:
        chrom = row.get('match_chrm', '')
        pos = int(row.get('match_position', 0))
        strand = row.get('match_strand', '+')
        sequence = row.get('match_sequence', '')
        distance_str = row.get('match_distance', '0')
        try:
            distance = float(distance_str)
        except ValueError:
            distance = float('inf')

        # Track best match (minimum distance) for expected cut site
        if best_match is None or distance < best_match['distance']:
            best_match = {'chrom': chrom, 'pos': pos, 'strand': strand, 'distance': distance}

        # Skip on-target matches (distance 0 and overlapping expected edit)
        if distance == 0:
            is_on_target = False
            for target_chr, target_start, target_end in on_target_regions:
                if chrom == target_chr and target_start <= pos <= target_end:
                    is_on_target = True
                    break
            if is_on_target:
                continue

        # Calculate BED start and end (0-based): adjust match_position by -1
        grna_length = 20
        pos_0 = pos - 1
        if strand == '-':
            start = pos_0 - grna_length
            end = pos_0
        else:
            start = pos_0
            end = pos_0 + grna_length

        # Ensure start is not negative
        start = max(0, start)

        # Write BED line (tab-delimited)
        outfile.write(f"{chrom}\t{start}\t{end}\t{sequence}\t{distance_str}\t{strand}\n")

# Compute expected cut site from best match (minimum mismatches)
# For + strand: cut site = match_position + 17
# For - strand: cut site = match_position + 3
if best_match:
    bm_strand = best_match['strand']
    cut_offset = 17 if bm_strand == '+' else 3
    cut_pos = best_match['pos'] + cut_offset
    expected_cut_site = f"{best_match['chrom']}:{cut_pos}"
    expected_cut_strand = bm_strand
else:
    expected_cut_site = ""
    expected_cut_strand = ""

with open("~{sample_id}_expected_cut_site.txt", 'w') as f:
    f.write(expected_cut_site)
with open("~{sample_id}_expected_cut_strand.txt", 'w') as f:
    f.write(expected_cut_strand)
EOF
  >>>

  output {
    File offtarget_bed = "~{sample_id}-offtarget.bed"
    String expected_cut_site = read_string("~{sample_id}_expected_cut_site.txt")
    String expected_cut_strand = read_string("~{sample_id}_expected_cut_strand.txt")
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87"
    cpu: 2
    memory: "4 GiB"
    disk: "20 GB"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
  }
}
