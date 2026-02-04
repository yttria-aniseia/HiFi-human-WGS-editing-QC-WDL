version 1.1

import "../wdl-common/wdl/structs.wdl"

task create_parts_query_fasta {
  meta {
    description: "Create FASTA file with edit parts (left_ha, payload, right_ha) from CRISPR edit JSON"
  }

  parameter_meta {
    crispr_edit_json: "JSON file describing CRISPR edit with left_ha, payload, right_ha sequences"
    include_payload_components: "Include payload_components as separate FASTA entries for detailed annotation"
  }

  input {
    File crispr_edit_json
    Boolean include_payload_components = false
    RuntimeAttributes runtime_attributes
  }

  Int threads = 1
  String out_prefix = basename(crispr_edit_json, ".json")

  command <<<
    set -euxo pipefail

    python3 <<'PYEOF'
import json

with open("~{crispr_edit_json}", 'r') as f:
    crispr_data = json.load(f)

include_components = ~{true="True" false="False" include_payload_components}

with open("~{out_prefix}_parts_query.fasta", 'w') as fasta_out:
    for edit in crispr_data['edits']:
        edit_id  = edit['id'].replace(" ", "_")
        left_ha  = edit['edit']['left_ha']
        payload  = edit['edit']['payload']
        right_ha = edit['edit']['right_ha']

        fasta_out.write(f">{edit_id}_left_HA\n{left_ha}\n")
        fasta_out.write(f">{edit_id}_payload\n{payload}\n")
        fasta_out.write(f">{edit_id}_right_HA\n{right_ha}\n")
        fasta_out.write(f">{edit_id}_HDR_template\n{left_ha}{payload}{right_ha}\n")
        fasta_out.write(f">{edit_id}_WT_template\n{left_ha}{right_ha}\n")

        if include_components and 'payload_components' in edit['edit']:
            for component in edit['edit']['payload_components']:
                comp_name = component['name'].replace(" ", "_")
                fasta_out.write(f">{comp_name}\n{component['seq']}\n")
PYEOF
>>>

  output {
    File parts_query_fasta = "~{out_prefix}_parts_query.fasta"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87"
    cpu: threads
    memory: "4 GB" 
    disk: "10 GB"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
  }
}

task minimap2_align_reads {
  meta {
    description: "Use minimap2 to identify ALL reads containing edit parts"
  }

  parameter_meta {
    fasta_reads: "FASTA converted reads"
    parts_query_fasta: "FASTA file with edit parts sequences"
    min_identity: "Minimum alignment identity (0-1)"
    threads: "Number of threads"
  }

  input {
    File fasta_reads
    File parts_query_fasta
    Float min_identity = 0.80
    Int threads = 8
    Int mem_gb = 48
    String sample_id
    RuntimeAttributes runtime_attributes
  }

  Float file_size = ceil(size(fasta_reads, "GB") * 2)
  String out_prefix = "~{sample_id}_edit_filter"

  command <<<
    set -euxo pipefail

    # Initial alignment to identify reads with any edit parts
    minimap2 -a -x map-hifi -p ~{min_identity} -t ~{threads} \
      ~{fasta_reads} ~{parts_query_fasta} | grep -v "^@" | cut -f3 > ~{out_prefix}_hit_names.txt

    # Count alignments
    echo "Total edit part hits: $(wc -l < ~{out_prefix}_hit_names.txt)"
  >>>

  output {
    File matched_read_names = "~{out_prefix}_hit_names.txt"
  }

  runtime {
    docker: "quay.io/biocontainers/minimap2@sha256:fdc9ef8bfbd31bab59a61b2e90dd226647deed971556175a4cd004f0bcdc7608"
    cpu: threads
    memory: "~{mem_gb} GB"
    disk: file_size + " GB" 
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
  }
}

task extract_reads {
  meta {
    description: "filter reads based on alignment results"
  }

  parameter_meta {
    fasta_reads: "FASTA converted reads (gzip compressed)"
    matched_read_names: "Read names from minimap2 alignment hits"
  }

  input {
    File fasta_reads
    File matched_read_names
    String sample_id
    RuntimeAttributes runtime_attributes
  }

  String out_prefix = "~{sample_id}_edit_filter"

  command <<<
    set -euxo pipefail

    # Extract matched reads to new FASTA using awk
    sort -u ~{matched_read_names} > matched_reads_uniq.txt
    awk 'NR==FNR {h[$1]; next} /^>/ {p = substr($0, 2) in h} p' matched_reads_uniq.txt <(gzip -dc ~{fasta_reads})  > ~{out_prefix}_filtered_reads.fasta

    # Count reads
    MATCHED_READS=$(grep -c '^>' ~{out_prefix}_filtered_reads.fasta || echo 0)
    echo "Reads matching edit parts: $MATCHED_READS"

    # Convert FASTA to gzip FASTQ with synthetic Q40 quality scores (for tools requiring FASTQ)
    python3 <<'PYEOF'
import gzip
name, buf = None, []
with open("~{out_prefix}_filtered_reads.fasta") as fh, \
     gzip.open("~{out_prefix}_filtered_reads.fastq.gz", 'wt') as out:
    for line in fh:
        line = line.rstrip()
        if line.startswith('>'):
            if name is not None:
                seq = ''.join(buf)
                out.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")
            name, buf = line[1:], []
        else:
            buf.append(line)
    if name is not None:
        seq = ''.join(buf)
        out.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")
PYEOF
  >>>

  output {
    File filtered_reads_fasta = "~{out_prefix}_filtered_reads.fasta"
    File filtered_reads_fastq = "~{out_prefix}_filtered_reads.fastq.gz"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87"
    cpu: 2
    memory: "4 GB"
    disk: "20 GB"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
  }
}


task extract_region_reads_fasta {
  meta {
    description: "Extract reads from a haplotagged BAM overlapping the edit target region; output FASTA"
  }

  parameter_meta {
    haplotagged_bam:       "Haplotagged BAM file from phasing"
    haplotagged_bam_index: "Index for haplotagged BAM"
    crispr_edit_json:      "JSON file describing CRISPR edit with target region"
    padding:               "Bases to add on each side of the target region"
  }

  input {
    File haplotagged_bam
    File haplotagged_bam_index
    File crispr_edit_json
    String sample_id
    Int padding = 10000
    RuntimeAttributes runtime_attributes
  }

  String out_prefix = "~{sample_id}_region_reads"

  command <<<
    set -euxo pipefail

    python3 <<'PYEOF'
import json
with open("~{crispr_edit_json}") as fh:
    data = json.load(fh)
edit = data['edits'][0]
chrom = edit['target']['chr']
start = max(1, edit['target']['start'] - ~{padding})
end   = edit['target']['end'] + ~{padding}
with open('region.txt', 'w') as fh:
    fh.write(f"{chrom}:{start}-{end}\n")
PYEOF

    REGION=$(cat region.txt)
    echo "Extracting reads from ${REGION}"

    samtools view -h ~{haplotagged_bam} "${REGION}" | \
      samtools fasta -F 2304 - \
      > ~{out_prefix}.fasta

    READ_COUNT=$(grep -c '^>' ~{out_prefix}.fasta || echo 0)
    echo "Extracted ${READ_COUNT} reads"

    # Synthetic Q40 FASTQ for tools that require FASTQ
    python3 <<'PYEOF'
import gzip
name, buf = None, []
with open("~{out_prefix}.fasta") as fh, \
     gzip.open("~{out_prefix}.fastq.gz", 'wt') as out:
    for line in fh:
        line = line.rstrip()
        if line.startswith('>'):
            if name is not None:
                seq = ''.join(buf)
                out.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")
            name, buf = line[1:], []
        else:
            buf.append(line)
    if name is not None:
        seq = ''.join(buf)
        out.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")
PYEOF
  >>>

  output {
    File region_reads_fasta = "~{out_prefix}.fasta"
    File region_reads_fastq = "~{out_prefix}.fastq.gz"
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


task blastn_remap_parts {
  meta {
    description: "Remap edit parts to filtered reads using blastn"
  }

  parameter_meta {
    filtered_reads_fasta: "FASTA with reads containing edit parts (subject)"
    parts_query_fasta: "FASTA file with edit parts sequences (query)"
  }

  input {
    File filtered_reads_fasta
    File parts_query_fasta
    String sample_id
    RuntimeAttributes runtime_attributes
  }

  String out_prefix = "~{sample_id}_edit_parts"

  # Tabular columns: query(part) coords, subject(read) coords, alignment stats, strand
  # Coordinates are 1-based; strand is "plus"/"minus".
  # word_size=32 prevents spurious 31 bp hits from XTEN-family repeat sequences
  # that share k-mers with the human genome (validated on czML1116 EEA1 dataset).
  # r=1 p=-2 go=2 ge=1: best CIGAR accuracy in ground-truth evaluation.
  String blastn_outfmt = "6 qseqid qlen qstart qend sseqid slen sstart send nident mismatch gapopen gaps length sstrand bitscore"

  command <<<
    set -euxo pipefail

    blastn \
      -query ~{parts_query_fasta} \
      -subject ~{filtered_reads_fasta} \
      -outfmt "~{blastn_outfmt}" \
      -task blastn \
      -word_size 32 \
      -reward 1 \
      -penalty -2 \
      -gapopen 2 \
      -gapextend 1 \
      -evalue 1 \
      -perc_identity 80 \
      -max_target_seqs 10000 \
      -parse_deflines \
      > ~{out_prefix}_blastn_crispr_parts.tsv

    echo "Total alignments: $(wc -l < ~{out_prefix}_blastn_crispr_parts.tsv)"
  >>>

  output {
    File parts_alignment_tsv = "~{out_prefix}_blastn_crispr_parts.tsv"
  }

  runtime {
    docker: "quay.io/biocontainers/blast@sha256:9dfc69f990c0aeb936276ee591ed32919a79f46dfa34060055e05a050a17959c"
    cpu: 2
    memory: "8 GiB"
    disk: "20 GB"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
  }
}



task clip_reads_to_parts {
  meta {
    description: "Clip filtered reads to the union of all blastn part-hit spans plus padding"
  }

  parameter_meta {
    parts_alignment_tsv: "Pass-1 blastn tabular output (from blastn_remap_parts)"
    filtered_reads_fasta: "Uncompressed FASTA with reads overlapping edit parts"
    clip_padding: "Bases to retain on each side of the union hit span (default 100)"
    min_part_hit_bp: "Minimum query-span (bp) for a blastn hit to be counted; filters repetitive short matches"
  }

  input {
    File parts_alignment_tsv
    File filtered_reads_fasta
    String sample_id
    Int clip_padding = 100
    Int min_part_hit_bp = 100
    RuntimeAttributes runtime_attributes
  }

  String out_prefix = "~{sample_id}_clipped"

  command <<<
    set -euxo pipefail

    python3 <<'PYEOF'
import sys
from collections import defaultdict

# Parse blastn TSV: collect hit t-spans per read (all parts, all hits)
# Columns: qseqid qlen qstart qend sseqid slen sstart send nident mismatch gapopen gaps length sstrand bitscore
hit_spans  = defaultdict(list)   # read_name -> [(t_start, t_end), ...]
read_lens  = {}

MIN_HIT = ~{min_part_hit_bp}

with open("~{parts_alignment_tsv}") as fh:
    for ln in fh:
        if not ln.strip() or ln.startswith('#'):
            continue
        f = ln.rstrip('\n').split('\t')
        if len(f) < 15:
            continue
        if abs(int(f[3]) - int(f[2])) < MIN_HIT:
            continue
        read_name = f[4]
        sstart, send = int(f[6]), int(f[7])
        t_start = min(sstart, send) - 1   # 0-based half-open
        t_end   = max(sstart, send)
        hit_spans[read_name].append((t_start, t_end))
        read_lens[read_name] = int(f[5])

padding = ~{clip_padding}

# Compute clip window per read
clip_regions = {}
for rname, spans in hit_spans.items():
    rlen = read_lens[rname]
    cs = max(0,    min(s[0] for s in spans) - padding)
    ce = min(rlen, max(s[1] for s in spans) + padding)
    clip_regions[rname] = (cs, ce)

# Stream FASTA and clip
n_written = 0
name = None
buf  = []

with open("~{filtered_reads_fasta}") as fh, \
     open("~{out_prefix}.fasta", 'w') as out:
    def flush(name, buf):
        global n_written
        if name not in clip_regions:
            return
        seq = ''.join(buf)
        cs, ce = clip_regions[name]
        out.write(f">{name} clip={cs}-{ce} orig_len={len(seq)}\n{seq[cs:ce]}\n")
        n_written += 1
    for ln in fh:
        ln = ln.rstrip('\n')
        if ln.startswith('>'):
            if name is not None:
                flush(name, buf)
            name = ln[1:].split()[0]
            buf  = []
        else:
            buf.append(ln.upper())
    if name is not None:
        flush(name, buf)

print(f"Clipped {n_written} reads (padding={padding} bp)")
PYEOF
  >>>

  output {
    File clipped_reads_fasta = "~{out_prefix}.fasta"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87"
    cpu: 1
    memory: "8 GiB"
    disk: "20 GB"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
  }
}


task structural_grouping {
  meta {
    description: "Derive structural signatures from blastn part hits and cluster reads into groups"
  }

  parameter_meta {
    parts_alignment_tsv:  "Pass-1 blastn tabular output"
    gap_merge_threshold:  "Max gap difference (bp) between two reads to merge into the same structural group"
    min_part_hit_bp:      "Minimum query-span (bp) for a blastn hit to be counted; filters repetitive short matches"
  }

  input {
    File parts_alignment_tsv
    String sample_id
    Int gap_merge_threshold = 5
    Int min_part_hit_bp     = 100
    RuntimeAttributes runtime_attributes
  }

  String out_prefix = "~{sample_id}_structural"

  command <<<
    set -euxo pipefail

    python3 <<'PYEOF'
from collections import defaultdict, Counter

# ── Parse blastn TSV ──────────────────────────────────────────────────────────
# Collect hits for individual parts only (left_HA, payload, right_HA).
# Skip HDR_template and WT_template — they are composite sequences and their
# hits are redundant with (and noisier than) the per-part hits.

hits_by_read = defaultdict(list)   # read_name -> [(part_base, strand, t_start, t_end)]

MIN_HIT = ~{min_part_hit_bp}

with open("~{parts_alignment_tsv}") as fh:
    for ln in fh:
        if not ln.strip() or ln.startswith('#'):
            continue
        f = ln.rstrip('\n').split('\t')
        if len(f) < 15:
            continue
        qname     = f[0]
        read_name = f[4]
        if abs(int(f[3]) - int(f[2])) < MIN_HIT:
            continue
        sstart, send = int(f[6]), int(f[7])
        t_start = min(sstart, send) - 1
        t_end   = max(sstart, send)
        strand  = '+' if f[13] == 'plus' else '-'

        if   qname.endswith('_left_HA'):    part_base = 'left_HA'
        elif qname.endswith('_right_HA'):   part_base = 'right_HA'
        elif qname.endswith('_payload'):    part_base = 'payload'
        else:                               continue   # composite or unknown

        hits_by_read[read_name].append((part_base, strand, t_start, t_end))

# ── Build structural signature for each read ─────────────────────────────────
# Signature: list of (part_base, strand) in t_start order; gaps between them.

def build_sig(hits):
    hits_s = sorted(hits, key=lambda h: h[2])   # sort by t_start
    parts_strands = [(h[0], h[1]) for h in hits_s]
    lengths = [h[3] - h[2] for h in hits_s]
    gaps = [hits_s[i+1][2] - hits_s[i][3] for i in range(len(hits_s)-1)]
    return parts_strands, gaps, lengths

def sig_str(ps, gs, ls):
    tokens = []
    for i, (p, s) in enumerate(ps):
        tokens.append(f"{p}:{s}:{ls[i]}")
        if i < len(gs):
            tokens.append(f"gap:{gs[i]}")
    return '|'.join(tokens) if tokens else 'empty'

# ── Cluster reads into structural groups ──────────────────────────────────────
# Prototype-based greedy clustering.
# Two reads join the same group iff:
#   1. Identical ordered (part, strand) sequence
#   2. |gap_i(A) - gap_i(B)| <= gap_merge_threshold for every gap position i

GAP_THR = ~{gap_merge_threshold}

prototypes  = defaultdict(list)   # pattern_key -> [(prototype_gaps, group_id)]
group_id_ctr = 0
read_groups  = {}                 # read_name -> (gid, sig_string)

for read_name, hits in sorted(hits_by_read.items()):
    ps, gs, ls = build_sig(hits)
    key = tuple(ps)

    matched_gid = None
    for proto_gs, gid in prototypes[key]:
        if len(proto_gs) == len(gs) and all(abs(pg - g) <= GAP_THR
                                            for pg, g in zip(proto_gs, gs)):
            matched_gid = gid
            break

    if matched_gid is None:
        group_id_ctr += 1
        matched_gid = group_id_ctr
        prototypes[key].append((gs, matched_gid))

    read_groups[read_name] = (matched_gid, sig_str(ps, gs, ls))

# ── Write output ──────────────────────────────────────────────────────────────
group_sizes = Counter(gid for gid, _ in read_groups.values())

with open("~{out_prefix}_groups.tsv", 'w') as out:
    out.write("read_name\tgroup_id\tsignature\tgroup_size\n")
    for rname, (gid, sig) in sorted(read_groups.items()):
        out.write(f"{rname}\t{gid}\t{sig}\t{group_sizes[gid]}\n")

with open("~{out_prefix}_hits.tsv", 'w') as out:
    out.write("read_name\tpart\tstrand\tt_start\tt_end\tmatch_len\n")
    for rname, hits in sorted(hits_by_read.items()):
        for part, strand, t_start, t_end in sorted(hits, key=lambda h: h[2]):
            out.write(f"{rname}\t{part}\t{strand}\t{t_start}\t{t_end}\t{t_end - t_start}\n")

print(f"Reads grouped : {len(read_groups)}")
print(f"Groups formed : {group_id_ctr}")
PYEOF
  >>>

  output {
    File groups_tsv = "~{out_prefix}_groups.tsv"
    File hits_tsv   = "~{out_prefix}_hits.tsv"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87"
    cpu: 1
    memory: "4 GiB"
    disk: "10 GB"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
  }
}


task prepare_group_fastas {
  meta {
    description: "Partition clipped reads into per-group input FASTAs ready for MAFFT; emit singletons separately"
  }

  input {
    File clipped_reads_fasta
    File groups_tsv
    String sample_id
    RuntimeAttributes runtime_attributes
  }

  String out_prefix = "~{sample_id}_msa"

  command <<<
    set -euxo pipefail

    python3 <<'PYEOF'
import os, tarfile
from collections import defaultdict

groups = defaultdict(list)

with open("~{groups_tsv}") as fh:
    fh.readline()
    for ln in fh:
        f = ln.rstrip('\n').split('\t')
        if len(f) < 4:
            continue
        rname, gid = f[0], int(f[1])
        groups[gid].append(rname)

reads = {}
name = None
buf  = []
with open("~{clipped_reads_fasta}") as fh:
    for ln in fh:
        ln = ln.rstrip('\n')
        if ln.startswith('>'):
            if name is not None:
                reads[name] = ''.join(buf)
            name = ln[1:].split()[0]
            buf  = []
        else:
            buf.append(ln.upper())
    if name is not None:
        reads[name] = ''.join(buf)

os.makedirs("msa_inputs", exist_ok=True)
singletons_out = open("~{out_prefix}_singletons.fasta", 'w')

for gid in sorted(groups):
    members = [r for r in groups[gid] if r in reads]
    if not members:
        continue
    if len(members) == 1:
        seq = reads[members[0]]
        singletons_out.write(f">group_{gid}_singleton:{members[0]}\n{seq}\n")
        continue
    with open(f"msa_inputs/group_{gid}_input.fasta", 'w') as gf:
        for r in members:
            gf.write(f">{r}\n{reads[r]}\n")

singletons_out.close()

with tarfile.open("~{out_prefix}_inputs.tar.gz", "w:gz") as tar:
    tar.add("msa_inputs", arcname="msa_inputs")

print(f"Groups written: {sum(1 for gid in groups if len([r for r in groups[gid] if r in reads]) > 1)} multi-read, "
      f"{sum(1 for gid in groups if len([r for r in groups[gid] if r in reads]) == 1)} singleton")
PYEOF
  >>>

  output {
    File group_inputs_tarball = "~{out_prefix}_inputs.tar.gz"
    File singletons_fasta     = "~{out_prefix}_singletons.fasta"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87"
    cpu: 1
    memory: "8 GiB"
    disk: "20 GB"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
  }
}


task run_group_msa {
  meta {
    description: "Run MAFFT --auto on each per-group input FASTA; pure bash, no Python"
  }

  input {
    File group_inputs_tarball
    String sample_id
    RuntimeAttributes runtime_attributes
  }

  String out_prefix = "~{sample_id}_msa"

  command <<<
    set -euxo pipefail

    mkdir msa_inputs msa_outputs
    tar -xzf ~{group_inputs_tarball} -C msa_inputs --strip-components=1

    shopt -s nullglob
    for input_fasta in msa_inputs/group_*_input.fasta; do
      gid=$(basename "$input_fasta" _input.fasta | sed 's/^group_//')
      mafft --auto --quiet "$input_fasta" > "msa_outputs/group_${gid}.fasta"
    done

    tar -czf ~{out_prefix}_msa.tar.gz -C msa_outputs .
    echo "MAFFT done: $(ls msa_outputs/*.fasta 2>/dev/null | wc -l) groups aligned"
  >>>

  output {
    File msa_tarball = "~{out_prefix}_msa.tar.gz"
  }

  runtime {
    docker: "quay.io/biocontainers/mafft@sha256:b8ccc402a9156868fefda5caa8693f847b47dcc69aae848557f7115850e60790"
    cpu: 4
    memory: "16 GiB"
    disk: "50 GB"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
  }
}


task build_group_consensus {
  meta {
    description: "Parse per-group MSA FASTAs; build error-corrected consensus; emit group metadata"
  }

  parameter_meta {
    mask_min_support:  "Mask consensus positions where plurality support <= this AND coverage >= mask_min_coverage"
    mask_min_coverage: "Coverage threshold for masking (see mask_min_support)"
  }

  input {
    File msa_tarball
    File singletons_fasta
    File groups_tsv
    String sample_id
    Int mask_min_support  = 2
    Int mask_min_coverage = 5
    RuntimeAttributes runtime_attributes
  }

  String out_prefix = "~{sample_id}_msa"

  command <<<
    set -euxo pipefail

    python3 <<'PYEOF'
import os, re, tarfile
from collections import defaultdict

# ── Load group metadata ───────────────────────────────────────────────────────
groups         = defaultdict(list)
group_strand   = {}   # gid -> '+' or '-'
group_category = {}   # gid -> 'HDR'|'WT'|'partial'|'misintegration'|'other'

_FWD_HDR_PARTS = [('left_HA', '+'), ('payload', '+'), ('right_HA', '+')]
_REV_HDR_PARTS = [('right_HA', '-'), ('payload', '-'), ('left_HA', '-')]

def _classify_sig(sig):
    """Classify a group signature string into HDR / WT / partial / misintegration.
    Signature format: '{part}:{strand}:{len}|gap:{bp}|...'. Structural only — does
    not check gap sizes or flanking (partial consensus may represent a real
    haplotype for inclusion decisions even if reads are truncated)."""
    pairs = []
    for tok in sig.split('|'):
        if tok.startswith('gap:') or tok == 'empty':
            continue
        parts = tok.split(':')
        if len(parts) >= 2:
            pairs.append((parts[0], parts[1]))
    if not pairs:
        return 'other'
    strands = {s for _, s in pairs}
    if len(strands) != 1:
        return 'misintegration'
    parts_only = [p for p, _ in pairs]
    strand = next(iter(strands))
    if strand == '+':
        if parts_only == ['left_HA', 'payload', 'right_HA']:
            return 'HDR'
        if parts_only == ['left_HA', 'right_HA']:
            return 'WT'
        if parts_only in (['left_HA', 'payload'], ['payload', 'right_HA'],
                          ['left_HA'], ['payload'], ['right_HA']):
            return 'partial'
    else:
        if parts_only == ['right_HA', 'payload', 'left_HA']:
            return 'HDR'
        if parts_only == ['right_HA', 'left_HA']:
            return 'WT'
        if parts_only in (['payload', 'left_HA'], ['right_HA', 'payload'],
                          ['left_HA'], ['payload'], ['right_HA']):
            return 'partial'
    return 'misintegration'

with open("~{groups_tsv}") as fh:
    fh.readline()
    for ln in fh:
        f = ln.rstrip('\n').split('\t')
        if len(f) < 4:
            continue
        rname, gid, sig = f[0], int(f[1]), f[2]
        groups[gid].append(rname)
        if gid not in group_strand:
            m = re.search(r':([+-]):', sig)
            if m:
                group_strand[gid] = m.group(1)
        if gid not in group_category:
            group_category[gid] = _classify_sig(sig)

def header_tags(gid):
    tags = []
    if gid in group_strand:
        tags.append(f"strand={group_strand[gid]}")
    if gid in group_category:
        tags.append(f"category={group_category[gid]}")
    return (' ' + ' '.join(tags)) if tags else ''

MASK_SUPP = ~{mask_min_support}
MASK_COV  = ~{mask_min_coverage}

os.makedirs("msa_groups", exist_ok=True)

# ── Helpers ───────────────────────────────────────────────────────────────────
def plurality_consensus(seqs_dict):
    """Plurality-vote consensus with N-masking over {name: msa_seq}.
    Gaps are treated as votes for deletion: if gaps are the plurality at a
    column, that position is omitted from the consensus (correct for a
    haplotype where most reads carry a deletion or the column is an insertion
    in a minority of reads).
    """
    width = len(next(iter(seqs_dict.values())))
    bases = []
    for col in range(width):
        col_bases = [s[col] for s in seqs_dict.values()]
        gap_count = col_bases.count('-')
        non_gap   = [b for b in col_bases if b != '-']
        if not non_gap or gap_count >= len(non_gap):
            continue  # all-gap or gap-plurality column → omit from consensus
        freq = {}
        for b in non_gap:
            freq[b] = freq.get(b, 0) + 1
        best_base, best_count = max(freq.items(), key=lambda x: x[1])
        bases.append('N' if best_count <= MASK_SUPP and len(non_gap) >= MASK_COV
                     else best_base)
    return ''.join(bases)

def hamming(v1, v2):
    # Gap vs base counts as a mismatch; only N positions are skipped
    return sum(1 for a, b in zip(v1, v2)
               if 'N' not in (a, b) and a != b)

def retained_columns(seqs_dict):
    """Return MSA column indices that plurality_consensus keeps (non-gap-plurality).
    Used to restrict haplotype clustering to error-corrected positions only."""
    width = len(next(iter(seqs_dict.values())))
    cols = []
    for col in range(width):
        col_bases = [s[col] for s in seqs_dict.values()]
        gap_count = col_bases.count('-')
        non_gap   = [b for b in col_bases if b != '-']
        if non_gap and gap_count < len(non_gap):
            cols.append(col)
    return cols

def cluster_haplotypes(msa_seqs, cols=None):
    """
    Detect variable sites (>= 2 distinct non-N alleles — including '-' for
    indels — each in >= 2 reads) and cluster reads by identical allele vector
    (Hamming = 0, skipping N positions only).  Returns {read_name: hap_label}
    where label is 'A', 'B', etc.  Returns all 'A' if no variable sites exist
    or any resulting cluster has fewer than 2 reads.

    cols: restrict variable-site search to these MSA column indices (pass-1
    error-corrected columns).  Defaults to all columns.
    """
    rnames = list(msa_seqs.keys())
    width  = len(next(iter(msa_seqs.values())))
    if cols is None:
        cols = range(width)

    var_cols = []
    for col in cols:
        freq = {}
        for r in rnames:
            b = msa_seqs[r][col]
            if b != 'N':                   # gaps count as a real allele
                freq[b] = freq.get(b, 0) + 1
        if sum(1 for cnt in freq.values() if cnt >= 2) >= 2:
            var_cols.append(col)

    if not var_cols:
        return {r: 'A' for r in rnames}

    vecs = {r: tuple(msa_seqs[r][col] for col in var_cols) for r in rnames}
    clusters = []
    for r in rnames:
        for cluster in clusters:
            if hamming(vecs[r], vecs[next(iter(cluster))]) == 0:
                cluster.add(r)
                break
        else:
            clusters.append({r})

    # Require every cluster to have >= 2 reads; otherwise the split is unreliable
    if any(len(c) < 2 for c in clusters):
        return {r: 'A' for r in rnames}

    labels = [chr(65 + i) for i in range(len(clusters))]
    return {r: labels[i] for i, cluster in enumerate(clusters) for r in cluster}

# ── Sub-haplotype assignments (written at end) ────────────────────────────────
subhap_assignments = {}   # read_name -> "group_{gid}_hap_{label}"

# ── Singletons pass-through ───────────────────────────────────────────────────
name = None
buf  = []

with open("~{out_prefix}_consensus.fasta", 'w') as consensus_out, \
     open("~{out_prefix}_group_meta.tsv",  'w') as group_meta_out:
    group_meta_out.write("group_id\tn_reads\thas_consensus\n")

    def flush_singleton(header, buf):
        tag = header.split()[0]   # group_{gid}_singleton:{read_name}
        try:
            gid = int(tag.split('_')[1])
        except ValueError:
            return
        seq       = ''.join(buf)
        read_name = tag.split(':')[1] if ':' in tag else tag
        consensus_out.write(f">group_{gid}_singleton n_reads=1{header_tags(gid)}\n{seq}\n")
        with open(f"msa_groups/group_{gid}.fasta", 'w') as mf:
            mf.write(f">{read_name}\n{seq}\n")
        group_meta_out.write(f"{gid}\t1\tFALSE\n")
        subhap_assignments[read_name] = f"group_{gid}_hap_A"

    with open("~{singletons_fasta}") as fh:
        for ln in fh:
            ln = ln.rstrip('\n')
            if ln.startswith('>'):
                if name is not None:
                    flush_singleton(name, buf)
                name = ln[1:]
                buf  = []
            else:
                buf.append(ln.upper())
        if name is not None:
            flush_singleton(name, buf)

    # ── Multi-read groups: parse MSA, sub-haplotype, build consensus ─────────
    with tarfile.open("~{msa_tarball}") as tar:
        for member in tar.getmembers():
            bn = os.path.basename(member.name)
            if not bn.endswith('.fasta'):
                continue
            try:
                gid = int(bn.replace('group_', '').replace('.fasta', ''))
            except ValueError:
                continue

            fobj = tar.extractfile(member)
            if not fobj:
                continue

            msa_seqs = {}
            mname, mbuf = None, []
            for ln in fobj.read().decode('utf-8').splitlines():
                if ln.startswith('>'):
                    if mname is not None:
                        msa_seqs[mname] = ''.join(mbuf)
                    mname = ln[1:]
                    mbuf  = []
                else:
                    mbuf.append(ln.upper())
            if mname is not None:
                msa_seqs[mname] = ''.join(mbuf)

            if not msa_seqs:
                group_meta_out.write(f"{gid}\t{len(groups[gid])}\tFALSE\n")
                continue

            with open(f"msa_groups/group_{gid}.fasta", 'w') as mf:
                for rname, seq in msa_seqs.items():
                    mf.write(f">{rname}\n{seq}\n")

            # Pass 1: identify error-corrected columns (indel noise filtered by
            # plurality thresholds), then cluster haplotypes only on those columns.
            kept_cols  = retained_columns(msa_seqs)
            hap_map    = cluster_haplotypes(msa_seqs, kept_cols)  # read -> 'A', 'B', ...
            hap_labels = sorted(set(hap_map.values()))
            for r, lbl in hap_map.items():
                subhap_assignments[r] = f"group_{gid}_hap_{lbl}"

            if len(hap_labels) > 1:
                # Emit one consensus per haplotype
                for lbl in hap_labels:
                    hap_seqs = {r: s for r, s in msa_seqs.items()
                                if hap_map[r] == lbl}
                    consensus = plurality_consensus(hap_seqs)
                    consensus_out.write(
                        f">group_{gid}_hap_{lbl} n_reads={len(hap_seqs)}{header_tags(gid)}\n"
                        f"{consensus}\n")
                group_meta_out.write(f"{gid}\t{len(msa_seqs)}\tTRUE\n")
                continue

            consensus = plurality_consensus(msa_seqs)
            consensus_out.write(
                f">group_{gid} n_reads={len(msa_seqs)}{header_tags(gid)}\n{consensus}\n")
            group_meta_out.write(f"{gid}\t{len(msa_seqs)}\tTRUE\n")

# ── Write sub-haplotype assignments ──────────────────────────────────────────
read_to_group = {r: gid for gid, rlist in groups.items() for r in rlist}
with open("~{out_prefix}_subhaps.tsv", 'w') as out:
    out.write("read_name\tgroup_id\tsub_haplotype\n")
    for r, hap in sorted(subhap_assignments.items()):
        out.write(f"{r}\t{read_to_group.get(r, 'NA')}\t{hap}\n")

# Tarball per-group MSA files (for diagnostic download)
with tarfile.open("~{out_prefix}_all_msa.tar.gz", "w:gz") as tar:
    tar.add("msa_groups", arcname="msa_groups")

n_multi = sum(1 for gid in set(read_to_group[r] for r in subhap_assignments)
              if len({subhap_assignments[r]
                      for r in subhap_assignments if read_to_group.get(r) == gid}) > 1)
print(f"Consensus build complete; {len(subhap_assignments)} reads assigned sub-haplotypes, "
      f"{n_multi} group(s) split by haplotype")
PYEOF
  >>>

  output {
    File all_msa_tarball  = "~{out_prefix}_all_msa.tar.gz"
    File consensus_fasta  = "~{out_prefix}_consensus.fasta"
    File group_meta_tsv   = "~{out_prefix}_group_meta.tsv"
    File subhaps_tsv      = "~{out_prefix}_subhaps.tsv"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87"
    cpu: 1
    memory: "8 GiB"
    disk: "50 GB"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
  }
}


task categorize_reads {
  meta {
    description: "Assign final edit-outcome category to each read using structural groups, sub-haplotypes, and Pass-2 blastn evidence"
  }

  parameter_meta {
    filtered_reads_fasta: "Original (unclipped) filtered reads FASTA — used to enumerate all reads"
    groups_tsv: "Structural group assignments (from structural_grouping)"
    subhaps_tsv: "Sub-haplotype assignments (from build_group_consensus)"
    pass2_alignment_tsv: "Pass-2 blastn tabular output run on unclipped filtered reads"
    gap_merge_threshold: "Max gap (bp) between abutting HDR parts to call HDR_candidate"
    wt_gap_threshold: "Max gap (bp) between left_HA and right_HA to call WT_candidate"
    min_flank_bp: "Minimum bases flanking all part hits required to call HDR or WT (vs inconclusive)"
    min_part_hit_bp: "Minimum query-span (bp) for a blastn hit to be counted; filters repetitive short matches"
  }

  input {
    File filtered_reads_fasta
    File groups_tsv
    File subhaps_tsv
    File pass2_alignment_tsv
    String sample_id
    Int gap_merge_threshold = 5
    Int wt_gap_threshold    = 20
    Int min_flank_bp = 50
    Int min_part_hit_bp = 100
    RuntimeAttributes runtime_attributes
  }

  String out_prefix = "~{sample_id}_edit_categories"

  command <<<
    set -euxo pipefail

    python3 <<'PYEOF'
from collections import defaultdict, Counter

HDR_THR = ~{gap_merge_threshold}
WT_THR  = ~{wt_gap_threshold}
MIN_HIT = ~{min_part_hit_bp}

# HDR part order constants (used for truncated-read detection)
_FWD_HDR = [('left_HA', '+'), ('payload', '+'), ('right_HA', '+')]
_REV_HDR = [('right_HA', '-'), ('payload', '-'), ('left_HA', '-')]

def is_partial_hdr(pairs, ordered):
    """Return True iff pairs is a non-empty strict contiguous proper subset of ordered."""
    n = len(pairs)
    if n == 0 or n >= len(ordered):
        return False
    for start in range(len(ordered) - n + 1):
        if list(ordered[start:start + n]) == pairs:
            return True
    return False

# ── Structural category from ordered part hits ────────────────────────────────
# HiFi reads can be in either orientation. A reverse-strand HDR read has all
# three parts on '-' in t_start order: right_HA:- → payload:- → left_HA:-.
# WT on '-': right_HA:- → left_HA:-.
# Mixed strands → misintegration.
def structural_category(ps, gs):
    """ps = [(part, strand), ...] sorted by t_start; gs = gaps between them."""
    if not ps:
        return 'no_hit'
    parts_list   = [p for p, s in ps]
    strands_list = [s for p, s in ps]

    if len(set(strands_list)) != 1:
        return 'misintegration'

    strand = strands_list[0]
    if strand == '+':
        if parts_list == ['left_HA', 'payload', 'right_HA'] and all(abs(g) <= HDR_THR for g in gs):
            return 'HDR_candidate'
        if parts_list == ['left_HA', 'right_HA'] and len(gs) == 1 and abs(gs[0]) <= WT_THR:
            return 'WT_candidate'
        if parts_list in (['left_HA', 'payload'], ['payload', 'right_HA']) and all(abs(g) <= HDR_THR for g in gs):
            return 'partial'
    else:
        if parts_list == ['right_HA', 'payload', 'left_HA'] and all(abs(g) <= HDR_THR for g in gs):
            return 'HDR_candidate'
        if parts_list == ['right_HA', 'left_HA'] and len(gs) == 1 and abs(gs[0]) <= WT_THR:
            return 'WT_candidate'
        if parts_list in (['payload', 'left_HA'], ['right_HA', 'payload']) and all(abs(g) <= HDR_THR for g in gs):
            return 'partial'

    return 'misintegration'

# ── Load structural groups ────────────────────────────────────────────────────
group_info = {}   # read_name -> {gid, sig, size}
with open("~{groups_tsv}") as fh:
    fh.readline()
    for ln in fh:
        f = ln.rstrip('\n').split('\t')
        if len(f) < 4:
            continue
        group_info[f[0]] = {'gid': int(f[1]), 'sig': f[2], 'size': int(f[3])}

# ── Load sub-haplotype assignments ────────────────────────────────────────────
subhap_info = {}
with open("~{subhaps_tsv}") as fh:
    fh.readline()
    for ln in fh:
        f = ln.rstrip('\n').split('\t')
        if len(f) >= 3:
            subhap_info[f[0]] = f[2]

# ── Load Pass-2 blastn hits ───────────────────────────────────────────────────
# Collect both span tuples (for flanking) and ordered (part, strand) lists (for category).
# Two-pass: parse all hits first, then apply proximity-aware size filter.
# A hit is filtered out only if ALL of:
#   (a) hit query-span < max(MIN_HIT, 0.40 * part_len)
#   (b) not within WT_THR bp of any other hit on the same read
pass2_hits    = defaultdict(list)   # read_name -> [(t_start, t_end)]
hits_by_read  = defaultdict(list)   # read_name -> [(part, strand, t_start, t_end)]
read_lens     = {}

_raw_hits = defaultdict(list)  # read_name -> [(part, strand, t_start, t_end, hit_size, qlen)]
_raw_slen = {}                 # read_name -> slen

with open("~{pass2_alignment_tsv}") as fh:
    for ln in fh:
        if not ln.strip() or ln.startswith('#'):
            continue
        f = ln.rstrip('\n').split('\t')
        if len(f) < 15:
            continue
        qname = f[0]
        if not (qname.endswith('_left_HA') or qname.endswith('_right_HA')
                or qname.endswith('_payload')):
            continue
        read_name = f[4]
        qlen      = int(f[1])
        hit_size  = abs(int(f[3]) - int(f[2]))
        sstart, send = int(f[6]), int(f[7])
        t_start = min(sstart, send) - 1
        t_end   = max(sstart, send)
        strand  = '+' if f[13] == 'plus' else '-'
        if   qname.endswith('_left_HA'):  part = 'left_HA'
        elif qname.endswith('_right_HA'): part = 'right_HA'
        else:                             part = 'payload'
        _raw_hits[read_name].append((part, strand, t_start, t_end, hit_size, qlen))
        _raw_slen[read_name] = int(f[5])

def _keep_hit(hit, others):
    """Return True if this hit should be kept (passes size OR is near a size-passing hit)."""
    part, strand, t_start, t_end, hit_size, qlen = hit
    min_size = max(MIN_HIT, int(0.40 * qlen))
    if hit_size >= min_size:
        return True
    # Keep if within WT_THR bp of another hit that itself passes the size threshold
    for o in others:
        o_hit_size, o_qlen = o[4], o[5]
        o_min_size = max(MIN_HIT, int(0.40 * o_qlen))
        if o_hit_size < o_min_size:
            continue
        o_start, o_end = o[2], o[3]
        gap = max(t_start, o_start) - min(t_end, o_end)  # <=0 means overlap
        if gap <= WT_THR:
            return True
    return False

for read_name, raw in _raw_hits.items():
    others_by_idx = [raw[:i] + raw[i+1:] for i in range(len(raw))]
    for hit, others in zip(raw, others_by_idx):
        if not _keep_hit(hit, others):
            continue
        part, strand, t_start, t_end, hit_size, qlen = hit
        pass2_hits[read_name].append((t_start, t_end))
        hits_by_read[read_name].append((part, strand, t_start, t_end))
    read_lens[read_name] = _raw_slen[read_name]

def get_ps_gs(read_name):
    hits = sorted(hits_by_read[read_name], key=lambda h: h[2])
    ps   = [(h[0], h[1]) for h in hits]
    gs   = [hits[i+1][2] - hits[i][3] for i in range(len(hits)-1)]
    return ps, gs

# ── Enumerate all reads from FASTA, collecting lengths ───────────────────────
all_reads  = []
fasta_lens = {}   # read_name -> sequence length
with open("~{filtered_reads_fasta}") as fh:
    cur_name = None
    cur_len  = 0
    for ln in fh:
        if ln.startswith('>'):
            if cur_name is not None:
                fasta_lens[cur_name] = cur_len
            cur_name = ln[1:].split()[0]
            cur_len  = 0
            all_reads.append(cur_name)
        else:
            cur_len += len(ln.rstrip('\n'))
    if cur_name is not None:
        fasta_lens[cur_name] = cur_len

MIN_FLANK = ~{min_flank_bp}

def flanking(read_name):
    hits = pass2_hits.get(read_name, [])
    rlen = read_lens.get(read_name, 0)
    if not hits:
        return 0, 0
    return min(h[0] for h in hits), rlen - max(h[1] for h in hits)

def categorize(read_name):
    info      = group_info.get(read_name)
    lf, rf    = flanking(read_name)
    has_flank = lf >= MIN_FLANK or rf >= MIN_FLANK

    if not hits_by_read.get(read_name):
        return 'inconclusive', 'no_parts_hit', ''

    ps, gs = get_ps_gs(read_name)
    cat    = structural_category(ps, gs)
    sig    = info['sig'] if info else '|'.join(f"{p}:{s}" for p, s in ps)
    gid    = info['gid']  if info else 'NA'
    size   = info['size'] if info else 1

    if cat == 'no_hit':
        return 'inconclusive', 'no_parts_hit', ''

    if cat == 'HDR_candidate':
        if not has_flank:
            return 'inconclusive', f'HDR_sig_no_flank|{sig}', ''
        subhap   = subhap_info.get(read_name, f'group_{gid}_hap_A')
        evidence = f'singleton|{sig}' if size == 1 else sig
        return 'HDR', evidence, subhap

    if cat == 'WT_candidate':
        if not has_flank:
            return 'inconclusive', f'WT_sig_no_flank|{sig}', ''
        subhap = subhap_info.get(read_name, f'group_{gid}_hap_A')
        return 'WT', sig, subhap

    # partial or misintegration
    # A partial signature where the read is truncated on the missing-part side
    # is inconclusive — the read simply ends before the edit is complete.
    if lf < MIN_FLANK or rf < MIN_FLANK:
        if is_partial_hdr(ps, _FWD_HDR) or is_partial_hdr(ps, _REV_HDR):
            return 'inconclusive', f'truncated_read|{sig}', ''
    return 'misintegration', sig, ''

# ── Write output ──────────────────────────────────────────────────────────────
with open("~{out_prefix}.tsv", 'w') as out:
    out.write("read_name\tcategory\tevidence\tsub_haplotype\t"
              "group_id\tgroup_size\tread_length\tleft_flank_bp\tright_flank_bp\n")
    cats = []
    for rname in all_reads:
        cat, evidence, subhap = categorize(rname)
        info  = group_info.get(rname)
        gid   = info['gid']  if info else 'NA'
        gsize = info['size'] if info else 0
        rlen  = fasta_lens.get(rname, 0)
        lf, rf = flanking(rname)
        out.write(f"{rname}\t{cat}\t{evidence}\t{subhap}\t{gid}\t{gsize}\t{rlen}\t{lf}\t{rf}\n")
        cats.append(cat)

counts = Counter(cats)
print(f"Total reads : {len(all_reads)}")
for cat, n in sorted(counts.items()):
    print(f"  {cat}: {n}")
PYEOF
  >>>

  output {
    File categories_tsv = "~{out_prefix}.tsv"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87"
    cpu: 1
    memory: "4 GiB"
    disk: "10 GB"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
  }
}


task make_wide_coverage_table {
  meta {
    description: "Build wide per-part coverage table from blastn part alignments (backwards-compatible replacement for parse_edit_parts)"
  }

  parameter_meta {
    parts_alignment_tsv: "Blastn tabular output from blastn_remap_parts"
    coverage_threshold: "Minimum per-hit coverage (n_match / query_len) for faithful alignments"
    min_part_hit_bp: "Minimum query-span (bp) for a blastn hit to be counted; filters repetitive short matches"
  }

  input {
    File parts_alignment_tsv
    Float coverage_threshold = 0.8
    Int min_part_hit_bp = 100
    String sample_id
    RuntimeAttributes runtime_attributes
  }

  String out_prefix = "~{sample_id}_edit_parts_coverage"

  command <<<
    set -euxo pipefail

    python3 <<'PYEOF'
import csv
from collections import defaultdict
from typing import Dict, List, Tuple

# Hit tuple: (q_start, q_end, t_start, t_end, strand, n_match)
Hit     = Tuple[int, int, int, int, str, int]
HitDict = Dict[str, Dict[str, List[Hit]]]
LenDict = Dict[str, int]

# ── blastn tabular columns ────────────────────────────────────────────────────
# 0 qseqid  1 qlen  2 qstart  3 qend  4 sseqid  5 slen  6 sstart  7 send
# 8 nident  9 mismatch  10 gapopen  11 gaps  12 length  13 sstrand  14 bitscore

PART_SUFFIXES = ('_left_HA', '_payload', '_right_HA')
MIN_HIT       = ~{min_part_hit_bp}
COV_THR       = ~{coverage_threshold}

faithful: HitDict = defaultdict(lambda: defaultdict(list))
tlen: LenDict = {}
qlen: LenDict = {}

with open("~{parts_alignment_tsv}") as fh:
    for ln in fh:
        ln = ln.rstrip("\n")
        if not ln or ln.startswith("#"):
            continue
        f = ln.split("\t")
        if len(f) < 15:
            continue
        qname = f[0]
        if not any(qname.endswith(s) for s in PART_SUFFIXES):
            continue
        if abs(int(f[3]) - int(f[2])) < MIN_HIT:
            continue
        n_match = int(f[8])
        n_gaps  = int(f[11])
        q_len   = int(f[1])
        cov_hit = n_match / q_len
        gap_ok  = n_gaps <= 0.05 * n_match
        if cov_hit < COV_THR or not gap_ok:
            continue
        sstart, send = int(f[6]), int(f[7])
        tname = f[4]
        tlen.setdefault(tname, int(f[5]))
        qlen.setdefault(qname, q_len)
        hit: Hit = (
            int(f[2]) - 1, int(f[3]),
            min(sstart, send) - 1, max(sstart, send),
            "+" if f[13] == "plus" else "-",
            n_match,
        )
        faithful[tname][qname].append(hit)

def _agg_cov(hit_list: List[Hit], q_len: int) -> float:
    return min(1.0, sum(h[5] for h in hit_list) / q_len) if q_len else 0.0

all_q  = sorted({q for t in faithful.values() for q in t})
header = ["target_name", "target_len", "all_one_hit"]
for q in all_q:
    header.extend([f"{q}_len", f"{q}_cov", f"{q}_num_hits", f"{q}_hit_coordinates"])

with open("~{out_prefix}_wide_faithful.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(header)
    for tgt in sorted(faithful):
        qdict    = faithful[tgt]
        num_hits = [len(qdict[q]) if q in qdict else 0 for q in all_q]
        row      = [tgt, tlen[tgt], "TRUE" if all(n == 1 for n in num_hits) else "FALSE"]
        for q, n in zip(all_q, num_hits):
            row.append(qlen[q])
            if n:
                cov      = _agg_cov(qdict[q], qlen[q])
                hits_str = ";".join(f"{qs}-{qe}:{ts}-{te}:{s}"
                                    for qs, qe, ts, te, s, _ in qdict[q])
                row.extend([f"{cov:.3f}", str(n), hits_str])
            else:
                row.extend(["", "", ""])
        w.writerow(row)

print(f"Wide faithful table: {len(faithful)} reads, {len(all_q)} parts")
PYEOF
  >>>

  output {
    File wide_faithful_table = "~{out_prefix}_wide_faithful.tsv"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87"
    cpu: 1
    memory: "4 GiB"
    disk: "10 GB"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
  }
}

task align_consensus_msa {
  meta {
    description: "Align child + parent consensus sequences plus the expected HDR template into a single MAFFT alignment. Reverse-strand consensus sequences are reverse-complemented so all inputs share a common orientation. Group consensus records tagged category=partial are omitted (partial-part signatures do not represent a complete observation of an allele)."
  }

  parameter_meta {
    child_consensus_fasta:  "Consensus FASTA from child sample (build_group_consensus); headers must contain strand=+/- and category=HDR|WT|partial|misintegration"
    parent_consensus_fasta: "Consensus FASTA from parent sample (build_group_consensus); headers must contain strand=+/- and category=..."
    child_sample_id:        "Sample ID used to prefix child sequence names in the MSA"
    parent_sample_id:       "Sample ID used to prefix parent sequence names in the MSA"
    crispr_edit_json:       "CRISPR edit JSON; left_ha + payload + right_ha are concatenated and included in the MSA as the expected HDR template"
  }

  input {
    File   child_consensus_fasta
    File   parent_consensus_fasta
    String child_sample_id
    String parent_sample_id
    File   crispr_edit_json
    RuntimeAttributes runtime_attributes
  }

  String out_prefix = "~{child_sample_id}_vs_~{parent_sample_id}_consensus"

  command <<<
    set -euxo pipefail

    OUT="~{out_prefix}"
    INPUT="${OUT}_input.fasta"

    # Extract left_ha/payload/right_ha from crispr_edit_json (DNA strings, no embedded quotes)
    extract_field() {
      grep -oE "\"$1\"[[:space:]]*:[[:space:]]*\"[^\"]*\"" "~{crispr_edit_json}" \
        | head -1 \
        | sed -E 's/^[^:]*:[[:space:]]*"([^"]*)"$/\1/'
    }
    LEFT_HA=$(extract_field left_ha)
    PAYLOAD=$(extract_field payload)
    RIGHT_HA=$(extract_field right_ha)
    HDR_TEMPLATE=$(printf '%s%s%s' "$LEFT_HA" "$PAYLOAD" "$RIGHT_HA" | tr 'acgtn' 'ACGTN')

    : > "$INPUT"
    if [[ -n "$HDR_TEMPLATE" ]]; then
      printf ">HDR_expected\n%s\n" "$HDR_TEMPLATE" >> "$INPUT"
    fi
    echo "HDR template: ${#HDR_TEMPLATE} bp" >&2

    # Parse a consensus FASTA: skip category=partial, RC strand=-, prefix names with label.
    parse_consensus() {
      local fasta=$1 label=$2
      awk -v label="$label" '
        function rc_char(c) {
          if (c=="A") return "T"; if (c=="T") return "A"
          if (c=="C") return "G"; if (c=="G") return "C"
          return "N"
        }
        function rc(s,   i,r) {
          r=""
          for (i=length(s); i>=1; i--) r = r rc_char(substr(s,i,1))
          return r
        }
        function flush(   out) {
          if (name == "" || seq == "") return
          if (cat == "partial") {
            print "  skip partial: " label "_" name > "/dev/stderr"
            return
          }
          out = toupper(seq)
          if (strand == "-") out = rc(out)
          print ">" label "_" name
          print out
          kept++
        }
        /^>/ {
          flush()
          name=""; strand=""; cat=""; seq=""
          header=substr($0, 2)
          if (match(header, /group_[0-9]+(_hap_[A-Za-z0-9_]+)?(_singleton)?/))
            name = substr(header, RSTART, RLENGTH)
          if (match(header, /strand=[+-]/))
            strand = substr(header, RSTART+7, 1)
          if (match(header, /category=[A-Za-z_]+/))
            cat = substr(header, RSTART+9, RLENGTH-9)
          if (name == "" || strand == "") name = ""
          next
        }
        { seq = seq $0 }
        END {
          flush()
          print label ": " kept+0 " non-partial records" > "/dev/stderr"
        }
      ' "$fasta"
    }

    parse_consensus "~{parent_consensus_fasta}" "~{parent_sample_id}" >> "$INPUT"
    parse_consensus "~{child_consensus_fasta}"  "~{child_sample_id}"  >> "$INPUT"

    OUTPUT_ALN="${OUT}.aln"

    nseq=$(grep -c '^>' "$INPUT" || true)

    if [[ "$nseq" -eq 0 ]]; then
      echo "No sequences to align" >&2
      printf "CLUSTAL W\n\n" > "$OUTPUT_ALN"
    elif [[ "$nseq" -eq 1 ]]; then
      name=$(grep '^>' "$INPUT" | head -1 | sed 's/^>//')
      seq=$(grep -v '^>' "$INPUT" | tr -d '\n')
      stars=$(printf '%*s' "${#seq}" '' | tr ' ' '*')
      printf "CLUSTAL W\n\n%-30s %s\n%-30s %s\n" \
        "${name:0:30}" "$seq" "" "$stars" > "$OUTPUT_ALN"
    else
      mafft --auto --quiet --clustalout "$INPUT" > "$OUTPUT_ALN"
    fi
    echo "${nseq} sequences → ${OUTPUT_ALN}"
  >>>

  output {
    File? consensus_msa_clustal = "~{out_prefix}.aln"
  }

  runtime {
    docker: "quay.io/biocontainers/mafft@sha256:b8ccc402a9156868fefda5caa8693f847b47dcc69aae848557f7115850e60790"
    cpu: 1
    memory: "8 GiB"
    disk: "20 GB"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
  }
}
