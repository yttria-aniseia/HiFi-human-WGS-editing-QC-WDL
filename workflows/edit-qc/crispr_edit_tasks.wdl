version 1.1

import "../wdl-common/wdl/structs.wdl"

task create_parts_query_fasta {
  meta {
    description: "Create FASTA file with edit parts (left_ha, payload, right_ha) from CRISPR edit JSON"
  }

  parameter_meta {
    crispr_edit_json: "JSON file describing CRISPR edit with left_ha, payload, right_ha sequences"
  }

  input {
    File crispr_edit_json
    RuntimeAttributes runtime_attributes
  }

  Int threads = 1
  String out_prefix = basename(crispr_edit_json, ".json")

  command <<<
    set -euxo pipefail

    python3 <<EOF
import json
import sys

# Read CRISPR edit JSON
with open("~{crispr_edit_json}", 'r') as f:
    crispr_data = json.load(f)

# Create parts query FASTA
with open("~{out_prefix}_parts_query.fasta", 'w') as fasta_out:
    for edit in crispr_data['edits']:
        edit_id = edit['id']
        left_ha = edit['edit']['left_ha']
        payload = edit['edit']['payload']  
        right_ha = edit['edit']['right_ha']
        
        # Write each part as separate FASTA entry
        fasta_out.write(f">{edit_id}_left_HA\n{left_ha}\n")
        fasta_out.write(f">{edit_id}_payload\n{payload}\n")
        fasta_out.write(f">{edit_id}_right_HA\n{right_ha}\n")
        fasta_out.write(f">{edit_id}_HDR_template\n{left_ha}{payload}{right_ha}\n")
        fasta_out.write(f">{edit_id}_WT_template\n{left_ha}{right_ha}\n")
EOF
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
    description: "Use minimap2 to align edit parts to reads"
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

task samtools_filter_reads {
  meta {
    description: "Use samtools to filter reads based on alignment results"
  }

  parameter_meta {
    fasta_reads: "FASTA converted reads"
    initial_sam: "SAM file from minimap2 alignment"
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
    awk 'NR==FNR {h[$1]; next} /^>/ {p = substr($0, 2) in h} p' matched_reads_uniq.txt ~{fasta_reads}  > ~{out_prefix}_filtered_reads.fasta

    # Count reads
    TOTAL_READS=$(grep -c '^>' ~{fasta_reads})
    MATCHED_READS=$(grep -c '^>' ~{out_prefix}_filtered_reads.fasta || echo 0)
    echo "Total reads in input: $TOTAL_READS"
    echo "Reads matching edit parts: $MATCHED_READS"
  >>>

  output {
    File filtered_reads_fasta = "~{out_prefix}_filtered_reads.fasta"
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

task minimap2_remap_parts {
  meta {
    description: "Remap edit parts to filtered reads with relaxed parameters"
  }

  parameter_meta {
    filtered_reads_fasta: "FASTA with reads containing edit parts"
    parts_query_fasta: "FASTA file with edit parts sequences"
  }

  input {
    File filtered_reads_fasta
    File parts_query_fasta
    String sample_id
    RuntimeAttributes runtime_attributes
  }

  String out_prefix = "~{sample_id}_edit_parts"

  command <<<
    set -euxo pipefail

    minimap2 -cx asm10 \
        -N 1000 \
        -p 0.5 \
        ~{filtered_reads_fasta} ~{parts_query_fasta} \
        > ~{out_prefix}_minimap2_crispr_parts.paf

    # Count alignments
    echo "Total alignments: $(wc -l < ~{out_prefix}_minimap2_crispr_parts.paf)"
  >>>

  output {
    File parts_alignment_paf = "~{out_prefix}_minimap2_crispr_parts.paf"
  }

  runtime {
    docker: "quay.io/biocontainers/minimap2@sha256:fdc9ef8bfbd31bab59a61b2e90dd226647deed971556175a4cd004f0bcdc7608"
    cpu: 4
    memory: "8 GB"
    disk: "20 GB"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
  }
}

task parse_edit_parts {
  meta {
    description: "Parse PAF alignments to create summary tables with coverage metrics"
  }

  parameter_meta {
    parts_alignment_paf: "PAF file from minimap2 alignment of edit parts to reads"
    coverage_threshold: "Minimum coverage threshold for faithful alignments"
  }

  input {
    File parts_alignment_paf
    Float coverage_threshold = 0.8
    String sample_id
    RuntimeAttributes runtime_attributes
  }

  String out_prefix = "~{sample_id}_edit_parts_summary"

  command <<<
    set -euxo pipefail

    cat > parse_paf_to_table.py <<'EOF'
#!/usr/bin/env python3
"""
parse_paf_to_table.py  –  2025-07-03

Adds per-query coverage columns computed as ∑n_match / query_len
----------------------------------------------------------------

Long table header:
  target_name target_len query_name query_len coverage num_hits hits

Wide table header:
  target_name target_len all_one_hit
              <q1>_len <q1>_cov <q1>_num_hits <q1>_hit_coordinates
              <q2>_len <q2>_cov <q2>_num_hits <q2>_hit_coordinates ...
"""

import argparse, csv
from collections import defaultdict
from typing import Dict, List, Tuple

# ──────────── type aliases ──────────── #
Hit      = Tuple[int, int, int, int, str, int]          # … + n_match
HitDict  = Dict[str, Dict[str, List[Hit]]]              # hits[target][query] = [Hit, …]
LenDict  = Dict[str, int]

# ──────────── parsing ──────────── #
def parse_paf_line(line: str) -> Dict:
    f = line.rstrip("\n").split("\t")
    if len(f) < 12:
        raise ValueError("PAF line has fewer than 12 mandatory columns")
    return dict(
        query_name=f[0],  query_len=int(f[1]),
        q_start=int(f[2]), q_end=int(f[3]), strand=f[4],
        target_name=f[5], target_len=int(f[6]),
        t_start=int(f[7]), t_end=int(f[8]),
        n_match=int(f[9]), n_align=int(f[10]), mapq=int(f[11]),
    )

def load_hits(paf_path: str, cov_thr: float):
    hits, faithful = defaultdict(lambda: defaultdict(list)), defaultdict(lambda: defaultdict(list))
    tlen, qlen = {}, {}
    with open(paf_path) as fh:
        for ln in fh:
            if not ln.strip(): continue
            r = parse_paf_line(ln)
            tlen.setdefault(r["target_name"], r["target_len"])
            qlen.setdefault(r["query_name"],  r["query_len"])
            hit: Hit = (
                r["q_start"], r["q_end"],
                r["t_start"], r["t_end"],
                r["strand"],  r["n_match"]
            )
            hits[r["target_name"]][r["query_name"]].append(hit)

            cov_hit = r["n_match"] / r["query_len"]
            gap_ok  = abs(r["n_match"] - r["n_align"]) <= 0.05 * r["n_match"]
            if cov_hit >= cov_thr and gap_ok:
                faithful[r["target_name"]][r["query_name"]].append(hit)
    return hits, faithful, tlen, qlen

# ──────────── writers ──────────── #
def _aggregate_coverage(hit_list: List[Hit], q_len: int) -> float:
    total_nm = sum(h[5] for h in hit_list)
    return min(1.0, total_nm / q_len) if q_len else 0.0

def write_long(h: HitDict, tlen: LenDict, qlen: LenDict, path: str):
    with open(path,"w",newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["target_name","target_len",
                    "query_name","query_len","coverage",
                    "num_hits","hits"])
        for tgt in sorted(h):
            for qry in sorted(h[tgt]):
                hl   = h[tgt][qry]
                cov  = _aggregate_coverage(hl, qlen[qry])
                hits_str = ";".join(f"{qs}-{qe}:{ts}-{te}:{s}"
                                    for qs,qe,ts,te,s,_ in hl)
                w.writerow([tgt, tlen[tgt],
                            qry, qlen[qry], f"{cov:.3f}",
                            len(hl), hits_str])

def write_wide(h: HitDict, tlen: LenDict, qlen: LenDict, path: str):
    all_q = sorted({q for t in h.values() for q in t})
    header=["target_name","target_len","all_one_hit"]
    for q in all_q:
        header.extend([f"{q}_len", f"{q}_cov", f"{q}_num_hits", f"{q}_hit_coordinates"])
    with open(path,"w",newline="") as f:
        w=csv.writer(f,delimiter="\t"); w.writerow(header)
        for tgt in sorted(h):
            qdict=h[tgt]
            num_hits=[len(qdict[q]) if q in qdict else 0 for q in all_q]
            all_one=all(n==1 for n in num_hits)
            row=[tgt, tlen[tgt], "TRUE" if all_one else "FALSE"]
            for q,n in zip(all_q,num_hits):
                row.append(qlen[q])                 # <query>_len
                if n:
                    cov=_aggregate_coverage(qdict[q], qlen[q])
                    hits_str=";".join(f"{qs}-{qe}:{ts}-{te}:{s}"
                                      for qs,qe,ts,te,s,_ in qdict[q])
                    row.extend([f"{cov:.3f}", str(n), hits_str])
                else:
                    row.extend(["", "", ""])        # cov, num_hits, coords blank
            w.writerow(row)

# ──────────── CLI ──────────── #
def main():
    ap=argparse.ArgumentParser(description="PAF summary tables with coverage")
    ap.add_argument("paf"); ap.add_argument("-o","--output",default="paf_target_summary.tsv")
    ap.add_argument("-c","--coverage_threshold",type=float,default=0.8)
    a=ap.parse_args()
    hits,faithful,tlen,qlen=load_hits(a.paf,a.coverage_threshold)
    write_long(hits,tlen,qlen,a.output)
    write_long(faithful,tlen,qlen,a.output.replace(".tsv","_faithful.tsv"))
    write_wide(hits,tlen,qlen,a.output.replace(".tsv","_wide.tsv"))
    write_wide(faithful,tlen,qlen,a.output.replace(".tsv","_wide_faithful.tsv"))
    print("✓ tables with coverage written")

if __name__=="__main__":
    main()
EOF

    chmod +x parse_paf_to_table.py

    # Run the parser
    python3 parse_paf_to_table.py \
      ~{parts_alignment_paf} \
      -o ~{out_prefix}.tsv \
      -c ~{coverage_threshold}
  >>>

  output {
    File summary_table = "~{out_prefix}.tsv"
    File faithful_table = "~{out_prefix}_faithful.tsv" 
    File wide_table = "~{out_prefix}_wide.tsv"
    File wide_faithful_table = "~{out_prefix}_wide_faithful.tsv"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87"
    cpu: 2
    memory: "4 GB"
    disk: "10 GB"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
  }
}
