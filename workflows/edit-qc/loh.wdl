version 1.1

import "../wdl-common/wdl/structs.wdl"

# Task to check for loss of heterozygosity (LOH) around a CRISPR cut site
task check_loh {
  meta {
    description: "Check for loss of heterozygosity around a CRISPR cut site using joint-called small variant VCF: count parental HET SNPs retained as HET in the edited sample"
  }

  parameter_meta {
    sample_id: "Sample identifier (for output file naming)"
    parent_sample_idx: "0-based index of the parent sample in joint_vcf"
    sample_idx: "0-based index of the edited sample in joint_vcf"
    joint_vcf: "Joint-called, phased small variant VCF containing both parent and sample"
    joint_vcf_index: "Tabix index for joint_vcf"
    expected_cut_site: "Expected cut site in chr:pos format (e.g. chr15:65869582)"
    cut_strand: "Strand of the guide at the cut site (+ or -)"
    window_size: "Window size in bp upstream and downstream of the cut site (default 200000)"
  }

  input {
    String sample_id
    Int parent_sample_idx
    Int sample_idx
    File joint_vcf
    File joint_vcf_index
    String expected_cut_site
    String cut_strand
    Int window_size = 200000
    Int mem_gb = 8
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(joint_vcf, "GB")) + 10

  command <<<
    set -euo pipefail

    python3 <<'PYEOF'
import subprocess
import sys

sample_id         = "~{sample_id}"
parent_sample_idx = ~{parent_sample_idx}
sample_idx        = ~{sample_idx}
joint_vcf         = "~{joint_vcf}"
cut_site          = "~{expected_cut_site}"
cut_strand        = "~{cut_strand}"
window_size       = ~{window_size}

# ---------------------------------------------------------------------------
# Parse cut site and compute upstream / downstream windows
# upstream:   opposite direction of cut strand
#   + strand: [cut_pos - window, cut_pos]
#   - strand: [cut_pos, cut_pos + window]
# downstream: further along in direction of cut strand
#   + strand: [cut_pos, cut_pos + window]
#   - strand: [cut_pos - window, cut_pos]
# ---------------------------------------------------------------------------
chrom, pos_str = cut_site.rsplit(":", 1)
cut_pos = int(pos_str)

if cut_strand == "+":
    up_start = max(1, cut_pos - window_size)
    up_end   = cut_pos
    dn_start = cut_pos
    dn_end   = cut_pos + window_size
else:
    up_start = cut_pos
    up_end   = cut_pos + window_size
    dn_start = max(1, cut_pos - window_size)
    dn_end   = cut_pos

upstream_region   = f"{chrom}:{up_start}-{up_end}"
downstream_region = f"{chrom}:{dn_start}-{dn_end}"

print(f"Cut site:          {cut_site} (strand {cut_strand})", file=sys.stderr)
print(f"Upstream region:   {upstream_region}", file=sys.stderr)
print(f"Downstream region: {downstream_region}", file=sys.stderr)


def bcftools_count(region, filter_expr):
    """Count variants in region passing filter_expr (bcftools -i expression)."""
    r = subprocess.run(
        [
            "bcftools", "view",
            "-r", region,
            "-i", filter_expr,
            "--output-type", "v",
            joint_vcf,
        ],
        capture_output=True, check=True,
    )
    return sum(1 for ln in r.stdout.decode().splitlines() if ln and not ln.startswith("#"))


def count_for_region(region, label):
    # Denominator: variants where parent is het
    denom = bcftools_count(region, f"GT[{parent_sample_idx}]='het'")
    # Numerator: variants where both parent AND sample are het
    numer = bcftools_count(region, f"GT[{parent_sample_idx}]='het' && GT[{sample_idx}]='het'")
    print(f"{label}: parent_het={denom}, both_het={numer}", file=sys.stderr)
    return numer, denom


up_num, up_den = count_for_region(upstream_region,   "upstream")
dn_num, dn_den = count_for_region(downstream_region, "downstream")

up_ratio_str = f"{up_num}/{up_den}"
dn_ratio_str = f"{dn_num}/{dn_den}"
up_pct = (up_num / up_den) if up_den > 0 else 0.0
dn_pct = (dn_num / dn_den) if dn_den > 0 else 0.0

print(f"Upstream   het intact: {up_ratio_str} = {up_pct:.4f}", file=sys.stderr)
print(f"Downstream het intact: {dn_ratio_str} = {dn_pct:.4f}", file=sys.stderr)

with open(f"{sample_id}_loh_het_intact_upstream_ratio.txt",   "w") as fh:
    fh.write(up_ratio_str + "\n")
with open(f"{sample_id}_loh_het_intact_downstream_ratio.txt", "w") as fh:
    fh.write(dn_ratio_str + "\n")
with open(f"{sample_id}_loh_het_intact_upstream_pct.txt",     "w") as fh:
    fh.write(f"{up_pct:.6f}\n")
with open(f"{sample_id}_loh_het_intact_downstream_pct.txt",   "w") as fh:
    fh.write(f"{dn_pct:.6f}\n")

with open(f"{sample_id}_loh_summary.tsv", "w") as fh:
    fh.write("field\tvalue\n")
    fh.write(f"upstream_region\t{upstream_region}\n")
    fh.write(f"downstream_region\t{downstream_region}\n")
    fh.write(f"upstream_parent_het_total\t{up_den}\n")
    fh.write(f"upstream_both_het\t{up_num}\n")
    fh.write(f"upstream_het_intact_ratio\t{up_ratio_str}\n")
    fh.write(f"upstream_het_intact_pct\t{up_pct:.6f}\n")
    fh.write(f"downstream_parent_het_total\t{dn_den}\n")
    fh.write(f"downstream_both_het\t{dn_num}\n")
    fh.write(f"downstream_het_intact_ratio\t{dn_ratio_str}\n")
    fh.write(f"downstream_het_intact_pct\t{dn_pct:.6f}\n")
PYEOF
  >>>

  output {
    String het_intact_upstream_ratio   = read_string("~{sample_id}_loh_het_intact_upstream_ratio.txt")
    String het_intact_downstream_ratio = read_string("~{sample_id}_loh_het_intact_downstream_ratio.txt")
    Float  het_intact_upstream_pct     = read_float("~{sample_id}_loh_het_intact_upstream_pct.txt")
    Float  het_intact_downstream_pct   = read_float("~{sample_id}_loh_het_intact_downstream_pct.txt")
    File   loh_summary                 = "~{sample_id}_loh_summary.tsv"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87"
    cpu: 2
    memory: "~{mem_gb} GiB"
    disk: disk_size + " GB"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
  }
}
