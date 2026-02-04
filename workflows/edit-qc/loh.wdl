version 1.1

import "../wdl-common/wdl/structs.wdl"

# Task to check for loss of heterozygosity (LOH) around a CRISPR cut site
task check_loh {
  meta {
    description: "Check for loss of heterozygosity around a CRISPR cut site by comparing parental HET SNPs to sample genotypes"
  }

  parameter_meta {
    sample_id: "Sample identifier"
    sample_bam: "Sample aligned BAM file"
    sample_bai: "Sample BAM index"
    parent_vcf: "Parent sample VCF (used to identify parental HET SNPs)"
    parent_vcf_index: "Parent VCF tabix index"
    ref_fasta: "Reference genome FASTA"
    ref_fasta_index: "Reference FASTA index (.fai)"
    expected_cut_site: "Expected cut site in chr:pos format (e.g. chr15:65869582)"
    cut_strand: "Strand of the guide at the cut site (+ or -)"
    window_size: "Window size in bp to check upstream and downstream of the cut site (default 100000)"
  }

  input {
    String sample_id
    File sample_bam
    File sample_bai
    File parent_vcf
    File parent_vcf_index
    File ref_fasta
    File ref_fasta_index
    String expected_cut_site
    String cut_strand
    Int window_size = 100000
    Int mem_gb = 8
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(sample_bam, "GB") + size(parent_vcf, "GB") + size(ref_fasta, "GB")) + 30

  command <<<
    set -euo pipefail

    python3 <<'EOF'
import subprocess
import sys

sample_bam    = "~{sample_bam}"
parent_vcf    = "~{parent_vcf}"
ref_fasta     = "~{ref_fasta}"
sample_id     = "~{sample_id}"
cut_site      = "~{expected_cut_site}"
cut_strand    = "~{cut_strand}"
window_size   = ~{window_size}

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


def get_parent_het_positions(region, label):
    """Return list of 0-based half-open BED lines for parental HET SNPs in region."""
    r_view = subprocess.run(
        ["bcftools", "view", "-g", "het", "-v", "snps", "-r", region, parent_vcf],
        capture_output=True, check=True,
    )
    r_query = subprocess.run(
        ["bcftools", "query", "-f", "%CHROM\t%POS0\t%END\n"],
        input=r_view.stdout, capture_output=True, check=True,
    )
    positions = [ln for ln in r_query.stdout.decode().splitlines() if ln.strip()]
    print(f"{label}: {len(positions)} parental HET SNP(s) in {region}", file=sys.stderr)
    return positions


def check_het_intact(positions, label):
    """
    For each parental HET position call genotypes in the sample BAM and return
    (pct_intact, n_het_in_sample, n_total_parental_het).
    """
    n_total = len(positions)
    if n_total == 0:
        return 0.0, 0, 0

    bed_file = f"{label}_het_positions.bed"
    with open(bed_file, "w") as fh:
        fh.write("\n".join(positions) + "\n")

    # mpileup over target positions, then call genotypes
    r_pileup = subprocess.run(
        [
            "bcftools", "mpileup",
            "-R", bed_file,
            "-f", ref_fasta,
            "-a", "AD",
            "--output-type", "u",
            sample_bam,
        ],
        capture_output=True, check=True,
    )
    r_call = subprocess.run(
        ["bcftools", "call", "-m", "--output-type", "v"],
        input=r_pileup.stdout, capture_output=True, check=True,
    )

    n_het = 0
    for line in r_call.stdout.decode().splitlines():
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        if len(fields) < 10:
            continue
        fmt = fields[8].split(":")
        smp = fields[9].split(":")
        if "GT" not in fmt:
            continue
        gt = smp[fmt.index("GT")].replace("|", "/")
        alleles = gt.split("/")
        if len(alleles) == 2 and "." not in alleles and alleles[0] != alleles[1]:
            n_het += 1

    pct = (n_het / n_total * 100) if n_total > 0 else 0.0
    print(
        f"{label}: {n_het}/{n_total} parental HET SNP(s) intact in sample = {pct:.2f}%",
        file=sys.stderr,
    )
    return pct, n_het, n_total


up_positions = get_parent_het_positions(upstream_region,   "upstream")
dn_positions = get_parent_het_positions(downstream_region, "downstream")

up_pct, up_het, up_total = check_het_intact(up_positions, "upstream")
dn_pct, dn_het, dn_total = check_het_intact(dn_positions, "downstream")

# Write float outputs (one value per file for WDL read_float())
with open(f"{sample_id}_het_intact_upstream_pct.txt", "w") as fh:
    fh.write(f"{up_pct:.4f}\n")
with open(f"{sample_id}_het_intact_downstream_pct.txt", "w") as fh:
    fh.write(f"{dn_pct:.4f}\n")

# Write human-readable summary
with open(f"{sample_id}_loh_summary.tsv", "w") as fh:
    fh.write("field\tvalue\n")
    fh.write(f"upstream_region\t{upstream_region}\n")
    fh.write(f"downstream_region\t{downstream_region}\n")
    fh.write(f"upstream_het_total\t{up_total}\n")
    fh.write(f"upstream_het_intact\t{up_het}\n")
    fh.write(f"het_intact_upstream_pct\t{up_pct:.4f}\n")
    fh.write(f"downstream_het_total\t{dn_total}\n")
    fh.write(f"downstream_het_intact\t{dn_het}\n")
    fh.write(f"het_intact_downstream_pct\t{dn_pct:.4f}\n")
EOF
  >>>

  output {
    Float het_intact_upstream_pct   = read_float("~{sample_id}_het_intact_upstream_pct.txt")
    Float het_intact_downstream_pct = read_float("~{sample_id}_het_intact_downstream_pct.txt")
    File  loh_summary               = "~{sample_id}_loh_summary.tsv"
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
