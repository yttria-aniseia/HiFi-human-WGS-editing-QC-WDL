version 1.1

import "../wdl-common/wdl/structs.wdl"

task annotate_consensus_genbank {
  meta {
    description: "Create GenBank files annotating edit parts on haplotype consensus sequences"
  }

  parameter_meta {
    parts_alignment_tsv: "Blastn tabular output from blastn_remap_parts (edit parts aligned to consensus sequences)"
    consensus_fasta: "FASTA file with haplotype consensus sequences"
    crispr_edit_json: "JSON file describing CRISPR edit structure"
    coverage_threshold: "Minimum coverage threshold for faithful alignments"
  }

  input {
    File parts_alignment_tsv
    File consensus_fasta
    File crispr_edit_json
    Float coverage_threshold = 0.8
    String sample_id
    RuntimeAttributes runtime_attributes
  }

  String out_prefix = "~{sample_id}_edit_annotations"

  command <<<
    set -euxo pipefail

    cat > create_genbank_annotations.py <<'EOF'
#!/usr/bin/env python3
"""
create_genbank_annotations.py - Annotate edit parts on haplotype consensus sequences

Creates GenBank files for each haplotype consensus sequence with proper source and gene
features indicating genomic origin, plus edit part annotations.
"""

import argparse
import json
import re
from collections import defaultdict
from datetime import datetime
from typing import Dict, List, Tuple
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

# Type aliases
Hit = Tuple[int, int, int, int, str, int]  # q_start, q_end, t_start, t_end, strand, n_match

def parse_blastn_line(line: str) -> Dict:
    """Parse a single blastn tabular line (outfmt 6 qseqid qlen qstart qend sseqid slen sstart send nident mismatch gapopen gaps length sstrand bitscore)"""
    f = line.rstrip("\n").split("\t")
    if len(f) < 15:
        raise ValueError(f"Expected >=15 columns, got {len(f)}")
    sstart, send = int(f[6]), int(f[7])
    return dict(
        query_name=f[0], query_len=int(f[1]),
        q_start=int(f[2]) - 1, q_end=int(f[3]),
        strand="+" if f[13] == "plus" else "-",
        target_name=f[4], target_len=int(f[5]),
        t_start=min(sstart, send) - 1, t_end=max(sstart, send),
        n_match=int(f[8]), n_align=int(f[12]),
    )

def load_faithful_alignments(tsv_path: str, cov_thr: float) -> Tuple[Dict[str, Dict[str, List[Hit]]], Dict[str, int]]:
    """
    Load only faithful alignments (high coverage, low gaps) from blastn tabular output.
    Returns: ({target_consensus_name: {query_part_name: [Hit, ...]}}, {query_name: length})
    """
    faithful = defaultdict(lambda: defaultdict(list))
    qlen = {}

    with open(tsv_path) as fh:
        for ln in fh:
            if not ln.strip() or ln.startswith('#'):
                continue
            try:
                r = parse_blastn_line(ln)
            except ValueError:
                continue
            qlen.setdefault(r["query_name"], r["query_len"])

            # Calculate coverage and gap metrics
            cov_hit = r["n_match"] / r["query_len"]
            gap_ok = r["n_align"] - r["n_match"] <= 0.05 * r["n_match"]

            # Only keep faithful alignments
            if cov_hit >= cov_thr and gap_ok:
                hit: Hit = (
                    r["q_start"], r["q_end"],
                    r["t_start"], r["t_end"],
                    r["strand"], r["n_match"]
                )
                faithful[r["target_name"]][r["query_name"]].append(hit)

    return faithful, qlen

def load_edit_metadata(json_path: str) -> Dict:
    """Load edit metadata from CRISPR edit JSON"""
    with open(json_path) as f:
        data = json.load(f)

    # Get first edit's metadata
    edit = data['edits'][0]
    return {
        'edit_id': edit['id'],
        'chromosome': edit['target']['chr'],
        'start': edit['target']['start'],
        'end': edit['target']['end'],
        'symbol': edit['target']['symbol']
    }

def parse_group_label(seq_id: str, description: str) -> str:
    """
    Extract a human-readable group label from MSA consensus FASTA headers.
    Header format: '>group_{N} n_reads={K}' or '>group_{N}_singleton n_reads=1'
    or '>group_{N}_hap_{L} n_reads={K}' for split haplotypes.
    Falls back gracefully for anything else.
    """
    import re
    m = re.match(r'(group_\d+(?:_(?:singleton|hap_\w+))?)', seq_id)
    if m:
        nr = ''
        for token in description.split():
            if token.startswith('n_reads='):
                nr = token[8:]
        parts = [m.group(1)]
        if nr:
            parts.append(f"n_reads={nr}")
        return ' '.join(parts)
    # Fallback
    return seq_id

def create_source_and_gene_features(seq_len: int, metadata: Dict, group_label: str) -> List[SeqFeature]:
    """
    Create source and gene features for the genomic context of the consensus sequence.
    """
    features = []

    source_feature = SeqFeature(
        location=FeatureLocation(0, seq_len, strand=+1),
        type="source",
        qualifiers={
            "organism": ["Homo sapiens"],
            "mol_type": ["genomic DNA"],
            "db_xref": ["taxon:9606"],
            "chromosome": [metadata['chromosome'].replace('chr', '')],
            "map": [f"{metadata['start']}-{metadata['end']}"],
            "note": [f"Consensus from ~{sample_id} {group_label}"]
        }
    )
    features.append(source_feature)

    gene_feature = SeqFeature(
        location=FeatureLocation(0, seq_len, strand=+1),
        type="gene",
        qualifiers={
            "label": [metadata['symbol']],
            "gene": [metadata['symbol']],
            "note": [f"color: #a6acb3"]
        }
    )
    features.append(gene_feature)

    return features

def create_edit_part_features(alignments: Dict[str, List[Hit]], query_lens: Dict[str, int]) -> List[SeqFeature]:
    """
    Create SeqFeature objects for each aligned edit part.

    Args:
        alignments: {query_part_name: [Hit, ...]}
        query_lens: {query_part_name: length}

    Returns:
        List of SeqFeature objects for edit parts
    """
    features = []

    for query_name, hits in alignments.items():
        # Parse part type from query name (e.g., "SYMBOL_exonX_left_HA")
        color = "#ff0000"
        if '_left_HA' in query_name:
            query_name = "5' Homology Arm"
            part_type = 'left_HA'
            color = "#00ffff"
        elif '_right_HA' in query_name:
            query_name = "3' Homology Arm"
            part_type = 'right_HA'
            color = "#00ffff"
        elif '_payload' in query_name:
            query_name = 'Payload'
            part_type = 'payload'
            color = "#75c6a9"
        elif '_HDR_template' in query_name:
            query_name = query_name.replace("_", " ")
            part_type = 'HDR_template'
            color = "#c0c0c0"
        elif '_WT_template' in query_name:
            query_name = query_name.replace("_", " ")
            part_type = 'WT_template'
            color = "#808080"
        else:
            query_name = query_name.replace("_", " ")
            part_type = 'unknown'

        for hit in hits:
            q_start, q_end, t_start, t_end, strand, n_match = hit

            # Calculate coverage for this hit
            coverage = n_match / query_lens[query_name] if query_name in query_lens else 0.0

            # Create feature location
            if strand == '+':
                location = FeatureLocation(t_start, t_end, strand=+1)
            else:
                location = FeatureLocation(t_start, t_end, strand=-1)

            # Create SeqFeature
            feature = SeqFeature(
                location=location,
                type="misc_feature",
                qualifiers={
                    "label": [query_name],
                    "note": [f"color: {color}"],
                    "part_type": [part_type],
                    "coverage": [f"{coverage:.3f}"],
                    "alignment_score": [str(n_match)]
                }
            )

            features.append(feature)

    return features

def annotate_consensus_sequences(faithful_alignments: Dict[str, Dict[str, List[Hit]]],
                                 fasta_path: str,
                                 query_lens: Dict[str, int],
                                 edit_metadata: Dict,
                                 output_prefix: str) -> Tuple[int, int]:
    """
    Annotate haplotype consensus sequences with source, gene, and edit part features.

    Returns: (num_sequences_annotated, num_edit_features_added)
    """
    # Load consensus sequences
    consensus_seqs = {rec.id: rec for rec in SeqIO.parse(fasta_path, "fasta")}

    annotated_records = []
    total_edit_features = 0

    for seq_id, alignments in faithful_alignments.items():
        if seq_id not in consensus_seqs:
            print(f"Warning: Sequence {seq_id} has alignments but no sequence in FASTA")
            continue

        # Get consensus sequence
        original_rec = consensus_seqs[seq_id]

        group_label = parse_group_label(seq_id, original_rec.description)

        # Create source and gene features
        source_gene_features = create_source_and_gene_features(
            len(original_rec.seq), edit_metadata, group_label
        )

        # Create edit part features
        edit_features = create_edit_part_features(alignments, query_lens)

        # Combine all features
        all_features = source_gene_features + edit_features

        # Create annotated SeqRecord
        annotated_rec = SeqRecord(
            seq=original_rec.seq,
            id=original_rec.id,
            name=original_rec.id[:16],  # GenBank name max 16 chars
            description=f"~{sample_id} {edit_metadata['symbol']} {group_label} consensus {edit_metadata['chromosome']}:{edit_metadata['start']}-{edit_metadata['end']}",
            features=all_features,
            annotations={
                "source": "~{sample_id}",
                "organism": "Homo sapiens",
                "molecule_type": "DNA",
                "topology": "linear",
                "date": datetime.now().strftime("%d-%b-%Y").upper()
            }
        )

        annotated_records.append(annotated_rec)
        total_edit_features += len(edit_features)

    # Write GenBank file
    output_gbk = f"{output_prefix}.gbk"
    SeqIO.write(annotated_records, output_gbk, "genbank")

    return len(annotated_records), total_edit_features

def main():
    parser = argparse.ArgumentParser(
        description="Create GenBank annotations for CRISPR edit parts on haplotype consensus sequences"
    )
    parser.add_argument("tsv", help="Blastn tabular alignment file")
    parser.add_argument("fasta", help="FASTA file with consensus sequences")
    parser.add_argument("json", help="CRISPR edit JSON file")
    parser.add_argument("-c", "--coverage_threshold", type=float, default=0.8,
                       help="Minimum coverage threshold for faithful alignments")
    parser.add_argument("-o", "--output_prefix", default="edit_annotations",
                       help="Output file prefix")

    args = parser.parse_args()

    print("Loading faithful alignments from blastn TSV...")
    faithful_alignments, query_lens = load_faithful_alignments(args.tsv, args.coverage_threshold)
    print(f"Found {len(faithful_alignments)} consensus sequences with faithful alignments")

    print("Loading edit metadata...")
    edit_metadata = load_edit_metadata(args.json)
    print(f"Edit: {edit_metadata['edit_id']} at {edit_metadata['chromosome']}:{edit_metadata['start']}-{edit_metadata['end']}")

    print("Annotating consensus sequences and creating GenBank files...")
    num_seqs, num_features = annotate_consensus_sequences(
        faithful_alignments,
        args.fasta,
        query_lens,
        edit_metadata,
        args.output_prefix
    )

    # Write summary
    summary_file = f"{args.output_prefix}_summary.txt"
    with open(summary_file, 'w') as f:
        f.write(f"CRISPR Edit Consensus Annotation Summary\n")
        f.write(f"=========================================\n")
        f.write(f"Edit: {edit_metadata['edit_id']}\n")
        f.write(f"Gene: {edit_metadata['symbol']}\n")
        f.write(f"Location: {edit_metadata['chromosome']}:{edit_metadata['start']}-{edit_metadata['end']}\n")
        f.write(f"Consensus sequences annotated: {num_seqs}\n")
        f.write(f"Edit part features added: {num_features}\n")
        if num_seqs > 0:
            f.write(f"Average features per sequence: {num_features/num_seqs:.2f}\n")
        f.write(f"Coverage threshold: {args.coverage_threshold}\n")

    print(f"✓ Created GenBank file with {num_seqs} annotated consensus sequences ({num_features} edit part features)")
    print(f"✓ Summary written to {summary_file}")

if __name__ == "__main__":
    main()
EOF

    chmod +x create_genbank_annotations.py

    # Run the annotation script
    python3 create_genbank_annotations.py \
      ~{parts_alignment_tsv} \
      ~{consensus_fasta} \
      ~{crispr_edit_json} \
      -c ~{coverage_threshold} \
      -o ~{out_prefix}
  >>>

  output {
    File genbank_annotations = "~{out_prefix}.gbk"
    File annotation_summary = "~{out_prefix}_summary.txt"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87"
    cpu: 2
    memory: "8 GB"
    disk: "20 GB"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
  }
}
