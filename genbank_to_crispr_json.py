#!/usr/bin/env python3
"""
Convert GenBank format CRISPR edit plasmid files to CRISPR edit JSON format.

Extracts homology arms, payload sequence, and payload components from GenBank
annotation and fetches target gene location from Ensembl.
"""

import argparse
import json
import re
import sys
from typing import Dict, List, Optional, Tuple

try:
    from Bio import SeqIO
    from Bio.SeqFeature import SeqFeature
except ImportError as e:
    print(f"Error: Biopython not found. Please install it:", file=sys.stderr)
    print(f"  pip install biopython", file=sys.stderr)
    sys.exit(1)

import urllib.request
import urllib.error


def find_feature_by_label(features: List[SeqFeature], patterns: List[str],
                          feature_name: str) -> Optional[SeqFeature]:
    """
    Find a feature by matching its label against a list of patterns (case-insensitive).

    Args:
        features: List of SeqFeature objects from the GenBank record
        patterns: List of regex patterns to match against the label
        feature_name: Name of the feature (for error messages)

    Returns:
        The matching SeqFeature or None if not found
    """
    for feature in features:
        if 'label' in feature.qualifiers:
            label = feature.qualifiers['label'][0]
            for pattern in patterns:
                if re.search(pattern, label, re.IGNORECASE):
                    return feature
    return None


def get_feature_sequence(feature: SeqFeature, record) -> str:
    """Extract the sequence for a feature from the record."""
    return str(feature.extract(record.seq))


def find_payload_components(features: List[SeqFeature], record, payload_start: int,
                            payload_end: int) -> List[Dict[str, str]]:
    """
    Find all sub-features within the payload region.

    Args:
        features: List of all SeqFeature objects
        record: The SeqRecord object containing the sequence
        payload_start: Start position of payload (0-based)
        payload_end: End position of payload (0-based)

    Returns:
        List of NamedSequence dicts with 'name' and 'seq' keys
    """
    components = []

    for feature in features:
        # Skip the payload feature itself
        if 'label' in feature.qualifiers:
            label = feature.qualifiers['label'][0]
            if re.search(r'^payload$', label, re.IGNORECASE):
                continue

        # Check if feature is within payload region
        feat_start = int(feature.location.start)
        feat_end = int(feature.location.end)

        # Feature must be completely contained within payload
        if feat_start >= payload_start and feat_end <= payload_end:
            if 'label' in feature.qualifiers:
                name = feature.qualifiers['label'][0]
                # Only process misc_feature types, skip CDS
                if feature.type == 'misc_feature':
                    # Get the sequence
                    seq = str(feature.extract(record.seq))

                    components.append({
                        'name': name,
                        'seq': seq
                    })

    # Sort by start position to maintain order
    components.sort(key=lambda x: int(features[
        [i for i, f in enumerate(features)
         if f.qualifiers.get('label', [''])[0] == x['name']][0]
    ].location.start))

    return components


def fetch_gene_location(transcript_id: str) -> Dict[str, any]:
    """
    Fetch gene location from Ensembl REST API for a given transcript ID.

    Args:
        transcript_id: ENSEMBL transcript ID (e.g., ENST00000356245)

    Returns:
        Dict with keys: chr, start, end, symbol, strand

    Raises:
        Exception if the transcript cannot be found or API fails
    """
    # Use Ensembl REST API
    server = "https://rest.ensembl.org"
    ext = f"/lookup/id/{transcript_id}?expand=1"

    try:
        req = urllib.request.Request(
            server + ext,
            headers={"Content-Type": "application/json"}
        )
        with urllib.request.urlopen(req) as response:
            if response.status != 200:
                raise Exception(f"Failed to fetch transcript {transcript_id}: HTTP {response.status}")
            data = json.loads(response.read().decode())

        # Extract chromosome (remove 'chr' prefix if present, add it back)
        chrom = data['seq_region_name']
        if not chrom.startswith('chr'):
            chrom = f'chr{chrom}'

        # Get gene symbol from display_name or gene name
        symbol = None
        if 'display_name' in data:
            symbol = data['display_name']
        elif 'Parent' in data:
            # Fetch parent gene
            gene_id = data['Parent']
            gene_ext = f"/lookup/id/{gene_id}"
            gene_req = urllib.request.Request(
                server + gene_ext,
                headers={"Content-Type": "application/json"}
            )
            try:
                with urllib.request.urlopen(gene_req) as gene_response:
                    if gene_response.status == 200:
                        gene_data = json.loads(gene_response.read().decode())
                        symbol = gene_data.get('display_name', gene_id)
            except urllib.error.URLError:
                pass  # Continue without symbol

        return {
            'chr': chrom,
            'start': data['start'],
            'end': data['end'],
            'symbol': symbol,
            'strand': '+' if data['strand'] == 1 else '-'
        }

    except urllib.error.URLError as e:
        raise Exception(f"Network error fetching transcript {transcript_id}: {e}")
    except (KeyError, json.JSONDecodeError) as e:
        raise Exception(f"Error parsing Ensembl response for {transcript_id}: {e}")


def parse_genbank_to_crispr_json(genbank_path: str, edit_id: str,
                                  transcript_id: str, grna_seq: str,
                                  edit_type: str = "INS",
                                  strand: Optional[str] = None) -> Dict:
    """
    Parse GenBank file and generate CRISPR edit JSON.

    Args:
        genbank_path: Path to GenBank format file
        edit_id: Identifier for this edit
        transcript_id: ENSEMBL transcript ID (e.g., ENST00000356245)
        grna_seq: gRNA sequence
        edit_type: Type of edit (SNV, INS, or DEL), default INS
        strand: Strand of edit (+ or -), if None will use gene strand

    Returns:
        Dict matching the CRISPR edit JSON schema
    """
    # Parse GenBank file
    try:
        record = SeqIO.read(genbank_path, "genbank")
    except Exception as e:
        print(f"Error reading GenBank file {genbank_path}: {e}", file=sys.stderr)
        sys.exit(1)

    features = list(record.features)

    # Find left homology arm
    left_ha_patterns = [r"5['\u2019\u0027]?[\s_]*(homology|ha)", r"left\s+ha\b", r"left\s+homology\s+arm"]
    left_ha_feature = find_feature_by_label(features, left_ha_patterns, "left homology arm")
    if not left_ha_feature:
        print("Error: Could not find left homology arm feature", file=sys.stderr)
        print("  Expected label matching: 5' Homology, Left HA, or Left Homology Arm", file=sys.stderr)
        sys.exit(1)

    # Find right homology arm
    right_ha_patterns = [r"3['\u2019\u0027]?[\s_]*(homology|ha)", r"right\s+ha\b", r"right\s+homology\s+arm"]
    right_ha_feature = find_feature_by_label(features, right_ha_patterns, "right homology arm")
    if not right_ha_feature:
        print("Error: Could not find right homology arm feature", file=sys.stderr)
        print("  Expected label matching: 3' Homology, Right HA, or Right Homology Arm", file=sys.stderr)
        sys.exit(1)

    # Find payload
    payload_patterns = [r"^payload$"]
    payload_feature = find_feature_by_label(features, payload_patterns, "payload")
    if not payload_feature:
        print("Error: Could not find payload feature", file=sys.stderr)
        print("  Expected label: Payload", file=sys.stderr)
        sys.exit(1)

    # Extract sequences
    left_ha_seq = get_feature_sequence(left_ha_feature, record)
    right_ha_seq = get_feature_sequence(right_ha_feature, record)
    payload_seq = get_feature_sequence(payload_feature, record)

    # Find payload components
    payload_start = int(payload_feature.location.start)
    payload_end = int(payload_feature.location.end)
    payload_components = find_payload_components(features, record, payload_start, payload_end)

    # Fetch gene location
    print(f"Fetching gene location for {transcript_id}...", file=sys.stderr)
    try:
        gene_info = fetch_gene_location(transcript_id)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    # Use gene strand if not specified
    if strand is None:
        strand = gene_info['strand']

    # Build the JSON structure
    crispr_json = {
        "ref": "hg38",
        "edits": [
            {
                "id": edit_id,
                "target": {
                    "chr": gene_info['chr'],
                    "start": gene_info['start'],
                    "end": gene_info['end'],
                },
                "edit": {
                    "type": edit_type,
                    "strand": strand,
                    "guide": grna_seq,
                    "left_ha": left_ha_seq,
                    "payload": payload_seq,
                    "right_ha": right_ha_seq,
                }
            }
        ]
    }

    # Add optional symbol if available
    if gene_info.get('symbol'):
        crispr_json['edits'][0]['target']['symbol'] = gene_info['symbol']

    # Add payload components if found
    if payload_components:
        crispr_json['edits'][0]['edit']['payload_components'] = payload_components

    return crispr_json


def main():
    parser = argparse.ArgumentParser(
        description='Convert GenBank CRISPR plasmid to CRISPR edit JSON format',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  %(prog)s -i plasmid.gb -e g3bp1_n_term_tag -t ENST00000356245 \\
           -g CGCCCGACCAGCAGGGGACT -o g3bp1_edit.json

  %(prog)s -i plasmid.gb -e my_edit -t ENST00000356245 \\
           -g CGCCCGACCAGCAGGGGACT --type INS --strand +
        """
    )

    parser.add_argument('-i', '--input', required=True,
                       help='Input GenBank format file')
    parser.add_argument('-o', '--output', required=True,
                       help='Output JSON file')
    parser.add_argument('-e', '--edit-id', required=True,
                       help='Edit identifier (e.g., g3bp1_n_term_tag)')
    parser.add_argument('-t', '--transcript', required=True,
                       help='ENSEMBL transcript ID (e.g., ENST00000356245)')
    parser.add_argument('-g', '--grna', required=True,
                       help='gRNA sequence (e.g., CGCCCGACCAGCAGGGGACT)')
    parser.add_argument('--type', choices=['SNV', 'INS', 'DEL'], default='INS',
                       help='Edit type (default: INS)')
    parser.add_argument('--strand', choices=['+', '-'],
                       help='Edit strand (default: use gene strand)')

    args = parser.parse_args()

    # Generate JSON
    crispr_json = parse_genbank_to_crispr_json(
        args.input,
        args.edit_id,
        args.transcript,
        args.grna,
        args.type,
        args.strand
    )

    # Write output
    with open(args.output, 'w') as f:
        json.dump(crispr_json, f, indent=2)

    print(f"Successfully wrote CRISPR edit JSON to {args.output}", file=sys.stderr)
    print(f"  Edit ID: {args.edit_id}", file=sys.stderr)
    print(f"  Target: {crispr_json['edits'][0]['target']['chr']}:"
          f"{crispr_json['edits'][0]['target']['start']}-"
          f"{crispr_json['edits'][0]['target']['end']}", file=sys.stderr)
    if 'symbol' in crispr_json['edits'][0]['target']:
        print(f"  Gene: {crispr_json['edits'][0]['target']['symbol']}", file=sys.stderr)
    print(f"  Payload components: {len(crispr_json['edits'][0]['edit'].get('payload_components', []))}",
          file=sys.stderr)


if __name__ == '__main__':
    main()
