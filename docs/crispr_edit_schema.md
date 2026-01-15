# CRISPR Edit JSON Schema

Expected edit files describe anticipated genomic changes from CRISPR experiments. These files are provided per edited sample via the `expected_edits` field in the sample struct.

## Schema Overview

```json
{
  "ref": "hg38",
  "edits": [
    {
      "id": "edit_identifier",
      "target": {
        "chr": "chr12",
        "start": 92770637,
        "end": 92929295,
        "symbol": "EEA1"
      },
      "edit": {
        "type": "INS",
        "strand": "+",
        "guide": "GTGGTGGTTAAACCATGTTA",
        "left_ha": "CAGGCTGCCTAGTC...",
        "payload": "TCTGGTGGTACTTC...",
        "right_ha": "TTAAGGAGGATTTA..."
      }
    }
  ]
}
```

## Top-level Fields

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `ref` | String | Yes | Reference assembly. Currently only `"hg38"` is supported. |
| `edits` | Array | Yes | Array of edit objects. |

## Edit Object

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `id` | String | Yes | Unique identifier for this edit. Used in output filenames and reports. |
| `target` | Object | Yes | Genomic region containing the edit site. |
| `edit` | Object | Yes | Edit-specific information. |

## Target Object

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `chr` | String | Yes | Chromosome name (e.g., `"chr1"`, `"chrX"`). |
| `start` | Integer | Yes | Start position (1-based, inclusive). |
| `end` | Integer | Yes | End position (1-based, inclusive). |
| `symbol` | String | No | Gene symbol overlapping the edit site. |

## Edit Details Object

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `type` | String | Yes | Edit type: `"INS"`, `"DEL"`, or `"SNV"`. |
| `strand` | String | Yes | Strand orientation: `"+"` or `"-"`. |
| `guide` | String | Yes | Guide RNA sequence (20bp). |
| `left_ha` | String | Yes | Left homology arm sequence. |
| `payload` | String | Yes | Payload sequence (inserted sequence for INS, empty for DEL). |
| `right_ha` | String | Yes | Right homology arm sequence. |
| `payload_components` | Array | No | Optional named subsequences within payload for annotation. |

**Note**: The `left_ha + payload + right_ha` represents the HDR template sequence.

## Edit Types

### INS (Insertion)

Large insertions, typically knock-ins of tags or functional domains.

```json
{
  "type": "INS",
  "strand": "+",
  "guide": "GTGGTGGTTAAACCATGTTA",
  "left_ha": "CAGGCTGCCTAGTCGTTCGC...",
  "payload": "TCTGGTGGTACTTCCGGGTA...",
  "right_ha": "TTAAGGAGGATTTACAGAGG..."
}
```

### DEL (Deletion)

Large deletions, typically exon knockouts.

```json
{
  "type": "DEL",
  "strand": "+",
  "guide": "AGCTGATCGATCGATCGATC",
  "left_ha": "AGCTGATCGATCGATCGATC...",
  "payload": "",
  "right_ha": "TCGATCGATCGATCGATCGA..."
}
```

### SNV (Single Nucleotide Variant)

Point mutations or small indels.

```json
{
  "type": "SNV",
  "strand": "+",
  "guide": "CGGCACATGCAGGCCGAGCT",
  "left_ha": "CTGCAGAAGCGCCTGGCAGT...",
  "payload": "C",
  "right_ha": "GACCTGGGCGACTCCGAGCC..."
}
```

## Payload Components (Optional)

For complex insertions with multiple domains, specify named subsequences for detailed annotation:

```json
{
  "type": "INS",
  "payload": "TCTGGTGGTACTTCCGGGTACC...GGAAGCCCTACTTCCACCGAAGAAGGCACGTCAACCGAACCAAGTG",
  "payload_components": [
    {
      "name": "3xHA",
      "seq": "TCTGGTGGTACTTCCGGGTACCCATACGATGTTCCAGATTACGCTTATCCGTATGACGTCCCGGACTATGCATACCCCACGATGTACCCGATTACGCC"
    },
    {
      "name": "mEGFP",
      "seq": "GCTGGAGCTGGCGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCGATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCAAGCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAG"
    },
    {
      "name": "XTEN80",
      "seq": "GGTGGCGGATTGGAAGTTTTGTTTCAAGGTCCAGGATCTAACGGCAGCGGAGGGCCGAGCTCTGGCGCACCCCCACCAAGTGGAGGGTCACCTGCCGGGTCCCCAACATCTACTGAAGAAGGCACCAGCGAATCCGCAACGCCCGAGTCAGGCCCTGGTACCTCCACAGAACCATCTGAAGGT"
    },
    {
      "name": "linker",
      "seq": "AGTGCGCCTGGTTCCCCAGCTGGAAGCCCTACTTCCACCGAAGAAGGCACGTCAACCGAACCAAGTGAAGGATCTGCCCCTGGGACCAGCACTGAACCATCTGAG"
    }
  ]
}
```

The pipeline will align each component to the sample genome for detailed reporting of which parts of the payload are present.

## Generating Edit Files from GenBank

For complex knock-ins, use `genbank_to_crispr_json.py` to convert annotated HDR template plasmids:

```bash
python3 genbank_to_crispr_json.py \
  -i plasmid.gb \
  -o edit.json \
  -e edit_name \
  -t ENST00000123456 \
  -g GGGCCACATCAGCGCGATCC \
  --type INS
```

### GenBank File Requirements

The input GenBank file must contain annotated features with specific labels:

**Required features**:
- **Left homology arm**: Label matching `5' Homology`, `5'HA`, `Left HA`, or `Left Homology Arm` (case-insensitive)
- **Right homology arm**: Label matching `3' Homology`, `3'HA`, `Right HA`, or `Right Homology Arm` (case-insensitive)
- **Payload**: Label exactly matching `Payload` (case-insensitive)

**Optional features**:
- **Payload components**: Any `misc_feature` annotations completely contained within the payload region will be extracted as named components (e.g., linker, fluorescent protein, epitope tag). The script ignores CDS features to focus on functional domains.

**Notes**:
- The script automatically fetches target gene coordinates from Ensembl using the provided transcript ID
- If strand is not specified, it uses the gene's strand from Ensembl
- Payload components are sorted by position and included if found

## Validation

The pipeline validates edit files against `workflows/edit-qc/crispr_edit.schema.json`. Common errors:

- Missing required fields (`type`, `strand`, `guide`, homology arms)
- Invalid edit type (must be INS, DEL, or SNV)
- Invalid strand (must be + or -)
- Invalid reference (must be hg38)
- Malformed coordinates (start >= end)

## See Also

- [example_expected_edit.json](../example_expected_edit.json) - Complete example
- [family.md](family.md#sample-struct) - Sample struct documentation
- [biohub-setup.md](biohub-setup.md) - Setup and edit file creation
