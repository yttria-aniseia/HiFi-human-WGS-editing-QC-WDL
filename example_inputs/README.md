# Example Input Configurations

Example input JSON files for common CRISPR editing QC workflows.

## Files

### single_edited_sample.json

Single edited sample without parent. Use when:
- Only have sequencing data for the edited clone
- Parent sample is unavailable or not needed
- Quick validation of a single edit

**Note**: Without a parent sample, the pipeline cannot perform parent-based variant filtering or inheritance analysis.

### parent_and_edited_clones.json

Parent sample with multiple edited clones. Use when:
- Multiple clones from the same CRISPR experiment
- Different edits in the same parent line (clone_A* vs clone_B*)
- Need parent-based filtering to identify edit-specific variants

**Best practice**: Always include the parental cell line when available for better variant filtering.

### Example edit file structure

See [example_expected_edit.json](../example_expected_edit.json) for a complete edit file example.

For multiple edits in one sample, the edit JSON can contain multiple entries:

```json
{
  "ref": "hg38",
  "edits": [
    {
      "id": "gene_A_n_terminal_tag",
      "target": {...},
      "edit": {...}
    },
    {
      "id": "gene_B_c_terminal_tag",
      "target": {...},
      "edit": {...}
    }
  ]
}
```

### Real publication examples (paths redacted)

The following subfolders hold **complete, real input files as run for publication**. All
absolute filesystem paths have been redacted to `/path/to/...` placeholders — substitute
your own paths before use. Agents should consult these for up-to-date, full-fidelity
examples of every input format:

- `family_examples/` — full `humanwgs_family` input configs (parent + edited clones,
  precomputed assemblies, ref/somatic/tertiary map references).
- `edit_json/` — expected-edit JSONs (the `-extracontext` variants include
  `donor_context_before`/`after` and `payload_components`).
- `genbank/` — source Benchling GenBank plasmid maps that the edit JSONs were derived
  from (via `genbank_to_crispr_json.py`).

## Usage

1. Copy the appropriate example file
2. Update file paths to your data:
   - HiFi read BAM files (`hifi_reads`)
   - Expected edit JSON files (`expected_edits`)
3. Update sample metadata (IDs, sex, relationships)
4. Run the workflow:

```bash
./scripts/launch.sh your_config.json --work-dir analysis_name
```

## See Also

- [example_input_config.json](../example_input_config.json) - Main example with comments
- [docs/crispr_edit_schema.md](../docs/crispr_edit_schema.md) - Edit JSON format documentation
- [docs/biohub-setup.md](../docs/biohub-setup.md) - Setup instructions
