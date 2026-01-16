#!/usr/bin/env python3
"""
process_input_config.py: Process workflow input config and copy files to work directory

This script processes a workflow input JSON file by:
- Merging multiple HiFi BAM files per sample into single files
- Stripping unnecessary tags from BAM files (fi,ri,fp,rp,ip,pw,HP,PS,PC)
- Copying reference map files to work directory
- Copying expected edits JSON to work directory
- Generating a new inputs JSON with local file paths
"""

import json
import os
import shutil
import sys
import tempfile
import filecmp
import subprocess
from pathlib import Path


def should_copy_file(dest_path):
    """Check if file should be created (doesn't exist or has missing EOF block)"""
    if not dest_path.exists():
        return True

    # Run samtools quickcheck to verify BAM integrity
    result = subprocess.run(
        ["samtools", "quickcheck", str(dest_path)],
        capture_output=True,
        text=True
    )

    # Only regenerate if output mentions "missing EOF block"
    if result.returncode != 0 and "missing EOF block" in result.stderr:
        print(f"  File {dest_path} has missing EOF block, will regenerate")
        return True

    return False


def main():
    if len(sys.argv) != 2:
        print("Usage: python3 process_input_config.py <input_config_json>")
        sys.exit(1)

    # Read input config
    with open(sys.argv[1], 'r') as f:
        config = json.load(f)

    # Create new config with local paths
    new_config = {}
    inputs_dir = Path("inputs")

    # Process family samples
    if "humanwgs_family.family" in config:
        family = config["humanwgs_family.family"]
        new_family = {
            "family_id": family["family_id"],
            "samples": []
        }

        for sample in family["samples"]:
            new_sample = sample.copy()

            # Merge all HiFi reads files into one BAM while stripping tags
            merged_filename = f"{sample['sample_id']}_hifi_reads_merged.bam"
            dest_path = inputs_dir / merged_filename

            # Check if all input files exist
            all_exist = all(os.path.exists(hifi_read_path) for hifi_read_path in sample["hifi_reads"])

            if all_exist and should_copy_file(dest_path):
                with tempfile.TemporaryDirectory() as tmpdirname:
                    tmpdir_path = Path(tmpdirname)
                    temp_path = tmpdir_path / f"{sample['sample_id']}_hifi_reads_temp.bam"
                    print(f"Merging {len(sample['hifi_reads'])} HiFi read files -> {dest_path} (stripping tags: fi,ri,fp,rp,ip,pw,HP,PS,PC)")

                    # Step 1: Concatenate BAM files
                    print(f"  Step 1/2: Concatenating BAM files...")
                    subprocess.run([
                        "samtools", "cat",
                        "-o", str(temp_path)
                    ] + sample["hifi_reads"], check=True)

                    # Step 2: Strip tags using samtools reset
                    print(f"  Step 2/2: Stripping tags...")
                    subprocess.run([
                        "samtools", "reset",
                        "--thread", "15",
                        "--remove-tag", "fi,ri,fp,rp,ip,pw,HP,PS,PC",
                        "-o", str(dest_path),
                        str(temp_path)
                    ], check=True)

                    # Clean up temporary files
                    temp_path.unlink()

                new_sample["hifi_reads"] = [str(dest_path.absolute())]
            elif all_exist:
                print(f"Skipping merge for {sample['sample_id']} -> {dest_path} (already exists)")
                new_sample["hifi_reads"] = [str(dest_path.absolute())]
            else:
                # Some files don't exist, keep original paths and warn
                print(f"Warning: Some HiFi reads files not found for {sample['sample_id']}")
                for hifi_read_path in sample["hifi_reads"]:
                    if not os.path.exists(hifi_read_path):
                        print(f"  Missing: {hifi_read_path}")
                new_sample["hifi_reads"] = sample["hifi_reads"]  # Keep original paths

            new_family["samples"].append(new_sample)

        new_config["humanwgs_family.family"] = new_family

    # Copy reference files if they exist locally
    ref_files = ["ref_map_file", "somatic_map_file", "tertiary_map_file"]
    for ref_key in ref_files:
        full_key = f"humanwgs_family.{ref_key}"
        if full_key in config:
            ref_path = config[full_key]
            if os.path.exists(ref_path):
                filename = f"{ref_key}.tsv"
                dest_path = inputs_dir / filename

                if should_copy_file(dest_path):
                    print(f"Copying {ref_path} -> {dest_path}")
                    shutil.copy2(ref_path, dest_path)
                else:
                    print(f"Skipping {ref_path} -> {dest_path} (already exists)")

                new_config[full_key] = str(dest_path.absolute())
            else:
                print(f"Warning: Reference file not found: {ref_path}")
                new_config[full_key] = ref_path  # Keep original path

    # Copy expected_edits file if it exists
    if "humanwgs_family.expected_edits" in config:
        expected_edits_path = config["humanwgs_family.expected_edits"]
        if os.path.exists(expected_edits_path):
            filename = "expected_edits.json"
            dest_path = inputs_dir / filename

            if dest_path.exists():
                print(f"Skipping {expected_edits_path} -> {dest_path} (already exists)")
            else:
                print(f"Copying {expected_edits_path} -> {dest_path}")
                shutil.copy2(expected_edits_path, dest_path)

            new_config["humanwgs_family.expected_edits"] = str(dest_path.absolute())
        else:
            print(f"Warning: Expected edits file not found: {expected_edits_path}")
            new_config["humanwgs_family.expected_edits"] = expected_edits_path

    # Copy other config values
    for key, value in config.items():
        if key not in new_config:
            new_config[key] = value

    # Write new inputs JSON
    with open("family.hpc.inputs.json", "w") as f:
        json.dump(new_config, f, indent=2)

    print("Generated family.hpc.inputs.json with local file paths")


if __name__ == "__main__":
    main()
