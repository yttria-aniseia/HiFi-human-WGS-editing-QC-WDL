version 1.1

import "../wdl-common/wdl/structs.wdl"

task knock_knock {
  meta {
    description: "Run knock-knock to classify editing outcomes at an HDR target site and generate visualizations"
  }

  parameter_meta {
    sample_id: "Sample identifier; used as knock-knock sample_name in the sample sheet"
    crispr_edit_json: "JSON describing the expected CRISPR edit (target coordinates, guide sequence, HA/payload sequences)"
    region_reads_fastq: "FASTQ of reads overlapping the edit region (optionally gzip-compressed). Should contain primary alignments only (samtools fastq default behavior satisfies this)."
    ref_fasta: "Reference genome FASTA used to extract local genomic context for knock-knock build-strategies"
    ref_fasta_index: "Reference genome FASTA index (.fai)"
    genome_name: "Genome label used in knock-knock strategy naming (e.g. 'hg38'). Must match the genome field in the sample sheet."
    flanking_bp: "Flanking bases extracted around the target region to build the local reference context. Must exceed 6000 (the buffer knock-knock uses when extracting the target window from the protospacer alignment)."
    threads: "Parallel worker processes passed to knock-knock parallel"
    mem_gb: "Memory in GiB"
  }

  input {
    String sample_id
    File   crispr_edit_json
    File   region_reads_fastq
    File   ref_fasta
    File   ref_fasta_index
    String genome_name = "hg38"
    Int    flanking_bp = 10000
    Int    threads     = 4
    Int    mem_gb      = 16
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(region_reads_fastq, "GB") * 3) + 30

  command <<<
    set -euxo pipefail

    # ── 1. Parse crispr_edit_json ──────────────────────────────────────────────
    python3 <<'PYEOF'
import json, sys

with open("~{crispr_edit_json}") as f:
    data = json.load(f)

edit      = data['edits'][0]
target    = edit['target']
edit_info = edit['edit']

symbol    = target['symbol']
chrom     = target['chr']
t_start   = target['start']
t_end     = target['end']
guide_seq = edit_info['guide']
left_ha   = edit_info.get('left_ha',  '')
payload   = edit_info.get('payload',  '')
right_ha  = edit_info.get('right_ha', '')
ctx_before = edit_info.get('donor_context_before', '')
ctx_after  = edit_info.get('donor_context_after',  '')
donor_seq = ctx_before + left_ha + payload + right_ha + ctx_after
donor_name = f"{symbol}_donor"

for fname, val in [
    ('symbol.txt',     symbol),
    ('chrom.txt',      chrom),
    ('t_start.txt',    str(t_start)),
    ('t_end.txt',      str(t_end)),
    ('donor_name.txt', donor_name),
]:
    open(fname, 'w').write(val)

# sgRNAs.csv – effector hardcoded to SpCas9; scaffold/extension blank (non-PE)
with open('sgRNAs.csv', 'w') as f:
    f.write("name,effector,protospacer,scaffold,extension\n")
    f.write(f"{symbol},SpCas9,{guide_seq},,\n")

# Donor FASTA: full HDR template = left_HA + payload + right_HA.
# Using FASTA (not GenBank) as the pipeline does not pass GenBank as input.
with open('donor.fasta', 'w') as f:
    f.write(f">{donor_name}\n{donor_seq}\n")

print(f"Target:  {symbol} at {chrom}:{t_start}-{t_end}", file=sys.stderr)
print(f"Guide:   {guide_seq}", file=sys.stderr)
print(f"Donor:   {len(donor_seq)} bp  ({len(ctx_before)} ctx_5p + {len(left_ha)} HA_L + {len(payload)} payload + {len(right_ha)} HA_R + {len(ctx_after)} ctx_3p)", file=sys.stderr)
PYEOF

    SYMBOL=$(cat symbol.txt)
    CHROM=$(cat chrom.txt)
    T_START=$(cat t_start.txt)
    T_END=$(cat t_end.txt)
    DONOR_NAME=$(cat donor_name.txt)

    # ── 2. Extract local genomic context from reference ────────────────────────
    # This region serves as the stand-in "genome" for knock-knock build-strategies.
    # build-strategies aligns the protospacer into this FASTA, locates the cut
    # site, then extracts a 6 kb window around it to write target.gb.
    #
    # NOTE: distal-insertion detection is limited to sequences within this window
    # because no full-genome minimap2 index is provided. This has not been tested
    # in WGS mode (knock-knock upstream note).
    FLANK=~{flanking_bp}
    CTX_START=$(python3 -c "print(max(1, int('${T_START}') - ${FLANK}))")
    CTX_END=$(python3   -c "print(int('${T_END}') + ${FLANK})")

    echo "Extracting genomic context ${CHROM}:${CTX_START}-${CTX_END}" >&2

    # Rename the record to the gene symbol to avoid colons/hyphens in the name
    # that could cause issues with pyfaidx when building the genome index.
    samtools faidx ~{ref_fasta} "${CHROM}:${CTX_START}-${CTX_END}" \
      | awk -v name="${SYMBOL}" 'NR==1 { print ">" name; next } { print }' \
      > ctx.fa

    # ── 3. Assemble knock-knock working directory ──────────────────────────────
    # Strategy name for WGS mode (no amplicon_primers) is {genome}_{sgRNAs}
    # (see experiment.sample_sheet_row_to_editing_strategy_name).
    BASE="kk_base"
    BATCH="~{sample_id}_kk"

    mkdir -p "${BASE}/data/${BATCH}"
    mkdir -p "${BASE}/strategies/additional_sequences"
    # knock-knock looks for genome FASTAs in indices/{genome}/fasta/
    mkdir -p "${BASE}/indices/~{genome_name}/fasta"

    cp ctx.fa "${BASE}/indices/~{genome_name}/fasta/~{genome_name}.fa"
    samtools faidx "${BASE}/indices/~{genome_name}/fasta/~{genome_name}.fa"

    cp donor.fasta "${BASE}/strategies/additional_sequences/${DONOR_NAME}.fasta"

    cp sgRNAs.csv "${BASE}/strategies/sgRNAs.csv"

    # Copy region reads, stripping SAM tags (e.g. HP, PS) from FASTQ headers.
    # samtools fastq -T HP,PS appends tags as tab-separated fields in the read
    # name line; knock-knock's read name parser doesn't handle these.
    FASTQ_FN=$(basename ~{region_reads_fastq})
    zcat ~{region_reads_fastq} \
      | awk '{if(NR%4==1){sub(/\t.*/,"")}; print}' \
      | gzip > "${BASE}/data/${BATCH}/${FASTQ_FN}"

    # Sample sheet – WGS mode: amplicon_primers and primers both present but empty.
    # build_strategies.py requires 'amplicon_primers' (mandatory column check), while
    # sample_sheet_row_to_editing_strategy_name() accesses 'primers' to decide the
    # strategy name (empty → WGS naming: {genome}_{sgRNAs}).
    cat > "${BASE}/data/${BATCH}/sample_sheet.csv" <<CSVEOF
sample_name,platform,CCS_fastq_fn,genome,donor,amplicon_primers,primers,sgRNAs
~{sample_id},pacbio,${FASTQ_FN},~{genome_name},${DONOR_NAME},,,${SYMBOL}
CSVEOF

    cat "${BASE}/data/${BATCH}/sample_sheet.csv" >&2

    # ── 3b. Patch knock-knock: guard against reads with zero alignments ────────
    # Upstream bug: ReadDiagram.draw_reference() in visualize/architecture.py
    # indexes alignment_coordinates[0] in the 'centered on primers' branch
    # without a length check, crashing with IndexError on unaligned reads.
    # The conda site-packages dir is read-only in the container, so patch at
    # runtime via a sitecustomize.py on PYTHONPATH that monkey-patches the
    # method when the module is imported.
    KK_PATCH_DIR="$(pwd)/kk_patch"
    mkdir -p "${KK_PATCH_DIR}"
    cat > "${KK_PATCH_DIR}/sitecustomize.py" <<'PYEOF'
# Runtime monkey-patch for knock_knock.visualize.architecture.ReadDiagram.
# Wraps draw_reference so that an IndexError from alignment_coordinates[0]
# on a read with no alignments is caught and that ref_name is simply skipped.
import sys

def _install():
    try:
        from knock_knock.visualize import architecture as _arch
    except Exception:
        return
    RD = getattr(_arch, "ReadDiagram", None)
    if RD is None or getattr(RD.draw_reference, "_kk_patched", False):
        return
    _orig = RD.draw_reference
    def draw_reference(self, ref_name, ref_y, *args, **kwargs):
        try:
            return _orig(self, ref_name, ref_y, *args, **kwargs)
        except IndexError as e:
            print(f"[kk_patch] skipping draw_reference for {ref_name!r}: {e}",
                  file=sys.stderr)
            return None
    draw_reference._kk_patched = True
    RD.draw_reference = draw_reference

def _install_num_examples():
    # Raise default num_examples for generate_all_outcome_example_figures
    # from 10 → 100 so more per-category example diagrams are emitted.
    try:
        from knock_knock import experiment as _exp
    except Exception:
        return
    Exp = getattr(_exp, "Experiment", None)
    if Exp is None:
        return
    fn = getattr(Exp, "generate_all_outcome_example_figures", None)
    if fn is None or getattr(fn, "_kk_patched", False):
        return
    def generate_all_outcome_example_figures(self, num_examples=100, **kwargs):
        return fn(self, num_examples=num_examples, **kwargs)
    generate_all_outcome_example_figures._kk_patched = True
    Exp.generate_all_outcome_example_figures = generate_all_outcome_example_figures

def _hook(modname, installer):
    if modname in sys.modules:
        installer()
        return
    import importlib.abc, importlib.util
    class _Hook(importlib.abc.MetaPathFinder):
        def find_spec(self, name, path, target=None):
            if name == modname:
                sys.meta_path.remove(self)
                spec = importlib.util.find_spec(name)
                if spec is not None:
                    loader = spec.loader
                    _exec = loader.exec_module
                    def exec_module(module):
                        _exec(module)
                        installer()
                    loader.exec_module = exec_module
                return None
            return None
    sys.meta_path.insert(0, _Hook())

# Install now if already imported; otherwise hook import.
_install()
_hook("knock_knock.visualize.architecture", _install)
_install_num_examples()
_hook("knock_knock.experiment", _install_num_examples)
PYEOF

    export PYTHONPATH="${KK_PATCH_DIR}${PYTHONPATH:+:${PYTHONPATH}}"

    # ── 4. Build strategies ────────────────────────────────────────────────────
    # Reads sample sheet, aligns protospacer into ctx.fa, extracts the target
    # window, and writes strategies/{genome}_{symbol}/target.gb + parameters.yaml.
    knock-knock build-strategies "${BASE}" "${BATCH}"

    # ── 5. Run knock-knock ─────────────────────────────────────────────────────
    knock-knock parallel "${BASE}" ~{threads} --batch "${BATCH}"

    # ── 6. Collect outputs ─────────────────────────────────────────────────────
    RESULTS="${BASE}/results/${BATCH}/~{sample_id}"

    cp "${RESULTS}/outcome_list.txt"   "~{sample_id}_kk_outcome_list.txt"
    cp "${RESULTS}/all_lengths.png"    "~{sample_id}_kk_all_lengths.png"
    cp "${RESULTS}/outcome_counts.csv" "~{sample_id}_kk_outcome_counts.csv"

    # Tarballs of per-category files (paths relative to outcomes/ inside the archive)
    (cd "${RESULTS}" && find outcomes -name "diagrams.html" 2>/dev/null | sort > /tmp/kk_diagrams.list)
    if [ -s /tmp/kk_diagrams.list ]; then
      tar -czf "~{sample_id}_kk_diagrams.tar.gz" -C "${RESULTS}" -T /tmp/kk_diagrams.list
    else
      tar -czf "~{sample_id}_kk_diagrams.tar.gz" -T /dev/null
    fi

    (cd "${RESULTS}" && find outcomes -name "first_examples.png" 2>/dev/null | sort > /tmp/kk_png.list)
    if [ -s /tmp/kk_png.list ]; then
      tar -czf "~{sample_id}_kk_first_examples.tar.gz" -C "${RESULTS}" -T /tmp/kk_png.list
    else
      tar -czf "~{sample_id}_kk_first_examples.tar.gz" -T /dev/null
    fi
  >>>

  output {
    File outcome_list_txt              = "~{sample_id}_kk_outcome_list.txt"
    File edit_lengths_png              = "~{sample_id}_kk_all_lengths.png"
    File outcome_counts_csv            = "~{sample_id}_kk_outcome_counts.csv"
    File outcome_diagrams_html_tarball = "~{sample_id}_kk_diagrams.tar.gz"
    File outcome_png_tarball           = "~{sample_id}_kk_first_examples.tar.gz"
  }

  runtime {
    docker: "knock-knock"
    cpu: threads
    memory: "~{mem_gb} GiB"
    disk: disk_size + " GB"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
  }
}
