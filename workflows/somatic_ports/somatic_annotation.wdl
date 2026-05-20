version 1.1



## task definitons
task make_splice_region_bed {
    input {
        File ref_gff
        Array[String] extra_genes = ["BCOR"]
        Int distance = 500
        Int mem_gb = 8
    }

    command <<<
        set -euxo pipefail

        printf '%s\n' ~{sep(' ', extra_genes)} > extra_genes.txt

        python3 <<'PYEOF'
import gzip, sys

distance = ~{distance}

extra_genes = set()
with open('extra_genes.txt') as f:
    for line in f:
        g = line.strip()
        if g:
            extra_genes.add(g)

cgc_genes = set()
cgc_path = '/app/Compendium_Cancer_Genes.tsv'
with open(cgc_path) as f:
    header = f.readline().rstrip('\n').split('\t')
    symbol_col = header.index('SYMBOL')
    for line in f:
        parts = line.rstrip('\n').split('\t')
        if len(parts) > symbol_col and parts[symbol_col]:
            cgc_genes.add(parts[symbol_col])

target_genes = cgc_genes | extra_genes
sys.stderr.write(f"Targeting {len(target_genes)} genes ({len(cgc_genes)} CGC + {len(extra_genes)} extra)\n")

gff_path = '~{ref_gff}'
opener = gzip.open if gff_path.endswith('.gz') else open
regions = []
with opener(gff_path, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.rstrip('\n').split('\t')
        if len(parts) < 9 or parts[2] != 'gene':
            continue
        attrs = {}
        for kv in parts[8].split(';'):
            if '=' in kv:
                k, v = kv.split('=', 1)
                attrs[k.strip()] = v.strip()
        gene_name = attrs.get('Name') or attrs.get('gene_name') or attrs.get('gene_id', '')
        if gene_name in target_genes:
            chrom = parts[0]
            start = max(0, int(parts[3]) - 1 - distance)  # GFF3 1-based → BED 0-based, plus padding
            end   = int(parts[4]) + distance
            regions.append((chrom, start, end, gene_name))

regions.sort(key=lambda x: (x[0], x[1]))
sys.stderr.write(f"Writing {len(regions)} gene regions to BED\n")
with open('splice_regions.bed', 'w') as f:
    for chrom, start, end, name in regions:
        f.write(f'{chrom}\t{start}\t{end}\t{name}\n')
PYEOF
    >>>

    output {
        File regions_bed = "splice_regions.bed"
    }

    runtime {
        docker: "quay.io/pacbio/somatic_general_tools@sha256:a25a2e62b88c73fa3c18a0297654420a4675224eb0cf39fa4192f8a1e92b30d6"
        memory: "~{mem_gb} GB"
        disk: ceil(size(ref_gff, "GB") + 1) + " GB"
        maxRetries: 2
        preemptible: 1
    }
}

task filter_vcf_by_regions {
    input {
        File input_vcf
        File input_vcf_index
        File regions_bed
        String out_prefix
        Int mem_gb = 4
    }

    Int disk_size = ceil(size(input_vcf, "GB") * 2 + 1)

    command <<<
        set -euxo pipefail
        bcftools --version

        bcftools view \
            --regions-file ~{regions_bed} \
            --output-type z \
            --output ~{out_prefix}.vcf.gz \
            ~{input_vcf}

        bcftools index --tbi ~{out_prefix}.vcf.gz
    >>>

    output {
        File filtered_vcf       = "~{out_prefix}.vcf.gz"
        File filtered_vcf_index = "~{out_prefix}.vcf.gz.tbi"
    }

    runtime {
        docker: "quay.io/pacbio/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87"
        memory: "~{mem_gb} GB"
        disk: disk_size + " GB"
        maxRetries: 2
        preemptible: 1
    }
}

task spliceai_annotate {
    input {
        File input_vcf
        File input_vcf_index
        File ref_fasta
        File ref_fasta_index
        String genome = "grch38"
        Int distance = 500
        Int threads = 64
        Int mem_gb = 64
    }

    String out_prefix = sub(basename(input_vcf), "\\.vcf\\.gz$", "")
    Float file_size = ceil(size(input_vcf, "GB") + size(ref_fasta, "GB") + 5)

    command <<<
        set -euxo pipefail

        spliceai \
            -I ~{input_vcf} \
            -O ~{out_prefix}.spliceai.vcf \
            -R ~{ref_fasta} \
            -A ~{genome} \
            -D ~{distance}

    >>>

    output {
        File spliceai_vcf = "~{out_prefix}.spliceai.vcf"
    }

    runtime {
        docker: "quay.io/biocontainers/spliceai@sha256:dcbe88faa015c5a92490b25d54a3fda981d0efaf73821f82ed088d8261754c78"
        cpu: threads
        memory: "~{mem_gb} GB"
        disk: file_size + " GB"
        maxRetries: 2
        preemptible: 1
    }
}

task bgzip_and_index_vcf {
    input {
        File input_vcf
        Int mem_gb = 4
    }

    String out_vcf = basename(input_vcf) + ".gz"
    Int disk_size  = ceil(size(input_vcf, "GB") * 2 + 1)

    command <<<
        set -euxo pipefail
        bgzip --version

        bgzip -c ~{input_vcf} > ~{out_vcf}
        tabix ~{out_vcf}
    >>>

    output {
        File vcf       = "~{out_vcf}"
        File vcf_index = "~{out_vcf}.tbi"
    }

    runtime {
        docker: "quay.io/pacbio/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87"
        memory: "~{mem_gb} GB"
        disk: disk_size + " GB"
        maxRetries: 2
        preemptible: 1
    }
}

task bcftools_annotate_spliceai {
    input {
        File input_vcf
        File input_vcf_index
        File spliceai_vcf
        File spliceai_vcf_index
        Int mem_gb = 8
    }

    String out_prefix = sub(basename(input_vcf), "\\.vcf\\.gz$", "")
    Int disk_size = ceil(size(input_vcf, "GB") * 2 + size(spliceai_vcf, "GB") + 2)

    command <<<
        set -euxo pipefail

        # Parse SpliceAI VCF: for each variant pick the entry with the highest
        # max DS score across all overlapping genes, output as annotation TSV
        python3 <<'PYEOF'
import gzip, sys

def best_spliceai(info):
    for field in info.split(';'):
        if not field.startswith('SpliceAI='):
            continue
        best, best_score = None, -1.0
        for entry in field[9:].split(','):
            parts = entry.split('|')
            if len(parts) < 6:
                continue
            try:
                score = max(float(x) for x in parts[2:6])
            except ValueError:
                continue
            if score > best_score:
                best_score, best = score, parts
        return best
    return None

opener = gzip.open if '~{spliceai_vcf}'.endswith('.gz') else open
with opener('~{spliceai_vcf}', 'rt') as fin, open('spliceai_annot.tsv', 'w') as fout:
    fout.write('#CHROM\tPOS\tREF\tALT\tSpliceAI_DS_AG\tSpliceAI_DS_AL\tSpliceAI_DS_DG\tSpliceAI_DS_DL\tSpliceAI_SYMBOL\n')
    for line in fin:
        if line.startswith('#'):
            continue
        cols = line.rstrip('\n').split('\t')
        chrom, pos, _, ref, alt, _, _, info = cols[:8]
        best = best_spliceai(info)
        if best:
            fout.write(f'{chrom}\t{pos}\t{ref}\t{alt}\t{best[2]}\t{best[3]}\t{best[4]}\t{best[5]}\t{best[1]}\n')
PYEOF

        bgzip spliceai_annot.tsv
        tabix -s1 -b2 -e2 spliceai_annot.tsv.gz

        printf '##INFO=<ID=SpliceAI_DS_AG,Number=1,Type=Float,Description="SpliceAI delta score: acceptor gain">\n'   >  spliceai_header.txt
        printf '##INFO=<ID=SpliceAI_DS_AL,Number=1,Type=Float,Description="SpliceAI delta score: acceptor loss">\n'   >> spliceai_header.txt
        printf '##INFO=<ID=SpliceAI_DS_DG,Number=1,Type=Float,Description="SpliceAI delta score: donor gain">\n'      >> spliceai_header.txt
        printf '##INFO=<ID=SpliceAI_DS_DL,Number=1,Type=Float,Description="SpliceAI delta score: donor loss">\n'      >> spliceai_header.txt
        printf '##INFO=<ID=SpliceAI_SYMBOL,Number=1,Type=String,Description="SpliceAI gene symbol">\n'                >> spliceai_header.txt

        bcftools annotate \
            --annotations spliceai_annot.tsv.gz \
            --header-lines spliceai_header.txt \
            --columns CHROM,POS,REF,ALT,INFO/SpliceAI_DS_AG,INFO/SpliceAI_DS_AL,INFO/SpliceAI_DS_DG,INFO/SpliceAI_DS_DL,INFO/SpliceAI_SYMBOL \
            --output-type z \
            --output ~{out_prefix}.spliceai_annot.vcf.gz \
            ~{input_vcf}

        tabix ~{out_prefix}.spliceai_annot.vcf.gz
    >>>

    output {
        File annotated_vcf       = "~{out_prefix}.spliceai_annot.vcf.gz"
        File annotated_vcf_index = "~{out_prefix}.spliceai_annot.vcf.gz.tbi"
    }

    runtime {
        docker: "quay.io/pacbio/somatic_general_tools@sha256:a25a2e62b88c73fa3c18a0297654420a4675224eb0cf39fa4192f8a1e92b30d6"
        memory: "~{mem_gb} GB"
        disk: disk_size + " GB"
        maxRetries: 2
        preemptible: 1
    }
}

task vep_annotate {
    input {
        File input_vcf
        File? vep_cache
        File ref_fasta
        File ref_fasta_index
        File? dbnsfp_file
        File? dbnsfp_file_index
        File? clinvar_vcf
        File? clinvar_vcf_index
        Int threads = 16
        Int mem_gb  = 16
    }

    Float file_size = ceil(size(input_vcf, "GB") + size(vep_cache, "GB") + size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(dbnsfp_file, "GB") + size(clinvar_vcf, "GB") + 10)
    String vep_annotated_vcf = sub(basename(input_vcf), "\\.vcf.gz$", "") + ".vep.vcf.gz"

    command <<<
        set -euxo pipefail

        mkdir -p vep_data/
        if [ ! -f ~{vep_cache} ]; then
            echo "VEP cache file not found. Please provide a valid cache file."
            exit 1
        fi

        vep --help

        tar -xzvf ~{vep_cache} -C vep_data/

        # Build optional plugin/custom args
        extra_args=""
        if [ -f "~{dbnsfp_file}" ]; then
            extra_args="${extra_args} --plugin dbNSFP,~{dbnsfp_file},BayesDel_noAF_score,REVEL_score,CADD_phred,AlphaMissense_score"
        fi
        if [ -f "~{clinvar_vcf}" ]; then
            extra_args="${extra_args} --custom ~{clinvar_vcf},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN"
        fi
        vep \
            --cache \
            --offline \
            --dir vep_data/ \
            --fasta ~{ref_fasta} \
            --format vcf \
            --fork ~{threads} \
            --species homo_sapiens \
            --assembly GRCh38 \
            --symbol \
            --hgvs \
            --refseq \
            --check_existing \
            --vcf \
            --pick \
            --flag_pick_allele_gene \
            --everything \
            --plugin SpliceRegion \
            --compress_output bgzip \
            ${extra_args} \
            -i ~{input_vcf} \
            -o ~{vep_annotated_vcf}

        rm -rf vep_data/

        tabix "~{vep_annotated_vcf}"
    >>>

    output {
        File annotated_vcf = "~{vep_annotated_vcf}"
        File annotated_vcf_index = "~{vep_annotated_vcf}.tbi"
    }

    runtime {
        docker: "ensemblorg/ensembl-vep@sha256:ff3c7e20d68e7e499c0bd79d398c7c121465db65f311ee87df425b2b600b853e"
        cpu: threads
        memory: "~{mem_gb} GB"
        disk: file_size + " GB"
        maxRetries: 2
        preemptible: 1
    }
}

task annotsv {
    input {
        File sv_vcf
        File sv_vcf_index
        File? annotsv_cache
        Int threads = 2
    }

    Float file_size = ceil(size(sv_vcf, "GB") + size(annotsv_cache, "GB"))
    #String annotsv_annotated_vcf = sub(basename(sv_vcf), "\\.vcf.gz$", "") + ".annotsv.vcf"
    String annotsv_annotated_tsv = sub(basename(sv_vcf), "\\.vcf.gz$", "") + ".annotsv.tsv"

    command <<<
        set -euxo pipefail

        # Process VCF to move BND alt to INFO field
        # Check if it ends in gz. Unzip it if it's the case
        if [[ ~{sv_vcf} == *.gz ]]; then
            gunzip -c ~{sv_vcf} > tmp.vcf
        else
            cp ~{sv_vcf} tmp.vcf
        fi

        awk -F'\t' -v OFS='\t' '
        { if (NR==2) {
                print "##INFO=<ID=SV_ALT,Number=1,Type=String,Description=\"Square bracketed notation for BND event\">"
            }
        }
        {
            if ($0 ~ /^#/) {
                print $0;  # Print header lines as is
            } else {
                if ($8 ~ /SVTYPE=BND/) {
                    $8 = $8 ";SV_ALT=" $5;  # Append ALT to INFO as SV_ALT
                    $5 = "<BND>";  # Change ALT to <BND>
                }
                print $0;
            }
        }' tmp.vcf > tmp_processed.vcf

        mkdir -p annotsv_cache_dir
        # If annotsv_cache is not provided, fail
        if [ ! -f ~{annotsv_cache} ]; then
            echo "AnnotSV cache file not found. Please provide a valid cache file."
            exit 1
        fi

        AnnotSV --version

        tar -xzvf ~{annotsv_cache} -C annotsv_cache_dir/
        AnnotSV_annotations=$(find annotsv_cache_dir/ -type d -name AnnotSV_annotations -print -quit)

        AnnotSV \
            -annotationsDir "${AnnotSV_annotations}" \
            -annotationMode both \
            -SVinputFile tmp_processed.vcf \
            -outputDir . \
            -outputFile "~{annotsv_annotated_tsv}" \
            -SVinputInfo 1 \
            -genomeBuild GRCh38

        # Delete cache after annotation
        rm -rf annotsv_cache_dir/
    >>>

    output {
        File annotated_tsv = "~{annotsv_annotated_tsv}"
    }

    runtime {
        docker: "quay.io/biocontainers/annotsv@sha256:0c73fef5fa529b11e10bea0355480f01b56d0feb21af54cb9bbbd1f9f4c862a7"
        cpu: threads
        memory: "~{threads * 2} GB"
        disk: file_size + " GB"
        maxRetries: 2
        preemptible: 1
    }
}

task chord_hrd {
    input {
        File small_variant_vcf
        File sv_vcf
        String pname
        Int threads = 4
    }

    Float file_size = ceil(size(small_variant_vcf, "GB") + size(sv_vcf, "GB"))

    command <<<
    set -euxo pipefail

    # Docker image uses GRIDSS as default. Change to Manta
    sed 's/gridss/manta/g' \
        /opt/chord/extractSigPredictHRD.R > ./extractSigPredictHRD.R
    chmod +x ./extractSigPredictHRD.R

    ./extractSigPredictHRD.R . ~{pname} ~{small_variant_vcf} ~{sv_vcf} 38 2>&1 | tee chord_hrd.log
    rm -f ./extractSigPredictHRD.R
    >>>

    output {
        File chord_log = "chord_hrd.log"
        File chord_prediction = pname + "_chord_prediction.txt"
        File chord_signature = pname + "_chord_signatures.txt"
    }

    runtime {
        docker: "scwatts/chord@sha256:9f6aa44ffefe3f736e66a0e2d7941d4f3e1cc6d848a9a11a17e85a6525e63a77"
        cpu: threads
        memory: "~{threads * 4} GB"
        disk: file_size + " GB"
        maxRetries: 2
        preemptible: 1
    }
}



task prioritize_sv_intogen {
    input {
        File annotSV_tsv
        Int threads = 2
    }

    Float file_size = ceil(size(annotSV_tsv, "GB") + 10)
    String CCG_tsv_name = "~{sub(basename(annotSV_tsv), "\\.tsv$", "")}_CCG.tsv"
    String CCG_ranked_tsv_name = "~{sub(basename(annotSV_tsv), "\\.tsv$", "")}_CCG.ranked.tsv"
    String ranked_tsv_name = "~{sub(basename(annotSV_tsv), "\\.tsv$", "")}.ranked.tsv"
    String ranked_concerning_tsv_name = "~{sub(basename(annotSV_tsv), "\\.tsv$", "")}.ranked_concerning.tsv"

    command <<<
    set -euxo pipefail

    csvtk version

    # Remove any quote from the file
    sed "s/\"//g" ~{annotSV_tsv} > ~{basename(annotSV_tsv)}_noquote.tsv

    # 'split' is gene-centric view; 'full' is variant-centric.
    csvtk filter2 -t -f '$Annotation_mode == "full"' \
      ~{basename(annotSV_tsv)}_noquote.tsv > ~{basename(annotSV_tsv)}_full.tsv
    csvtk filter2 -t -f '$Annotation_mode == "split"' \
      ~{basename(annotSV_tsv)}_noquote.tsv > ~{basename(annotSV_tsv)}_split.tsv

    csvtk join --left-join -t \
      ~{basename(annotSV_tsv)}_split.tsv \
      /app/Compendium_Cancer_Genes.tsv \
      -f "Gene_name;SYMBOL" |\
        csvtk summary -t -g "$(csvtk headers -t ~{annotSV_tsv} | tr '\n' ',' | sed 's/,$//g')" \
            -f CANCER_TYPE:collapse,COHORT:collapse,TRANSCRIPT:collapse,MUTATIONS:collapse,ROLE:collapse,CGC_GENE:collapse,CGC_CANCER_GENE:collapse,DOMAINS:collapse,2D_CLUSTERS:collapse,3D_CLUSTERS:collapse -s ";" |\
        sed 's/:collapse//g' > ~{CCG_tsv_name}

    rm -f ~{basename(annotSV_tsv)}_noquote.tsv

    # sort by consequence score and prune columns
    # TODO: use named column finding
    (head -1 ~{basename(annotSV_tsv)}_full.tsv &&
     tail -n +2 ~{basename(annotSV_tsv)}_full.tsv | sort -t$'\t' -s -k119,119nr -k2,2n -k3,3n
    ) > ~{ranked_tsv_name}

    (head -1 ~{CCG_tsv_name} &&
     tail -n +2 ~{CCG_tsv_name} | sort -t$'\t' -s -k119,119nr -k2,2n -k3,3n
    ) > ~{CCG_ranked_tsv_name}

    rm -f ~{basename(annotSV_tsv)}_full.tsv
    rm -f ~{basename(annotSV_tsv)}_split.tsv

    # ACMG class > VOUS (3) (likely pathogenic (4) or pathogenic (5))
    keep_columns='{print $1,$2,$3,$5,$6,$9,$10,$14,$15,$18,$27,$28,$29,$97,$98,$105,$106,$107,$108,$111,$112,$116,$119,$120,$121}'
    acmg_col='$121'
    (head -1 ~{ranked_tsv_name} | awk -F '\t' -vOFS='\t' "${keep_columns}" &&
     tail -n +2 ~{ranked_tsv_name} | awk -F '\t' -vOFS='\t' "${acmg_col} != \"NA\" && ${acmg_col} > 3 ${keep_columns}"
    ) > ~{ranked_concerning_tsv_name}
    >>>

    output {
        File annotSV_ranked_tsv = "~{ranked_tsv_name}"
        File annotSV_concerning_tsv = "~{ranked_concerning_tsv_name}"
        File annotSV_ccg_tsv = "~{CCG_tsv_name}"
        File annotSV_ccg_ranked_tsv = "~{CCG_ranked_tsv_name}"
    }

    runtime {
        docker: "quay.io/pacbio/somatic_general_tools@sha256:a25a2e62b88c73fa3c18a0297654420a4675224eb0cf39fa4192f8a1e92b30d6"
        cpu: threads
        memory: "~{threads * 2} GB"
        disk: file_size + " GB"
        maxRetries: 2
        preemptible: 1
    }
}

task prioritize_small_variants {
    input {
        File vep_annotated_vcf
        Int threads = 2
        String sample
        String pname="sample"
    }

    Float file_size = ceil(size(vep_annotated_vcf, "GB") + 10)
    String fname = sub(basename(vep_annotated_vcf), "\\.vcf.gz", "") + ".tsv"
    String fname2 = sub(basename(vep_annotated_vcf), "\\.vcf.gz", "") + "_intogenCCG.tsv"
    String fname_ranked = sub(basename(vep_annotated_vcf), "\\.vcf.gz", "") + "_ranked.tsv"
    String fname_filtered = sub(basename(vep_annotated_vcf), "\\.vcf.gz", "") + "_filtered.tsv"

    command <<<
    set -euxo pipefail

    csvtk version

    echo -e "CHROM\tPOS\tREF\tALT\tFORMAT\t~{pname}\t$(bcftools +split-vep ~{vep_annotated_vcf} -l | cut -f2 | tr '\n' '\t' | sed 's/\t$//g')\tSpliceAI_DS_AG\tSpliceAI_DS_AL\tSpliceAI_DS_DG\tSpliceAI_DS_DL\tSpliceAI_SYMBOL" > ~{fname}
    bcftools view -s ~{sample} -e 'GT=".|." || GT="./." || GT="." || GT="ref"' ~{vep_annotated_vcf} \
      | bcftools +split-vep -A tab -f '%CHROM\t%POS\t%REF\t%ALT\t%FORMAT\t%CSQ\t%INFO/SpliceAI_DS_AG\t%INFO/SpliceAI_DS_AL\t%INFO/SpliceAI_DS_DG\t%INFO/SpliceAI_DS_DL\t%INFO/SpliceAI_SYMBOL\n' >> ~{fname}

    #RANKING AND FILTERING
    [[ -f "~{fname}" ]] || { echo "ERROR: ${fname} not found" >&2; exit 1; }

    # 2) find IMPACT column index (1‐based)
    IFS=$'\t' read -r -a hdr < <(head -n1 "~{fname}")
    impactcol=0
    for i in "${!hdr[@]}"; do
      [[ "${hdr[i]}" == "IMPACT" ]] && { impactcol=$((i+1)); break; }
    done
    [[ $impactcol -gt 0 ]] || { echo "ERROR: IMPACT column not found" >&2; exit 1; }

    # 3) write header unchanged
    head -n1 "~{fname}" > "~{fname_ranked}"

    # 4) rank & sort, then strip the rank in awk (no cut!)
    tail -n +2 "~{fname}" \
      | awk -F $'\t' -v OFS=$'\t' -v impact="$impactcol" '
          {
            # compute numeric rank
            r = ($impact=="HIGH"     ? 1 : $impact=="MODERATE" ? 2 : $impact=="LOW"      ? 3 : $impact=="MODIFIER" ? 4 : 5)
            # prepend it to the line
            print r, $0
          }
        ' \
      | sort -t $'\t' -k1,1n \
      | awk -F $'\t' -v OFS=$'\t' '
          {
            # drop the first field (the rank) and re-print the rest
            for(i=2; i<=NF; i++){
              printf("%s%s", $i, i<NF?OFS:ORS)
            }
          }
        ' \
      >> "~{fname_ranked}"
    # 5) find MAX_AF column index (1‐based)
    IFS=$'\t' read -r -a hdr2 < <(head -n1 "~{fname_ranked}")
    maxafcol=0
    for i in "${!hdr2[@]}"; do
      [[ "${hdr2[i]}" == "MAX_AF" ]] && { maxafcol=$((i+1)); break; }
    done
    [[ $maxafcol -gt 0 ]] || { echo "ERROR: MAX_AF column not found" >&2; exit 1; }

    # 6) filter HIGH or MODERATE impact and MAX_AF < 0.001 (Cite: 10.1038/s41525-021-00227-3)
    awk -F $'\t' -v OFS=$'\t' -v impact="$impactcol" -v maxaf="$maxafcol" '
      NR==1 { print; next }
      ($impact=="HIGH" || $impact=="MODERATE") && ($maxaf == "" || $maxaf < 0.001)
    ' "~{fname_ranked}" > "~{fname_filtered}"
    >>>

    output {
        File vep_annotated_tsv = fname
        File vep_annotated_ranked_tsv = fname_ranked
        File vep_annotated_filtered_tsv = fname_filtered
        #File vep_annotated_tsv_intogenCCG = fname2
    }

    runtime {
        docker: "quay.io/pacbio/somatic_general_tools@sha256:a25a2e62b88c73fa3c18a0297654420a4675224eb0cf39fa4192f8a1e92b30d6"
        cpu: threads
        memory: "~{threads * 2} GB"
        disk: file_size + " GB"
        maxRetries: 2
        preemptible: 1
    }
}
