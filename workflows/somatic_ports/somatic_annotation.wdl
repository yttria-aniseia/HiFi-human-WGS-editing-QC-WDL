version 1.1



## task definitons
task vep_annotate {
    input {
        File input_vcf
        File? vep_cache
        File ref_fasta
        File ref_fasta_index
        Int threads = 16
        Int mem_gb  = 8
    }

    Float file_size = ceil(size(input_vcf, "GB") + size(vep_cache, "GB") + size(ref_fasta, "GB") + size(ref_fasta_index, "GB"))
    String vep_annotated_vcf = sub(basename(input_vcf), "\\.vcf.gz$", "") + ".vep.vcf.gz"

    command <<<
        set -euxo pipefail

        mkdir -p vep_data/
        # If vep_cache is not provided, fail
        if [ ! -f ~{vep_cache} ]; then
            echo "VEP cache file not found. Please provide a valid cache file."
            exit 1
        fi

        vep --help

        tar -xzvf ~{vep_cache} -C vep_data/
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
            --compress_output bgzip \
            -i ~{input_vcf} \
            -o ~{vep_annotated_vcf}

        # Delete cache after annotation
        rm -rf vep_data/

        tabix "~{vep_annotated_vcf}"
    >>>

    output {
        File annotated_vcf = "~{vep_annotated_vcf}"
        File annotated_vcf_index = "~{vep_annotated_vcf}.tbi"
    }

    runtime {
        docker: "ensemblorg/ensembl-vep@sha256:e7612ab7c2923f2b9a78592b939e74874cd29f7494d70ee7135c8303841b03a8"
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
        docker: "quay.io/biocontainers/annotsv@sha256:d814365a1d81bca2b4690a6a8af316881587cadb38b3bb524ec8747380bc7a36"
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

    echo -e "CHROM\tPOS\tREF\tALT\tFORMAT\t~{pname}\t$(bcftools +split-vep ~{vep_annotated_vcf} -l | cut -f2 | tr '\n' '\t' | sed 's/\t$//g')" > ~{fname}
     bcftools view -s ~{sample} -e 'GT="ref"||GT="mis"' ~{vep_annotated_vcf} | bcftools +split-vep -A tab -f '%CHROM\t%POS\t%REF\t%ALT\t%FORMAT\t%CSQ\n' >> ~{fname}

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
    # 5) filter HIGH or MODERATE into the final file
    awk -F $'\t' -v OFS=$'\t' -v impact="$impactcol" '
      NR==1 { print; next }
      ($impact=="HIGH" || $impact=="MODERATE")
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
