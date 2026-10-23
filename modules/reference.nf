// REPLACED (was: transcript-header based tables, ref in data/ref). Now: ref/ folder, md5 check, tables from the GTF.
// Genome = filtered primary-chromosome hg38 FASTA (user decision 2026-10-08). Written from tool documentation (Salmon "Preparing transcriptome indices (mapping-based mode)", GENCODE FTP layout).
process REF_FETCH {
    label 'process_low'
    executor 'local'
    publishDir "${params.ref_dir}/source", mode: 'copy'

    output:
    path "gencode.${params.gencode_release}.transcripts.fa.gz", emit: transcripts
    path "gencode.${params.gencode_release}.primary_assembly.annotation.gtf.gz", emit: gtf
    path "md5_check.txt"

    script:
    def base = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${params.gencode_release.replace('v','')}"
    """
    curl -sS -L --retry 5 --retry-delay 10 -o MD5SUMS ${base}/MD5SUMS
    : > md5_check.txt
    for f in gencode.${params.gencode_release}.transcripts.fa.gz gencode.${params.gencode_release}.primary_assembly.annotation.gtf.gz; do
        curl -sS -L --retry 5 --retry-delay 10 -C - -o \$f ${base}/\$f
        gzip -t \$f
        exp=\$(awk -v f=\$f '\$2==f || \$2=="./"f {print \$1}' MD5SUMS | head -1)
        got=\$(md5sum \$f | cut -d' ' -f1)
        if [ -z "\$exp" ]; then echo "\$f md5_not_listed got=\$got" >> md5_check.txt
        elif [ "\$exp" = "\$got" ]; then echo "\$f md5_OK \$got" >> md5_check.txt
        else echo "\$f md5_MISMATCH expected=\$exp got=\$got" >> md5_check.txt; cat md5_check.txt >&2; exit 1; fi
    done
    cat md5_check.txt
    """

    stub:
    """
    touch gencode.${params.gencode_release}.transcripts.fa.gz gencode.${params.gencode_release}.primary_assembly.annotation.gtf.gz
    echo "stub" > md5_check.txt
    """
}

process REF_FILTER {
    label 'process_low'
    executor 'local'
    publishDir "${params.ref_dir}", mode: 'copy'

    input:
    path gtf
    path transcripts

    output:
    path "gencode.${params.gencode_release}.primary_chr.annotation.gtf.gz", emit: gtf
    path "gencode.${params.gencode_release}.transcripts.primary_chr.fa.gz", emit: transcripts
    path "ref_filter_counts.txt"
    path "removed_par_y_genes.txt"

    script:
    // One gene set for STAR, featureCounts and Salmon (tx2gene is built from this GTF):
    //  1. GTF restricted to chr1-22,X,Y,M   2. chrY genes lying entirely inside PAR1/PAR2 removed (chrX copies carry the PAR counts;
    //  the genome has these chrY regions hard-masked)   3. transcript FASTA restricted to the transcript IDs of the resulting GTF
    """
    zcat ${gtf} | awk -F'\\t' '/^#/ || \$1 ~ /^chr([0-9]+|X|Y|M)\$/' > chr.gtf
    awk -F'\\t' '\$1=="chrY" && \$3=="gene" && ((\$4>=${params.par1_start} && \$5<=${params.par1_end}) || (\$4>=${params.par2_start} && \$5<=${params.par2_end})) {match(\$9,/gene_id "[^"]+"/); print substr(\$9,RSTART+9,RLENGTH-10)}' chr.gtf | sort -u > removed_par_y_genes.txt
    awk -F'\\t' 'BEGIN{while((getline l < "removed_par_y_genes.txt")>0) r[l]=1} /^#/{print; next} {match(\$9,/gene_id "[^"]+"/); g=substr(\$9,RSTART+9,RLENGTH-10); if(!(g in r)) print}' chr.gtf | gzip > gencode.${params.gencode_release}.primary_chr.annotation.gtf.gz
    zcat gencode.${params.gencode_release}.primary_chr.annotation.gtf.gz | awk -F'\\t' '\$3=="transcript"' | sed -E 's/.*transcript_id "([^"]+)".*/\\1/' | sort -u > keep_ids.txt
    zcat ${transcripts} | awk 'BEGIN{while((getline l < "keep_ids.txt")>0) k[l]=1} /^>/{split(substr(\$0,2),a,"|"); w=(a[1] in k)} w' | gzip > gencode.${params.gencode_release}.transcripts.primary_chr.fa.gz
    {
      echo "GTF lines in: \$(zcat ${gtf} | grep -vc '^#'); after chr filter: \$(grep -vc '^#' chr.gtf); after PAR-Y removal: \$(zcat gencode.${params.gencode_release}.primary_chr.annotation.gtf.gz | grep -vc '^#')"
      echo "PAR-Y genes removed: \$(wc -l < removed_par_y_genes.txt)"
      echo "transcript FASTA in: \$(zcat ${transcripts} | grep -c '^>'); kept: \$(zcat gencode.${params.gencode_release}.transcripts.primary_chr.fa.gz | grep -c '^>'); GTF transcripts: \$(wc -l < keep_ids.txt)"
    } | tee ref_filter_counts.txt
    """

    stub:
    """
    touch gencode.${params.gencode_release}.primary_chr.annotation.gtf.gz gencode.${params.gencode_release}.transcripts.primary_chr.fa.gz removed_par_y_genes.txt
    echo stub > ref_filter_counts.txt
    """
}

process REF_GENOME_FILTER {
    label 'process_low'
    executor 'local'
    publishDir "${params.ref_dir}", mode: 'copy'

    output:
    path "GRCh38_primary_chr.fa", emit: genome

    script:
    // read-only use of params.genome_source; keeps chr1-22, X, Y, M and rewrites headers to the bare contig name
    """
    awk '/^>/{keep=(\$1 ~ /^>chr([0-9]+|X|Y|M)\$/); if(keep) print \$1; next} keep' ${params.genome_source} > GRCh38_primary_chr.fa
    test \$(grep -c '^>' GRCh38_primary_chr.fa) -eq 25
    """

    stub:
    """
    printf ">chr1\\nACGT\\n" > GRCh38_primary_chr.fa
    """
}

process REF_TABLES {
    label 'process_low'
    publishDir "${params.ref_dir}", mode: 'copy'

    input:
    path gtf

    output:
    path "tx2gene.tsv", emit: tx2gene
    path "gene_id_to_symbol.tsv", emit: symbols

    script:
    // transcript lines -> (transcript_id, gene_id); gene lines -> (gene_id, gene_name, gene_type); versions kept (same IDs as featureCounts/Salmon)
    """
    zcat ${gtf} | awk -F'\\t' '\$3=="transcript"' | \\
      sed -E 's/.*gene_id "([^"]+)".*transcript_id "([^"]+)".*/\\2\\t\\1/' > tx2gene.tsv
    printf "gene_id\\tgene_name\\tgene_type\\n" > gene_id_to_symbol.tsv
    zcat ${gtf} | awk -F'\\t' '\$3=="gene"' | \\
      sed -E 's/.*gene_id "([^"]+)".*gene_type "([^"]+)".*gene_name "([^"]+)".*/\\1\\t\\3\\t\\2/' >> gene_id_to_symbol.tsv
    echo "tx2gene rows: \$(wc -l < tx2gene.tsv); genes: \$(( \$(wc -l < gene_id_to_symbol.tsv) - 1 ))"
    """

    stub:
    """
    printf "ENST1\\tENSG1\\n" > tx2gene.tsv
    printf "gene_id\\tgene_name\\tgene_type\\nENSG1\\tG1\\tprotein_coding\\nGENE1\\tGENE1\\tprotein_coding\\nGENE2\\tGENE2\\tprotein_coding\\n" > gene_id_to_symbol.tsv
    """
}

process SALMON_INDEX {
    label 'process_high'
    conda "${params.install_dir}/env/prefix/salmon"
    publishDir "${params.ref_dir}", mode: 'copy'

    input:
    path transcripts
    path genome

    output:
    path "salmon_index", emit: index
    path "decoys.txt"

    script:
    // decoy-aware (Salmon docs): decoys = genome sequence names, gentrome = transcripts followed by genome
    """
    grep '^>' ${genome} | cut -d ' ' -f 1 | sed 's/>//' > decoys.txt
    ( zcat ${transcripts}; cat ${genome} ) > gentrome.fa
    salmon index -t gentrome.fa -d decoys.txt -i salmon_index -p ${task.cpus} --gencode -k 31
    rm gentrome.fa
    """

    stub:
    """
    mkdir salmon_index && touch salmon_index/info.json decoys.txt
    """
}
