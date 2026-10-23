#!/usr/bin/env bash
# Runs every tool once with its version flag / loads every R package, writes the REAL installed versions to docs/tool_versions.tsv
B=/workspace/Biomarker_pipeline; P=$B/env/prefix; OUT=$B/docs/tool_versions.tsv
export MAMBA_ROOT_PREFIX=$B/env/mamba TMPDIR=$B/tmp
v() { tool=$1; how=$2; shift 2; ver=$("$@" 2>&1 | grep -v '^\s*$' | head -${N:-1} | tr '\n' ' ' | sed 's/[[:space:]]*$//'); printf "%s\t%s\t%s\n" "$tool" "$ver" "$how"; }
{
printf "tool\tversion_reported_by_tool\thow_installed\n"
v micromamba "static binary env/bin" $B/env/bin/micromamba --version
v java "micromamba env nf (conda-forge)" $B/env/bin/micromamba run -n nf java -version
v nextflow "micromamba env nf (bioconda)" $B/env/bin/micromamba run -n nf nextflow -version
v fastp "env/prefix/fastp from env/fastp.yml" $P/fastp/bin/fastp --version
v fastqc "env/prefix/fastqc from env/fastqc.yml (openjdk 25.0.2 from the same env)" env PATH=$P/fastqc/bin:$PATH $P/fastqc/bin/fastqc --version
v multiqc "env/prefix/multiqc from env/multiqc.yml (pip multiqc==1.35 inside python 3.12 env)" $P/multiqc/bin/multiqc --version
N=1 v samtools "env/prefix/samtools from env/samtools.yml" $P/samtools/bin/samtools --version
v STAR "env/prefix/star from env/star.yml" $P/star/bin/STAR --version
v featureCounts "env/prefix/subread from env/subread.yml" $P/subread/bin/featureCounts -v
v salmon "env/prefix/salmon from env/salmon.yml" $P/salmon/bin/salmon --version
v python "env/prefix/pyqc from env/pyqc.yml" $P/pyqc/bin/python --version
v pandas "env/prefix/pyqc" $P/pyqc/bin/python -c 'import pandas; print(pandas.__version__)'
v R "env/prefix/r_dge from env/r_dge.yml" $P/r_dge/bin/Rscript -e 'cat(R.version.string)'
for pk in DESeq2 tximport apeglm ggplot2 pheatmap ggrepel jsonlite readr; do
  v "R:$pk" "env/prefix/r_dge" $P/r_dge/bin/Rscript -e "suppressPackageStartupMessages(library($pk)); cat(as.character(packageVersion('$pk')))"; done
v R "env/prefix/r_pathway from env/r_pathway.yml" $P/r_pathway/bin/Rscript -e 'cat(R.version.string)'
for pk in clusterProfiler fgsea msigdbr org.Hs.eg.db enrichplot ggplot2; do
  v "R:$pk" "env/prefix/r_pathway" $P/r_pathway/bin/Rscript -e "suppressPackageStartupMessages(library($pk)); cat(as.character(packageVersion('$pk')))"; done
printf "GENCODE annotation\trelease 50 (v50), GRCh38.p14, Ensembl 116 (from GTF header)\tref/ (md5 checked against release MD5SUMS)\n"
} > $OUT
cat $OUT
