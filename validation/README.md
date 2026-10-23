# validation/ (not part of QUANT or ANALYSIS)
Paper-comparison code moved out of the product path in the refactor (design section 8). Nothing here feeds back into any cutoff.
- `bin/check_vs_paper.R <deseq2_results_all.tsv> <ora_dir> <gsea_dir> <padj>`: compares our DE genes / pathways with the expected lists of Hosseini et al. 2023 (Hosseini-specific gene and pathway lists inside the script).
- `check_vs_paper.nf`: the old Nextflow process wrapper (no longer included by nextflow/main.nf).
Run by hand, e.g. (r_dge env): `env/prefix/r_dge/bin/Rscript validation/bin/check_vs_paper.R results/06_dge/deseq2_results_all.tsv results/07_pathway/ora results/07_pathway/gsea 0.05`
Backups of the originals: modules/check_vs_paper.nf.bak_pre_refactor, nextflow/bin/check_vs_paper.R.bak_pre_refactor.
