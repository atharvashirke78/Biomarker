# shared helper: gene sets from msigdbr (human), tolerant to the old (category/subcategory) and new (collection/subcollection) API
suppressPackageStartupMessages(library(msigdbr))
get_sets <- function(coll, sub = NULL) {
  r <- tryCatch({ if (is.null(sub)) msigdbr(species = "Homo sapiens", collection = coll) else msigdbr(species = "Homo sapiens", collection = coll, subcollection = sub) },
                error = function(e) { if (is.null(sub)) msigdbr(species = "Homo sapiens", category = coll) else msigdbr(species = "Homo sapiens", category = coll, subcategory = sub) })
  r
}
kegg_sets <- function() {   # C2 CP:KEGG (legacy KEGG set names such as KEGG_PPAR_SIGNALING_PATHWAY)
  cols <- tryCatch(msigdbr_collections(), error = function(e) NULL)
  sub <- "CP:KEGG_LEGACY"
  if (!is.null(cols)) { k <- unique(cols[[grep("sub", names(cols))[1]]]); k <- k[grepl("^CP:KEGG", k)]; if (length(k)) sub <- if ("CP:KEGG_LEGACY" %in% k) "CP:KEGG_LEGACY" else k[1] }
  tryCatch(get_sets("C2", sub), error = function(e) get_sets("C2", "CP:KEGG"))
}
sets_to_list <- function(df, col = "gene_symbol") split(df[[col]], df$gs_name)
