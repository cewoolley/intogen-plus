suppressMessages({library(dndscv); library(data.table)})
m <- fread(cmd = "zcat coadread_mutations.tsv.gz")
groups <- fread("tumour_groups.tsv")
sets <- list(ALL = groups$PATIENT, MSS = groups[GROUP == "MSS"]$PATIENT, MSI = groups[GROUP == "MSI"]$PATIENT)
for (s in names(sets)) {
  d <- file.path("data", s); dir.create(d, recursive = TRUE, showWarnings = FALSE)
  muts <- unique(as.data.frame(m[PATIENT %in% sets[[s]], .(sampleID = PATIENT, chr, pos, ref, mut)]))
  t0 <- Sys.time()
  out <- suppressWarnings(dndscv(muts, outmats = TRUE))
  cat(s, "tumours", length(unique(muts$sampleID)), "excluded", length(out$exclsamples), "theta", out$nbreg$theta,
      "drivers q<0.1:", sum(out$sel_cv$qglobal_cv < 0.1), format(Sys.time() - t0), "\n")
  fwrite(out$genemuts, file.path(d, "genemuts.tsv"), sep = "\t")
  fwrite(out$sel_cv, file.path(d, "sel_cv.tsv"), sep = "\t")
  fwrite(out$annotmuts[, c("sampleID", "chr", "pos", "ref", "mut", "gene", "impact")], file.path(d, "annotmuts.tsv.gz"), sep = "\t")
  fwrite(out$globaldnds, file.path(d, "globaldnds.tsv"), sep = "\t")
  if (!is.null(out$sel_loc)) fwrite(out$sel_loc, file.path(d, "sel_loc.tsv"), sep = "\t")
  saveRDS(list(N = out$N, L = out$L, genes = as.vector(out$genemuts$gene_name)), file.path(d, "mats.rds"))
}

# gene sets (MSigDB 7.5.1 hallmarks and Reactome, from msigdbr 7.5.1) and the epigenomic covariates of dNdScv
e <- new.env(); load(Sys.getenv("MSIGDBR_SYSDATA", "msigdbr/R/sysdata.rda"), envir = e)
sets <- as.data.table(e$msigdbr_genesets)[(gs_cat == "H") | (gs_cat == "C2" & gs_subcat == "CP:REACTOME")]
links <- as.data.table(e$msigdbr_geneset_genes)[gs_id %in% sets$gs_id]
genes <- as.data.table(e$msigdbr_genes)[, .(gene_id, SYMBOL = human_gene_symbol)]
tab <- unique(merge(merge(links, sets[, .(gs_id, gs_name, gs_cat)], by = "gs_id"), genes, by = "gene_id")[
  , .(SET = gs_name, SOURCE = ifelse(gs_cat == "H", "MSigDB_hallmarks", "Reactome"), NAME = gs_name, SYMBOL)])
fwrite(tab, "data/gene_sets.tsv", sep = "\t")
data("covariates_hg19_hg38_epigenome_pcawg", package = "dndscv")
write.table(data.frame(gene = rownames(covs), covs), "data/covariates_pcs.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
