# Real-data inputs for the case study: dNdScv on TCGA MC3 exomes and MSigDB gene sets
suppressMessages({library(dndscv); library(data.table)})
args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]
cohorts <- args[-1]
dir.create(out_dir, showWarnings = FALSE)

# Gene sets: MSigDB 7.5.1 hallmarks and Reactome (as bundled in msigdbr 7.5.1)
gs_file <- file.path(out_dir, "gene_sets.tsv.gz")
if (!file.exists(gs_file)) {
  e <- new.env(); load(Sys.getenv("MSIGDBR_SYSDATA", "msigdbr/R/sysdata.rda"), envir = e)
  sets <- as.data.table(e$msigdbr_genesets)[(gs_cat == "H") | (gs_cat == "C2" & gs_subcat == "CP:REACTOME")]
  links <- as.data.table(e$msigdbr_geneset_genes)[gs_id %in% sets$gs_id]
  genes <- as.data.table(e$msigdbr_genes)[, .(gene_id, SYMBOL = human_gene_symbol)]
  tab <- unique(merge(merge(links, sets[, .(gs_id, gs_name, gs_cat)], by = "gs_id"), genes, by = "gene_id")[
    , .(SET = gs_name, SOURCE = ifelse(gs_cat == "H", "MSigDB_hallmarks", "Reactome"), NAME = gs_name, SYMBOL)])
  fwrite(tab, gs_file, sep = "\t")
  cat("gene sets:", length(unique(tab$SET)), "\n")
}

for (cohort in cohorts) {
  cat("==", cohort, "\n")
  m <- readRDS(file.path(Sys.getenv("TCGAMUTATIONS", "tcgamutations"), "inst", "extdata", "MC3", paste0(cohort, ".RDs")))
  maf <- rbind(attr(m, "data"), attr(m, "maf.silent"), fill = TRUE)
  keep <- c("Missense_Mutation", "Nonsense_Mutation", "Silent", "Splice_Site", "Frame_Shift_Del", "Frame_Shift_Ins",
            "In_Frame_Del", "In_Frame_Ins", "Nonstop_Mutation", "Translation_Start_Site")
  maf <- maf[Variant_Classification %in% keep]
  muts <- unique(data.frame(sampleID = substr(as.character(maf$Tumor_Sample_Barcode), 1, 12),
                            chr = as.character(maf$Chromosome), pos = maf$Start_Position,
                            ref = as.character(maf$Reference_Allele), mut = as.character(maf$Tumor_Seq_Allele2)))
  cat("tumours:", length(unique(muts$sampleID)), " mutations:", nrow(muts), "\n")
  t0 <- Sys.time()
  out <- suppressWarnings(dndscv(muts, outmats = TRUE))
  cat("dNdScv:", format(Sys.time() - t0), " theta:", out$nbreg$theta, " excluded samples:", length(out$exclsamples), "\n")
  d <- file.path(out_dir, cohort); dir.create(d, showWarnings = FALSE)
  fwrite(out$genemuts, file.path(d, "genemuts.tsv"), sep = "\t")
  fwrite(out$sel_cv, file.path(d, "sel_cv.tsv"), sep = "\t")
  fwrite(out$annotmuts[, c("sampleID", "chr", "pos", "ref", "mut", "gene", "impact")], file.path(d, "annotmuts.tsv.gz"), sep = "\t")
  fwrite(out$globaldnds, file.path(d, "globaldnds.tsv"), sep = "\t")
  writeLines(c(sprintf("theta\t%.6g", out$nbreg$theta), sprintf("excluded_samples\t%d", length(out$exclsamples)),
               sprintf("tumours\t%d", length(unique(out$annotmuts$sampleID)))), file.path(d, "info.tsv"))
  saveRDS(list(N = out$N, L = out$L, genes = as.vector(out$genemuts$gene_name)), file.path(d, "mats.rds"))
}

# epigenomic covariates of dNdScv (as exported by the build of the IntOGen datasets)
data("covariates_hg19_hg38_epigenome_pcawg", package = "dndscv")
write.table(data.frame(gene = rownames(covs), covs), file.path(out_dir, "covariates_pcs.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)
data("cancergenes_cgc81", package = "dndscv")
writeLines(known_cancergenes, file.path(out_dir, "cgc81.txt"))
