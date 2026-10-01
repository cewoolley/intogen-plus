# dNdScv steps of the replication (README.md)
#   Rscript dnds.R regions  <hg19|hg38> <regions.bed>                       coding and splice-site regions of RefCDS
#   Rscript dnds.R annotate <hg19|hg38> <mutations.tsv.gz> <outdir>         annotation of all tumours (groups)
#   Rscript dnds.R strata   <hg19|hg38> <mutations.tsv.gz> <groups.tsv> <outdir>   dNdScv in MSS and MSI tumours
suppressMessages({library(dndscv); library(data.table)})
args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]; build <- args[2]
stopifnot(build %in% c("hg19", "hg38"))

read_mutations <- function(path) {
  m <- fread(cmd = paste("zcat -f", shQuote(path)), colClasses = list(character = c("SAMPLE", "CHROM", "REF", "ALT")))
  unique(as.data.frame(m[, .(sampleID = SAMPLE, chr = sub("^chr", "", CHROM), pos = POS, ref = REF, mut = ALT)]))
}

annotation <- function(out) {
  a <- as.data.table(out$annotmuts)
  a[, .(sampleID, chr, pos, ref, mut, gene, aachange, impact)]
}

if (mode == "regions") {
  e <- new.env(); data(list = ifelse(build == "hg19", "refcds_hg19", "refcds_GRCh38_hg38"), package = "dndscv", envir = e)
  bed <- rbindlist(lapply(e$RefCDS, function(x) {
    cds <- matrix(x$intervals_cds, ncol = 2)
    rbind(data.table(chr = x$chr, start = cds[, 1] - 1, end = cds[, 2]),
          if (length(x$intervals_splice)) data.table(chr = x$chr, start = x$intervals_splice - 1, end = x$intervals_splice))
  }))
  setorder(bed, chr, start)
  fwrite(unique(bed), args[3], sep = "\t", col.names = FALSE)
  cat("regions:", nrow(bed), "\n")
} else if (mode == "annotate") {
  muts <- read_mutations(args[3])
  out <- suppressWarnings(dndscv(muts, refdb = build, max_coding_muts_per_sample = Inf, max_muts_per_gene_per_sample = Inf))
  dir.create(args[4], recursive = TRUE, showWarnings = FALSE)
  fwrite(annotation(out), file.path(args[4], "annotmuts.tsv.gz"), sep = "\t")
} else if (mode == "strata") {
  muts <- read_mutations(args[3])
  groups <- fread(args[4])
  for (s in c("MSS", "MSI")) {
    d <- file.path(args[5], s); dir.create(d, recursive = TRUE, showWarnings = FALSE)
    m <- muts[muts$sampleID %in% groups[GROUP == s]$SAMPLE, ]
    if (length(unique(m$sampleID)) < 20) { cat(s, ": fewer than 20 tumours, skipped\n"); next }
    out <- suppressWarnings(dndscv(m, refdb = build))
    cat(s, "tumours", length(unique(m$sampleID)), "excluded", length(out$exclsamples),
        "drivers q<0.1:", sum(out$sel_cv$qglobal_cv < 0.1), "\n")
    fwrite(out$genemuts, file.path(d, "genemuts.tsv"), sep = "\t")
    fwrite(out$sel_cv, file.path(d, "sel_cv.tsv"), sep = "\t")
    fwrite(annotation(out), file.path(d, "annotmuts.tsv.gz"), sep = "\t")
    fwrite(data.table(SAMPLE = as.character(out$exclsamples)), file.path(d, "excluded.tsv"), sep = "\t")
  }
  data("covariates_hg19_hg38_epigenome_pcawg", package = "dndscv")
  write.table(data.frame(gene = rownames(covs), covs), file.path(args[5], "covariates_pcs.tsv"), sep = "\t", quote = FALSE,
              row.names = FALSE)
} else {
  stop("mode must be regions, annotate or strata")
}
