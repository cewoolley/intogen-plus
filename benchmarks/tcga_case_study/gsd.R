# dNdScv genesetdnds (Martincorena et al. 2017) for many gene sets.
# Same model and fitting code as dndscv::genesetdnds; the matrices of the rest of the
# exome are computed as total - set instead of summing all the other genes each time.
# Usage: Rscript gsd.R <mats.rds> <units.tsv> <output.tsv> [validate]
suppressMessages({library(dndscv); library(data.table)})
args <- commandArgs(trailingOnly = TRUE)
mats <- readRDS(args[1]); units <- fread(args[2]); out_file <- args[3]
data("submod_192r_3w", package = "dndscv")
syneqs <- substmodel[, 1]
allg <- mats$genes
Ltot <- rowSums(mats$L, dims = 2); Ntot <- rowSums(mats$N, dims = 2)

fit_substmodel <- function(N, L, substmodel, testpar) {   # verbatim from dndscv::genesetdnds
  l = c(L); n = c(N); r = c(substmodel)
  n = n[l!=0]; r = r[l!=0]; l = l[l!=0]
  params = unique(base::strsplit(x=paste(r,collapse="*"), split="\\*")[[1]])
  indmat = as.data.frame(array(0, dim=c(length(r),length(params))))
  colnames(indmat) = params
  for (j in 1:length(r)) { indmat[j, base::strsplit(r[j], split="\\*")[[1]]] = 1 }
  model = glm(formula = n ~ offset(log(l)) + . -1, data=indmat, family=poisson(link=log))
  mle = exp(coefficients(model))
  model.lrt = drop1(model, test="LRT", scope=testpar)
  pvals.lrt = setNames(model.lrt[[5]], row.names(model.lrt))
  data.frame(name=names(mle), mle=mle, pval.lrt=pvals.lrt[names(mle)])
}

geneset <- function(genes) {
  gi <- which(allg %in% genes)
  Ls <- rowSums(mats$L[, , gi, drop = FALSE], dims = 2); Ns <- rowSums(mats$N[, , gi, drop = FALSE], dims = 2)
  L <- rbind(Ls, Ltot - Ls); N <- rbind(Ns, Ntot - Ns)
  rm2 <- array("", dim = dim(L))
  rm2[,1] = c(paste(syneqs,"*r_rel",sep=""), syneqs)
  rm2[,2] = c(paste(syneqs,"*r_rel*wmis_geneset",sep=""), paste(syneqs,"*wmis_rest",sep=""))
  rm2[,3] = c(paste(syneqs,"*r_rel*wtru_geneset",sep=""), paste(syneqs,"*wtru_rest",sep=""))
  rm2[,4] = c(paste(syneqs,"*r_rel*wtru_geneset",sep=""), paste(syneqs,"*wtru_rest",sep=""))
  p2 <- fit_substmodel(N, L, rm2, testpar = c("wmis_geneset", "wtru_geneset"))
  rm3 <- array("", dim = dim(L))
  rm3[,1] = c(paste(syneqs,"*r_rel",sep=""), syneqs)
  for (k in 2:4) rm3[,k] = c(paste(syneqs,"*r_rel*wall_geneset",sep=""), paste(syneqs,"*wall_rest",sep=""))
  p3 <- fit_substmodel(N, L, rm3, testpar = c("wall_geneset"))
  p <- rbind(p2, p3); p <- p[grepl("_geneset", p$name), ]; rownames(p) <- p$name
  c(n_genes = length(gi),
    wall = p["wall_geneset", "mle"], p_wall = p["wall_geneset", "pval.lrt"],
    wmis = p["wmis_geneset", "mle"], p_wmis = p["wmis_geneset", "pval.lrt"],
    wtru = p["wtru_geneset", "mle"], p_wtru = p["wtru_geneset", "pval.lrt"])
}

if (length(args) > 3 && args[4] == "validate") {
  # compare with the original function on the first units
  data("submod_192r_3w", package = "dndscv")
  fake <- list(N = mats$N, L = mats$L, genemuts = data.frame(gene = allg))
  for (i in 1:3) {
    genes <- intersect(strsplit(units$GENES[i], ";")[[1]], allg)
    orig <- genesetdnds(fake, genes)$globaldnds_geneset
    mine <- geneset(genes)
    cat(units$UNIT[i], " original wall:", orig["wall_geneset", "mle"], orig["wall_geneset", "pval.lrt"],
        " fast:", mine["wall"], mine["p_wall"], "\n")
  }
  quit(save = "no")
}

t0 <- Sys.time()
res <- rbindlist(lapply(seq_len(nrow(units)), function(i) {
  genes <- strsplit(units$GENES[i], ";")[[1]]
  as.data.table(c(list(UNIT = units$UNIT[i]), as.list(geneset(genes))))
}))
fwrite(res, out_file, sep = "\t")
cat("genesetdnds:", nrow(units), "units in", format(Sys.time() - t0), "\n")
