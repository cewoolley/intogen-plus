# TCGA MC3 COAD + READ coding mutations (one tumour per patient)
suppressMessages(library(data.table))
keep <- c("Missense_Mutation", "Nonsense_Mutation", "Silent", "Splice_Site", "Frame_Shift_Del", "Frame_Shift_Ins",
          "In_Frame_Del", "In_Frame_Ins", "Nonstop_Mutation", "Translation_Start_Site")
maf <- rbindlist(lapply(c("COAD", "READ"), function(c) {
  m <- readRDS(file.path(Sys.getenv("TCGAMUTATIONS", "tcgamutations"), "inst", "extdata", "MC3", paste0(c, ".RDs")))
  d <- rbind(attr(m, "data"), attr(m, "maf.silent"), fill = TRUE)
  d[Variant_Classification %in% keep, .(SAMPLE = as.character(Tumor_Sample_Barcode), PROJECT = c,
       chr = as.character(Chromosome), pos = Start_Position, ref = as.character(Reference_Allele),
       mut = as.character(Tumor_Seq_Allele2), gene = as.character(Hugo_Symbol),
       class = as.character(Variant_Classification), type = as.character(Variant_Type), protein = as.character(HGVSp_Short))]
}))
maf[, PATIENT := substr(SAMPLE, 1, 12)]
# one tumour sample per patient (the first primary tumour aliquot)
first <- maf[, .(SAMPLE = sort(unique(SAMPLE))[1]), by = PATIENT]
maf <- unique(maf[first, on = c("PATIENT", "SAMPLE")])
cat("patients:", length(unique(maf$PATIENT)), " mutations:", nrow(maf), "\n")
fwrite(maf, "coadread_mutations.tsv.gz", sep = "\t")

out <- rbindlist(lapply(c("COAD", "READ"), function(c) {
  m <- readRDS(file.path(Sys.getenv("TCGAMUTATIONS", "tcgamutations"), "inst", "extdata", "MC3", paste0(c, ".RDs")))
  cl <- as.data.table(attr(m, "clinical.data"))
  idcol <- intersect(c("Tumor_Sample_Barcode", "bcr_patient_barcode", "patient_id"), names(cl))[1]
  cl[, .(PATIENT = substr(get(idcol), 1, 12), PROJECT = c, SITE = anatomic_neoplasm_subdivision,
         MMR_IHC_LOSS = loss_expression_of_mismatch_repair_proteins_by_ihc, HISTOLOGY = histological_type,
         STAGE = pathologic_stage, AGE = age_at_initial_pathologic_diagnosis, SEX = gender)]
}))
fwrite(out, "clinical.tsv", sep = "\t")
print(table(out$PROJECT)); print(table(out$SITE, useNA = "ifany")); print(table(out$MMR_IHC_LOSS, out$PROJECT, useNA = "ifany")); print(table(out$HISTOLOGY))
