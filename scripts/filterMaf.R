# RNA Editing filtering
library(data.table)

args <- commandArgs(trailingOnly=TRUE)
maf_files <- list.files("/scratch", pattern="norm.vep.coding.maf", recursive=TRUE)

mafs <- lapply(maf_files, fread, stringsAsFactors=FALSE)
maf <- do.call(rbind, mafs)
maf <- subset(maf, BIOTYPE == "protein_coding")
maf <- subset(maf, t_DP >= 10)
maf <- subset(maf, t_GQ >= 20)
maf <- subset(maf, (is.na(gnomAD_AF)) | (gnomAD_AF <= 0.05))
maf <- subset(maf, (is.na(AF)) | (AF <= 0.05))
maf <- subset(maf, Variant_Classification %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del",
						"In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation",
						"Nonstop_Mutation", "Splice_Site"))
maf <- subset(maf, !grepl("benign", PolyPhen))
maf <- subset(maf, !grepl("tolerated", SIFT))
maf$pathogenic <- ifelse(grepl("benign", maf$CLIN_SIG), no=c("pathogenic"), yes=strsplit(maf$CLIN_SIG, ","))
maf <- subset(maf, "pathogenic" %in% pathogenic)
write.table(maf, file="menieres.531.filtered.08222022.tsv", row.names = FALSE, quote = FALSE, sep = "\t")
