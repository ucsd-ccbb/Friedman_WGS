library(data.table)

args <- commandArgs(trailingOnly = TRUE)
maf <- args[1]
outmaf <- gsub("vep", "filt", maf)
maf <- fread(maf, sep="\t", stringsAsFactors=FALSE)
maf$VarID <- paste0(maf$Chromosome, ":", maf$Start_Position, maf$Reference_Allele, ">", maf$Tumor_Seq_Allele2)
maf$VAF <- maf$t_alt_count/maf$t_depth
filt <- subset(maf, Variant_Classification %in% c("Frame_Shift_Ins", "Frame_Shift_Del", "In_Frame_Ins", "In_Frame_Del",
                        "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Region",
                        "Splice_Site", "Translation_Start_Site"))
filt <- subset(filt, (is.na(gnomAD_AF) | gnomAD_AF < 0.01))
filt <- subset(filt, (is.na(ExAC_AF) | ExAC_AF < 0.01))
filt <- subset(filt, (is.na(AF) | AF < 0.01)) # 1000 Genomes
filt <- subset(filt, !grepl("benign", CLIN_SIG))

write.table(filt, file=outmaf, sep="\t", row.names=FALSE, quote=FALSE)
