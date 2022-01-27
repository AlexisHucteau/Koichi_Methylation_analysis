library(dplyr)
library(bsseq)
library(DSS)

setwd("~/GitHub/Koichi_Methylation_analysis/Scripts/")

DATA <- "/media/alexis/DATA/WGBS_IDH/"

# meth.bw = bigwig files of methylation ratios at individual CpGs
# .coverage.bw = bigwig files of CpG read coverage

Meth_ratio <- read.csv(paste0(DATA, "Bed_file/aml111_412761_tcga_2835.meth.bw.wig.bed"), sep = "\t", check.names = F, header = F)
Cov_ratio <- read.csv(paste0(DATA, "Cov_bed_file/aml111_412761_tcga_2835.coverage.bw.wig.bed"), sep = "\t", check.names = F, header = F)

Methylation_count <- data.frame("chr" = Cov_ratio$V1,
                                "pos" = Cov_ratio$V2,
                                "N" = Cov_ratio$V5,
                                "X" = round(Meth_ratio$V5*Cov_ratio$V5),
                                "id" = Cov_ratio$V4)
write.csv(Methylation_count, paste0(DATA, "Merged_files/aml111_412761_tcga_2835.csv"), col.names = F, row.names = T)

