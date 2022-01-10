library(dplyr)
library(ChAMP)

setwd("~/GitHub/Koichi_Methylation_analysis/Scripts/")

# samples <- list.files("/media/alexis/DATA/Koichi_methylation_dat/IDAT_files/") %>% stringr::str_split(., pattern = "_") %>% sapply(function(x){x[1:3]})
# Clinical_data <- read.csv("~/GitHub/Koichi_gene_expression_git/Koichi_gene_expression_analyses/DATA/Clinical_patient_data.csv")

# Samples2 <- data.frame("A" = samples[seq(from = 1, to = 633, by = 3)],
#                        "B" = samples[seq(from = 2, to = 633, by = 3)],
#                        "C" = samples[seq(from = 3, to = 633, by = 3)]) %>%
#   unique()
# 
# write.csv(Samples2, "/media/alexis/DATA/Koichi_methylation_dat/IDAT_files/Samples2.csv")

DATA_loaded <- champ.load("/media/alexis/DATA/Koichi_methylation_dat/IDAT_files/", arraytype = "EPIC")
