library(dplyr)
library(ChAMP)
library(GenomicRanges)

"%ni%" <- Negate("%in%")

setwd("~/GitHub/Koichi_Methylation_analysis/Scripts/")

# samples <- list.files("/media/alexis/DATA/Koichi_methylation_dat/IDAT_files/") %>% stringr::str_split(., pattern = "_") %>% sapply(function(x){x[1:3]})
# Clinical_data <- read.csv("~/GitHub/Koichi_gene_expression_git/Koichi_gene_expression_analyses/DATA/Clinical_patient_data.csv")

# Samples2 <- data.frame("A" = samples[seq(from = 1, to = 633, by = 3)],
#                        "B" = samples[seq(from = 2, to = 633, by = 3)],
#                        "C" = samples[seq(from = 3, to = 633, by = 3)]) %>%
#   unique()
#
# write.csv(Samples2, "/media/alexis/DATA/Koichi_methylation_dat/IDAT_files/Samples2.csv")

pheno <- read.csv("/media/alexis/DATA/Koichi_methylation_dat/samplesheet.csv")

# DATA_loaded <- champ.load("/media/alexis/DATA/Koichi_methylation_dat/IDAT_files/", arraytype = "EPIC")

#
# DATA_loaded$pd <- merge(DATA_loaded$pd, pheno, by.x = "Sample_Name", by.y = "Sample")
#
# champ.QC(beta = DATA_loaded$beta, pheno = DATA_loaded$pd$Pheno)
#
# BMIQ_norm_Koichi_samples <- champ.norm(DATA_loaded$beta, resultsDir = "./BMIQ_Normalization/", arraytype = "EPIC", cores = 8, method = "BMIQ")
# saveRDS(BMIQ_norm_Koichi_samples, "/media/alexis/DATA/Koichi_methylation_dat/BMIQ_norm_Koichi_samples.Rdata")
BMIQ_norm_Koichi_samples <- readRDS("/media/alexis/DATA/Koichi_methylation_dat/BMIQ_norm_Koichi_samples.Rdata")

Control_sample <- colnames(BMIQ_norm_Koichi_samples)[colnames(BMIQ_norm_Koichi_samples) %ni% pheno$Sample]
Control_df <- data.frame("Sample" = Control_sample, "Pheno" = rep("Control", 8))
pheno <- rbind(pheno, Control_df) %>% dplyr::filter(Sample != "") %>% unique()

plot.beta.densities <- function(beta, title) {
  if (!is.null(dim(beta))) {
    densities <- apply(beta, 2, function(x) {
      density(x, na.rm = TRUE)
    })
    xmax <- max(sapply(densities, function(d) {
      max(d$x)
    }))
    xmin <- min(sapply(densities, function(d) {
      min(d$x)
    }))
    ymax <- max(sapply(densities, function(d) {
      max(d$y)
    }))
    plot(NA, xlim = c(xmin, xmax), ylim = c(0, ymax), main = title, ylab = "")
    colors <- rainbow(10)
    for (i in 1:ncol(beta)) {
      lines(densities[[i]], col = colors[i %% 10 + 1])
    }
  } else if (length(beta) > 1) {
    plot(density(beta, na.rm = TRUE), main = title)
  }
}

Baseline_Sample <- pheno[stringr::str_detect(pheno$Pheno, "Baseline"), "Sample"]
Post_response_Sample <- pheno[stringr::str_detect(pheno$Pheno, "Post_treatment"), "Sample"]
Good_responder_Sample <- pheno[stringr::str_detect(pheno$Pheno, "CR") | stringr::str_detect(pheno$Pheno, "CRi"), "Sample"]
Bad_responder_Sample <- pheno[stringr::str_detect(pheno$Pheno, "PD") | stringr::str_detect(pheno$Pheno, "SD"), "Sample"]
IDH1_Sample <- pheno[stringr::str_detect(pheno$Pheno, "IDH1") & !stringr::str_detect(pheno$Pheno, "IDH1_IDH2"), "Sample"]
IDH2_Sample <- pheno[stringr::str_detect(pheno$Pheno, "IDH2") & !stringr::str_detect(pheno$Pheno, "IDH1_IDH2"), "Sample"]
Control_Sample <- pheno[stringr::str_detect(pheno$Pheno, "Control"), "Sample"]

############ GOOD Baseline vs Control

Good_baseline_Control <- c(intersect(Good_responder_Sample, Baseline_Sample), Control_Sample)
Factor_Good_Baseline_Control <- colnames(BMIQ_norm_Koichi_samples)[colnames(BMIQ_norm_Koichi_samples) %in% Good_baseline_Control]
Factor_Good_Baseline_Control <- ifelse(Factor_Good_Baseline_Control %in% Good_responder_Sample, "Good_Responder_Baseline", "Control")

DMR_Good_Baseline_vs_Control <- champ.DMR(BMIQ_norm_Koichi_samples[,colnames(BMIQ_norm_Koichi_samples) %in% Good_baseline_Control],
                                          pheno = Factor_Good_Baseline_Control,
                                          cores = 6,
                                          arraytype = "EPIC")

write.csv(DMR_Good_Baseline_vs_Control$BumphunterDMR, "../Results_DMR/DMR_Good_Baseline_vs_Control.csv")

############ Bad Baseline vs Control

Bad_baseline_Control <- c(intersect(Bad_responder_Sample, Baseline_Sample), Control_Sample)
Factor_Bad_Baseline_Control <- colnames(BMIQ_norm_Koichi_samples)[colnames(BMIQ_norm_Koichi_samples) %in% Bad_baseline_Control]
Factor_Bad_Baseline_Control <- ifelse(Factor_Bad_Baseline_Control %in% Bad_responder_Sample, "Bad_Responder_Baseline", "Control")

DMR_Bad_Baseline_vs_Control <- champ.DMR(BMIQ_norm_Koichi_samples[,colnames(BMIQ_norm_Koichi_samples) %in% Bad_baseline_Control],
                                         pheno = Factor_Bad_Baseline_Control,
                                         cores = 6,
                                         arraytype = "EPIC")

write.csv(DMR_Bad_Baseline_vs_Control$BumphunterDMR, "../Results_DMR/DMR_Bad_Baseline_vs_Control.csv")


############ Response Specific DMR

### Overlapping DMRs

DMR_Good_Baseline_vs_Control_GRanges <- GRanges(
  seqnames = DMR_Good_Baseline_vs_Control$BumphunterDMR$seqnames,
  ranges = IRanges(start = DMR_Good_Baseline_vs_Control$BumphunterDMR$start, end = DMR_Good_Baseline_vs_Control$BumphunterDMR$end),
  DMR_name = rownames(DMR_Good_Baseline_vs_Control$BumphunterDMR),
  chr = DMR_Good_Baseline_vs_Control$BumphunterDMR$seqnames,
  start_chr = DMR_Good_Baseline_vs_Control$BumphunterDMR$start,
  end_chr = DMR_Good_Baseline_vs_Control$BumphunterDMR$end,
  pvalue = DMR_Good_Baseline_vs_Control$BumphunterDMR$p.value
)

DMR_Bad_Baseline_vs_Control_GRanges <- GRanges(
  seqnames = DMR_Bad_Baseline_vs_Control$BumphunterDMR$seqnames,
  ranges = IRanges(start = DMR_Bad_Baseline_vs_Control$BumphunterDMR$start, end = DMR_Bad_Baseline_vs_Control$BumphunterDMR$end),
  DMR_name = rownames(DMR_Bad_Baseline_vs_Control$BumphunterDMR),
  chr = DMR_Bad_Baseline_vs_Control$BumphunterDMR$seqnames,
  start_chr = DMR_Bad_Baseline_vs_Control$BumphunterDMR$start,
  end_chr = DMR_Bad_Baseline_vs_Control$BumphunterDMR$end,
  pvalue = DMR_Bad_Baseline_vs_Control$BumphunterDMR$p.value
)

overlaps_GC_BC <- findOverlaps(DMR_Good_Baseline_vs_Control_GRanges, DMR_Bad_Baseline_vs_Control_GRanges)
overlaps_GC_BC_df <- data.frame(mcols(DMR_Good_Baseline_vs_Control_GRanges[queryHits(overlaps_GC_BC),]), data.frame(mcols(DMR_Bad_Baseline_vs_Control_GRanges[subjectHits(overlaps_GC_BC),])))

Specific_Bad_response <- DMR_Bad_Baseline_vs_Control$BumphunterDMR[rownames(DMR_Bad_Baseline_vs_Control$BumphunterDMR) %ni% overlaps_GC_BC_df$DMR_name.1,]

write.csv(Specific_Bad_response, "../Results_DMR/Specific_Bad_response.csv")

############ GOOD vs BAD RESPONDER AT BASELINE

Good_Bad_Baseline_sample <- intersect(c(Good_responder_Sample, Bad_responder_Sample), Baseline_Sample)
Factor_Good_Bad_Baseline <- colnames(BMIQ_norm_Koichi_samples)[colnames(BMIQ_norm_Koichi_samples) %in% Good_Bad_Baseline_sample]
Factor_Good_Bad_Baseline <- ifelse(Factor_Good_Bad_Baseline %in% Good_responder_Sample, "Good_Responder", "Bad_Responder")

DMR_Good_vs_Bad_Baseline <- champ.DMR(BMIQ_norm_Koichi_samples[,colnames(BMIQ_norm_Koichi_samples) %in% Good_Bad_Baseline_sample],
                                      pheno = Factor_Good_Bad_Baseline,
                                      cores = 6,
                                      arraytype = "EPIC")

write.csv(DMR_Good_vs_Bad_Baseline$BumphunterDMR, "../Results_DMR/DMR_Good_vs_Bad_Baseline.csv")

############ BASELINE vs CD34+

Baseline_Control_sample <- c(Baseline_Sample, Control_Sample)
Factor_Baseline_Control <- colnames(BMIQ_norm_Koichi_samples)[colnames(BMIQ_norm_Koichi_samples) %in% Baseline_Control_sample]
Factor_Baseline_Control <- ifelse(Factor_Baseline_Control %in% Baseline_Sample, "Baseline", "Control")

DMR_Baseline_vs_Control <- champ.DMR(BMIQ_norm_Koichi_samples[,colnames(BMIQ_norm_Koichi_samples) %in% Baseline_Control_sample],
                                     pheno = Factor_Baseline_Control,
                                     cores = 6,
                                     arraytype = "EPIC")

write.csv(DMR_Baseline_vs_Control$BumphunterDMR, "../Results_DMR/DMR_Baseline_vs_Control.csv")

############ GOOD Post vs Control

Good_Post_Control_sample <- c(intersect(Good_responder_Sample, Post_response_Sample), Control_Sample)
Factor_Good_Post_Control <- colnames(BMIQ_norm_Koichi_samples)[colnames(BMIQ_norm_Koichi_samples) %in% Good_Post_Control_sample]
Factor_Good_Post_Control <- ifelse(Factor_Good_Post_Control %in% Good_responder_Sample, "Good_Responder_Post", "Control")

DMR_Good_Post_vs_Control <- champ.DMR(BMIQ_norm_Koichi_samples[,colnames(BMIQ_norm_Koichi_samples) %in% Good_Post_Control_sample],
                                      pheno = Factor_Good_Post_Control,
                                      cores = 6,
                                      arraytype = "EPIC")

write.csv(DMR_Good_Post_vs_Control$BumphunterDMR, "../Results_DMR/DMR_Good_Post_vs_Control.csv")

############ Bad Post vs Control

Bad_Post_Control <- c(intersect(Bad_responder_Sample, Post_response_Sample), Control_Sample)
Factor_Bad_Post_Control <- colnames(BMIQ_norm_Koichi_samples)[colnames(BMIQ_norm_Koichi_samples) %in% Bad_Post_Control]
Factor_Bad_Post_Control <- ifelse(Factor_Bad_Post_Control %in% Bad_responder_Sample, "Bad_Responder_Post", "Control")

DMR_Bad_Post_vs_Control <- champ.DMR(BMIQ_norm_Koichi_samples[,colnames(BMIQ_norm_Koichi_samples) %in% Bad_Post_Control],
                                     pheno = Factor_Bad_Post_Control,
                                     cores = 6,
                                     arraytype = "EPIC")

write.csv(DMR_Bad_Post_vs_Control$BumphunterDMR, "../Results_DMR/DMR_Bad_Post_vs_Control.csv")

############ Response Specific DMR Post treatment

### Overlapping DMRs

DMR_Good_Post_vs_Control_GRanges <- GRanges(
  seqnames = DMR_Good_Post_vs_Control$BumphunterDMR$seqnames,
  ranges = IRanges(start = DMR_Good_Post_vs_Control$BumphunterDMR$start, end = DMR_Good_Post_vs_Control$BumphunterDMR$end),
  DMR_name = rownames(DMR_Good_Post_vs_Control$BumphunterDMR),
  chr = DMR_Good_Post_vs_Control$BumphunterDMR$seqnames,
  start_chr = DMR_Good_Post_vs_Control$BumphunterDMR$start,
  end_chr = DMR_Good_Post_vs_Control$BumphunterDMR$end,
  pvalue = DMR_Good_Post_vs_Control$BumphunterDMR$p.value
)

DMR_Bad_Post_vs_Control_GRanges <- GRanges(
  seqnames = DMR_Bad_Post_vs_Control$BumphunterDMR$seqnames,
  ranges = IRanges(start = DMR_Bad_Post_vs_Control$BumphunterDMR$start, end = DMR_Bad_Post_vs_Control$BumphunterDMR$end),
  DMR_name = rownames(DMR_Bad_Post_vs_Control$BumphunterDMR),
  chr = DMR_Bad_Post_vs_Control$BumphunterDMR$seqnames,
  start_chr = DMR_Bad_Post_vs_Control$BumphunterDMR$start,
  end_chr = DMR_Bad_Post_vs_Control$BumphunterDMR$end,
  pvalue = DMR_Bad_Post_vs_Control$BumphunterDMR$p.value
)

overlaps_GP_BP <- findOverlaps(DMR_Good_Post_vs_Control_GRanges, DMR_Bad_Post_vs_Control_GRanges)
overlaps_GP_BP_df <- data.frame(mcols(DMR_Good_Post_vs_Control_GRanges[queryHits(overlaps_GP_BP),]), data.frame(mcols(DMR_Bad_Post_vs_Control_GRanges[subjectHits(overlaps_GP_BP),])))

Specific_Bad_response_Post <- DMR_Good_Post_vs_Control$BumphunterDMR[rownames(DMR_Good_Post_vs_Control$BumphunterDMR) %ni% overlaps_GP_BP_df$DMR_name.1,]

write.csv(Specific_Bad_response_Post, "../Results_DMR/Specific_Bad_response_Post.csv")

############ GOOD vs BAD RESPONDER AT Post response

Good_Bad_Post_sample <- intersect(c(Good_responder_Sample, Bad_responder_Sample), Post_response_Sample)
Factor_Good_Bad_Post_treatment <- colnames(BMIQ_norm_Koichi_samples)[colnames(BMIQ_norm_Koichi_samples) %in% Good_Bad_Post_sample]
Factor_Good_Bad_Post_treatment <- ifelse(Factor_Good_Bad_Post_treatment %in% Good_responder_Sample, "Good_Responder", "Bad_Responder")

DMR_Good_vs_Bad_Post_treatment <- champ.DMR(BMIQ_norm_Koichi_samples[,colnames(BMIQ_norm_Koichi_samples) %in% Good_Bad_Post_sample],
                                            pheno = Factor_Good_Bad_Post_treatment,
                                            cores = 6,
                                            arraytype = "EPIC")

write.csv(DMR_Good_vs_Bad_Post_treatment$BumphunterDMR, "../Results_DMR/DMR_Good_vs_Bad_Post_treatment.csv")

############ Post treatment vs Baseline

Baseline_Post_sample <- c(Baseline_Sample, Post_response_Sample)
Factor_Baseline_Post_treatment <- colnames(BMIQ_norm_Koichi_samples)[colnames(BMIQ_norm_Koichi_samples) %in% Baseline_Post_sample]
Factor_Baseline_Post_treatment <- ifelse(Factor_Baseline_Post_treatment %in% Baseline_Sample, "Baseline", "Post_response")

DMR_Baseline_vs_Post_treatment <- champ.DMR(BMIQ_norm_Koichi_samples[,colnames(BMIQ_norm_Koichi_samples) %in% Baseline_Post_sample],
                                            pheno = Factor_Baseline_Post_treatment,
                                            cores = 6,
                                            arraytype = "EPIC")

write.csv(DMR_Baseline_vs_Post_treatment$BumphunterDMR, "../Results_DMR/DMR_Baseline_vs_Post_treatment.csv")

############ Wilson DMR WGBS analysis

IDH1_DMR_WGBS <- read.csv("../WGBS_Wilson_Tables/IDH1_DMR_WGBS.csv")
IDH2_DMR_WGBS <- read.csv("../WGBS_Wilson_Tables/IDH2_DMR_WGBS.csv")

DMR_specific_response_GRanges <- GRanges(
  seqnames = Specific_Bad_response$chrom,
  ranges = IRanges(start = Specific_Bad_response$chromStart, end = Specific_Bad_response$chromEnd),
  DMR_number = rownames(Specific_Bad_response)
)

IDH1_DMR_WGBS_GRanges <- GRanges(
  seqnames = IDH1_DMR_WGBS$seqnames,
  ranges = IRanges(start = IDH1_DMR_WGBS$start, end = IDH1_DMR_WGBS$end),
  DMR_number = rownames(IDH1_DMR_WGBS)
)

IDH2_DMR_WGBS_GRanges <- GRanges(
  seqnames = IDH2_DMR_WGBS$seqnames,
  ranges = IRanges(start = IDH2_DMR_WGBS$start, end = IDH2_DMR_WGBS$end),
  DMR_number = rownames(IDH2_DMR_WGBS)
)

A <- IDH1_DMR_WGBS_GRanges
B <- DMR_BR_C

C <- findOverlaps(A, B)
overlaps_Baseline_Response_IDH1_df <- data.frame(mcols(A[queryHits(C),]),
                                                 data.frame(mcols(B[subjectHits(C),])))

A <- DMR_BR_C
B <- IDH2_DMR_WGBS_GRanges

C <- findOverlaps(A, B)
overlaps_Baseline_Response_IDH2_df <- data.frame(mcols(A[queryHits(C),]),
                                                 data.frame(mcols(B[subjectHits(C),])))

Specific_Bad_response_Baseline_IDH1 <- Specific_Bad_response[rownames(Specific_Bad_response) %in% overlaps_GP_BP_df$DMR_name.1,]

############ Change Threshold of DMR analysis

DMR_analysis_threshold <- function(DATA, 
                                   Samples,
                                   Pheno_A_samples, 
                                   Pheno_name_A, 
                                   Pheno_name_B, 
                                   adjPvalDmr_param = 0.05, 
                                   DMR_gap_param = 300, 
                                   minProbes_param = 7, 
                                   file_name_addition = ""){
  Factor <- colnames(DATA)[colnames(DATA) %in% Samples]
  Factor <- ifelse(Factor %in% Pheno_A_samples, "Pheno_name_A", "Pheno_name_B")
  
  DMR <- champ.DMR(DATA[,colnames(DATA) %in% Samples], 
                   adjPvalDmr = adjPvalDmr_param, 
                   maxGap = DMR_gap_param, 
                   minProbes = minProbes_param, 
                   pheno = Factor,
                   cores = 7,
                   arraytype = "EPIC")
  
  write.csv(DMR$BumphunterDMR, paste0("../Results_DMR/DMR_", Pheno_name_A, "_vs_", Pheno_name_B, file_name_addition, ".csv"))
  
  return(DMR)
}


DMR_Bad_Baseline_vs_Control_adjPval_0.1 <- DMR_analysis_threshold(DATA = BMIQ_norm_Koichi_samples, 
                                                                  Samples = Baseline_Control_sample, 
                                                                  Pheno_A_samples = Bad_responder_Sample, 
                                                                  Pheno_name_A = "Bad responder", 
                                                                  Pheno_name_B = "Control", 
                                                                  adjPvalDmr_param = 0.1, 
                                                                  file_name_addition = "_adjPval_0.1")


DMR_Bad_Baseline_vs_Control_gap_3000 <- DMR_analysis_threshold(DATA = BMIQ_norm_Koichi_samples, 
                                                               Samples = Baseline_Control_sample, 
                                                               Pheno_A_samples = Bad_responder_Sample, 
                                                               Pheno_name_A = "Bad responder", 
                                                               Pheno_name_B = "Control",  
                                                               DMR_gap_param = 3000,
                                                               file_name_addition = "_gap_3000")


############ Change overlap parameters

DMR_BR_C <- read.csv("../Results_DMR/DMR_Bad_Baseline_vs_Control.csv")
DMR_BR_C_GRanges <- GRanges(
  seqnames = DMR_BR_C$seqnames,
  ranges = IRanges(start = DMR_BR_C$start, end = DMR_BR_C$end),
  DMR_name = rownames(DMR_BR_C),
  chr = DMR_BR_C$seqnames,
  start_chr = DMR_BR_C$start,
  end_chr = DMR_BR_C$end,
  pvalue = DMR_BR_C$p.value
)

Specific_Bad_response_GRanges <- GRanges(
  seqnames = Specific_Bad_response$seqnames,
  ranges = IRanges(start = Specific_Bad_response$start, end = Specific_Bad_response$end),
  DMR_name = rownames(Specific_Bad_response),
  chr = Specific_Bad_response$seqnames,
  start_chr = Specific_Bad_response$start,
  end_chr = Specific_Bad_response$end,
  pvalue = Specific_Bad_response$p.value
)


A <- IDH1_DMR_WGBS_GRanges
B <- Specific_Bad_response_GRanges

# maxgap	
# A single integer >= -1.
# If type is set to "any", maxgap is interpreted as the maximum gap that is allowed between 2 ranges for the ranges to be considered as overlapping. The gap between 2 ranges is the number of positions that separate them. The gap between 2 adjacent ranges is 0. By convention when one range has its start or end strictly inside the other (i.e. non-disjoint ranges), the gap is considered to be -1.
# If type is set to anything else, maxgap has a special meaning that depends on the particular type. See type below for more information.


for(i in seq(from = 0, to = 1000, by = 50)){
  C <- findOverlaps(A, B, maxgap = i)
  overlaps_ <- data.frame(mcols(A[queryHits(C),]),data.frame(mcols(B[subjectHits(C),])))
  message(paste0("number of overlap for gap: ", i, "\n", nrow(overlaps_)))
}

C <- findOverlaps(A, B, maxgap = 100)
overlaps_ <- data.frame(mcols(A[queryHits(C),]),data.frame(mcols(B[subjectHits(C),])))

A <- DMR_BR_C
B <- IDH2_DMR_WGBS_GRanges

C <- findOverlaps(A, B)
overlaps_Baseline_Response_IDH2_df <- data.frame(mcols(A[queryHits(C),]),
                                                 data.frame(mcols(B[subjectHits(C),])))

Specific_Bad_response_Baseline_IDH1 <- Specific_Bad_response[rownames(Specific_Bad_response) %in% overlaps_GP_BP_df$DMR_name.1,]




############ Associate genes to DMRs







############ Analyse IDH1 & IDH2 separately









############ Investigate DMRs location and compare to their DMR WGBS
