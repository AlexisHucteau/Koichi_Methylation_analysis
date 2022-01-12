library("dplyr")

setwd("~/GitHub/Koichi_Methylation_analysis/")

SRAtable <- read.csv("WGBS_Wilson.txt")
Primary_Sample_information <- read.csv("Primary_Sample_information.xlsx")

table(SRAtable$Assay.Type)

SRAtable_Usefull <- SRAtable[which(SRAtable[,"Assay.Type"] == "Bisulfite-Seq" & 
                                                    SRAtable[,"body_site"] == "bone_marrow" &  
                                                    SRAtable[,"Instrument"] == "Illumina HiSeq 2500"),]
