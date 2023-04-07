#!/usr/bin/env Rscript
rm(list = ls()); gc()
cat("[AirPol-DNAme] Running Rscript 0-DNAme-Preprocessing...\n")
################################################################################
# Pooled-DNAme (EPIC) - Preprocessing methylation data (pooled data set)
################################################################################
# Author: Lucile
# Date of creation: 04/11/2022
# Check: 06/03/2023
################################################################################
# Steps:
# 1. Normalisation: InterpolatedXY adjusted funnorm (already performed)
# 2. Removal of CpGs close to known SNPs
# 3. Removal of cross-reactive CpGs 
# 4. Transformation to M-values
# --> Save meth data set for QC and technical var estimation 
# 5. Removal of outlying probes (cf: Abraham2018)
# 6. Removal of technical replicates
# 7. Split autosomes (F+M) and autosomes+chrX (F) & autosomes + chrXY (M)
#------------------------------------------------------------------------------#
# Notes:
################################################################################
# Paths and parameters
#------------------------------------------------------------------------------#

# Normalized data 
methFile <- "~/Data/Pooled-DNAm/METHYL/EPI+SEP+PEL_adjFunnorm.rds"

designFile <- "~/Data/Pooled-DNAm/METHYL/FDF+EPI+SEP+PEL_design.csv"

# Metadata, with missing value imputation (cf Rscript: 0-CollectMetaData)
metaFile <- "~/Work/AirPollution-DNAme/PooledAnalysis/Data/Pool_metaData_imp.rds"

# Probe annotations
annotFile.EPIC <- "~/Work/Misc/Data/IlluminaHumanMethylationEPICanno+Homo_sapiens.GRCh38.104.rds"

# Data with duplicates (eg: to compute CpG dispersion, technical variation)
outFile0 <- "~/Data/Pooled-DNAm//METHYL/EPI+SEP+PEL_adjFunnorm_Mval.rds"  # used for QC
# Final data set, ready for EWAS on autosomes
outFile1 <- "~/Data/Pooled-DNAm//METHYL/EPI+SEP+PEL_adjFunnorm_Mval_noOutliers_autosomes.rds"
# Final data set, ready for sex-stratified EWAS on all chromosomes
outFile.F <- "~/Data/Pooled-DNAm//METHYL/EPI+SEP+PEL_adjFunnorm_Mval_noOutliers_females.rds"
outFile.M <- "~/Data/Pooled-DNAm//METHYL/EPI+SEP+PEL_adjFunnorm_Mval_noOutliers_males.rds"

################################################################################
# R packages
#------------------------------------------------------------------------------#

library(Gmisc, quietly = TRUE)
library(glue)
library(htmlTable)
library(grid)
library(magrittr)

################################################################################
# R functions
#------------------------------------------------------------------------------#

# Exact same as the function used for paper (Abraham2018)
removeOutliers <- function(probes){
  
  require(matrixStats)
  
  if(nrow(probes) < ncol(probes)) warning("expecting probes are rows (long dataset)")
  
  rowIQR <- rowIQRs(probes, na.rm = T)
  row2575 <- rowQuantiles(probes, probs = c(0.25, 0.75), na.rm = T)
  maskL <- probes < row2575[,1] - 3 * rowIQR 
  maskU <- probes > row2575[,2] + 3 * rowIQR 
  initial_NAs <- rowSums(is.na(probes))
  probes[maskL] <- NA
  removed_lower <- rowSums(is.na(probes))-initial_NAs
  probes[maskU] <- NA
  removed_upper <- rowSums(is.na(probes))-removed_lower-initial_NAs
  N_for_probe <-rowSums(!is.na(probes))
  Log <- data.frame(initial_NAs,removed_lower,removed_upper,N_for_probe)
  
  return( list(probes, Log) )
}

################################################################################
# 1. Get normalized beta values (custom) and design
#------------------------------------------------------------------------------#

meth <- readRDS(methFile)
dim(meth)
# [1] 803343    928

if(nrow(meth) < ncol(meth)) meth <- t( meth )

Nprobes <- c(norm = nrow(meth))

targets <- read.csv(designFile)
head(targets); dim(targets)
# [1] 1686    5
setdiff(targets$Sample_Id, colnames(meth))
# [1] "9588557124_R01C01"
# -> sample removed by QC
targets <- targets %>% dplyr::filter(Sample_Id %in% colnames(meth))
targets$idnew <- paste0(stringr::str_sub(targets$Dataset, end = 3), "_", targets$Sample_Name)

# Match samples
meth <- meth[,targets$Sample_Id]

stopifnot(identical(colnames(meth), targets$Sample_Id))

colnames(meth) <- targets$idnew

length(unique(colnames(meth))) == ncol(meth) # Check ids are unique, they are not.

################################################################################
# 2. Removal of CpGs close to SNPs (dist<2 bp)
#------------------------------------------------------------------------------#

meth <- DMRcate::rmSNPandCH(meth, 
                            dist = 2, 
                            mafcut = 0.05,
                            rmcrosshyb = FALSE, 
                            rmXY = FALSE)
dim(meth)
# 384686   1685

Nprobes <- c(Nprobes, snps = nrow(meth))

################################################################################
# 3. Removal of known cross-reactive probes
#------------------------------------------------------------------------------#

xloci <- maxprobes::xreactive_probes(array_type = "450K")
length(xloci)
# 38941 known cross-reactive probes in 450K arrays

x <- which(rownames(meth) %in% xloci)
length(x)
# 32905 known cross-reactive probes still in the dataset

meth <- meth[-x,]
dim(meth)
# 758415    928 

xloci <- maxprobes::xreactive_probes(array_type = "EPIC") %>% unlist()
length(xloci)
# 87464 known cross-reactive probes in EPIC arrays

x <- which(rownames(meth) %in% xloci)
length(x)
# 14716 known cross-reactive probes still in the dataset
if(length(x)>0){
  meth <- meth[-x,]
  dim(meth)
}

Nprobes <- c(Nprobes, crossh = nrow(meth))

################################################################################
# 4. Transformation into M-values
# Slit due to memory shortage
#------------------------------------------------------------------------------#

meth <- cbind(minfi::logit2(meth[,1:200]),
              minfi::logit2(meth[,201:400]),
              minfi::logit2(meth[,401:600]),
              minfi::logit2(meth[,601:ncol(meth)]))

stopifnot(identical(colnames(meth), targets$idnew))

dim(meth)
#[1] 743699    928

################################################################################
cat("Pre-processed M-values with duplicated samples saved in file:", outFile0, "\n")
#------------------------------------------------------------------------------#

saveRDS(meth, file = outFile0)

################################################################################
# 5. Removal of technical duplicates 
#------------------------------------------------------------------------------#

set.seed(38) # Seed for random picking (EPIMEX, SEPAGES, PELAGIE)

# From EDEN-EPIMEX
duplicates <- list(c("4280", "4280_1","4280_2","4280_3","4280_4","4280_5", "4280_6"),
                   c("3779", "3779_1","3779_2","3779_3","3779_4","3779_5"),
                   c("7645_1b", "7645_2", "7645_3b", "7645_4"),
                   c("2812_1_3", "2812_1", "2812_2","2812_3"),  
                   c("2789_1", "2789_2"))

# Keep replicates (random picking)
keptSamples <- c()
## Iterate on duplicates list
for(i in seq_along(duplicates)){
  x <- targets %>% dplyr::filter(idnew %in% paste0("EPI_", duplicates[[i]]))
  x <- sample(x = x$Sample_Name, size = 1)
  keptSamples <- c(keptSamples,x)
}
rm(x)
keptSamples

duplicates <- unlist(duplicates)

targets <- targets %>% dplyr::filter(!(idnew %in% setdiff(paste0("EPI_", duplicates),paste0("EPI_", keptSamples))))
table(targets$Dataset)
#EPIMEX PELAGIE SEPAGES 
#382     96     432

# SEPAGES and PELAGIE: duplicated ids
duplicates <- targets %>% dplyr::group_by(idnew) %>% dplyr::summarise(count = dplyr::n()) %>% data.frame()
duplicates <- duplicates %>% dplyr::filter(count>1)
duplicates <- duplicates$idnew

## Iterate on duplicates list
for(i in seq_along(duplicates)){
  x <- which(targets$idnew == duplicates[i])
  x <- sample(x = x, size = length(x)-1)
  targets <- targets[-x,]
}
table(targets$Dataset)
#EPIMEX   PELAGIE SEPAGES 
# 382       94     395 

meth <- meth[,targets$idnew]
dim(meth)
# [1] 743699    871

# RM replicate number
colnames(meth) <- tapply(colnames(meth), seq_along(colnames(meth)), 
                         FUN = function(x){
                           sp <- stringr::str_split(x, pattern = "_", simplify = T)
                           sp <- paste(sp[1:2], collapse = "_")
                           return(sp)
                         })

################################################################################
# 6. Removal of outlying probes
# Split due to memory shortage
#------------------------------------------------------------------------------#

x1 <- removeOutliers(meth[1:100000,]); x1 <- x1[[1]]
x2 <- removeOutliers(meth[100001:200000,]); x2 <- x2[[1]] 
x <- rbind(x1,x2); rm(x1,x2)
x3 <- removeOutliers(meth[200001:300000,]); x3 <- x3[[1]] 
x <- rbind(x,x3); rm(x3)
x4 <- removeOutliers(meth[300001:400000,]); x4 <- x4[[1]] 
x <- rbind(x,x4); rm(x4)
x5 <- removeOutliers(meth[400001:500000,]); x5 <- x5[[1]] 
x <- rbind(x,x5); rm(x5)
x6 <- removeOutliers(meth[500001:600000,]); x6 <- x6[[1]] 
x <- rbind(x,x6); rm(x6)
x7 <- removeOutliers(meth[600001:nrow(meth),]); x7 <- x7[[1]] 

meth <- rbind(x,x7); rm(x7, x)
dim(meth)
#[1] 743699    871

#Filtrage pour exclusion des CpGs qui ont plus de 25% de données manquantes
nbdmcpg <- apply(meth,1, function(x) sum(is.na(x))) #nb missing per CpG
nbdmcpg <- data.frame(nbdmcpg)
exclucpg <- nbdmcpg %>% dplyr::filter(nbdmcpg>=25*ncol(meth)/100) #liste des Cpg exclus car plus de 25% de données manquantes
rownames(exclucpg)

if( nrow(exclucpg)>0 ) meth <- meth[, which(!(rownames(meth)%in%exclucpg))]
dim(meth)
# 
# 0 CpGs were removed as none contained >25% of missing values
sum(nbdmcpg$nbdmcpg)/(ncol(meth)*nrow(meth))*100
# 0.226% of the values in the data set were outliers. 

Nprobes <- c(Nprobes, outliers = nrow(meth))

################################################################################
# 7. Split data set according to child sex
#------------------------------------------------------------------------------#

metaData <- readRDS(metaFile)

# recode metadat$id
ids.M <- metaData$Sample_Name[which(metaData$ChildSex == "Male")]
meth.M <- meth[,colnames(meth) %in% ids.M]
dim(meth.M)
# [1] 743699    470
saveRDS(meth.M, file = outFile.M)
rm(meth.M, ids.M)

ids.F <- metaData$Sample_Name[which(metaData$ChildSex == "Female")]
meth.F <- meth[,colnames(meth) %in% ids.F]
dim(meth.F)
#  [1] 743699    401

################################################################################
# 7. Removal of Chr X and Chr Y (F data set, joint data set - M+F)
#------------------------------------------------------------------------------#

annot <- data.frame(readRDS(annotFile.EPIC), array_type = "EPIC")
annot <- annot %>% dplyr::distinct(CpG, chr, pos, .keep_all = TRUE)

probesX <- annot$CpG[annot$chr == "chrX"]; length(probesX)
#[1] 19090
probesY <- annot$CpG[annot$chr == "chrY"]; length(probesY)
# [1] 537

meth.F <- meth.F[!(rownames(meth.F) %in% probesY),]
dim(meth.F)
# [1] 743668    401
saveRDS(meth.F, file = outFile.F)
rm(meth.F, ids.F)

meth <- meth[!(rownames(meth) %in% c(probesX, probesY)),]
dim(meth)
# 728335    871

Nprobes <- c(Nprobes, chrXY = nrow(meth))

################################################################################
cat("Pre-processed M-values after outlier exclusion saved in file:", outFile1, "\n")
#------------------------------------------------------------------------------#

saveRDS(meth, file = outFile1)

################################################################################
rm(list = ls())
################################################################################
