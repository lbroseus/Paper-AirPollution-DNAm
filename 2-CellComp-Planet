#!/usr/bin/env Rscript
rm(list = ls()); gc()
cat("[AirP-DNAme] Running Rscript 1-Pool-CellComp-Planet...\n")
################################################################################
# Pool EPIC - Estimation of reference cell type composition using planet
# Zero imputation pooling all cohorts
################################################################################
# Author: Lucile
# Date of Creation: 04/11/2022
# Check: 06/03/2023
################################################################################
# Procedure:
# 1. Import Yuan's reference as in Planet
# 2. Filter probes not in the processed methylation dataset
# 3. Use planet and the RPC method to estimate CC
# 4. Impute zeros
# 5. Apply ilr transformation
#------------------------------------------------------------------------------#
# Notes:
################################################################################
# Paths and parameters
#------------------------------------------------------------------------------#

methFile <- "~/Data/Pooled-DNAm/METHYL/EPI+SEP+PEL_adjFunnorm_Mval_noOutliers_autosomes.rds"

saveFile.prop <- "~/Work/AirPollution-DNAme/Data/EPIC_CC.planet_rpc.rds"
saveFile.ilr <- "~/Work/AirPollution-DNAme/Data/EPIC_CC.planet_rpc_ilr.rds"

################################################################################
# R packages
#------------------------------------------------------------------------------#

library(magrittr)
library(planet)
library(EpiDISH)

################################################################################
# R function
#------------------------------------------------------------------------------#

applyPlanet <- function(methFile, plCellCpGsThird, method = "RPC", betaValues = FALSE){
  
  cat(methFile, "\n")
  betaValues <- readRDS( methFile )
  
  if(!betaValues) betaValues <- minfi::ilogit2(betaValues)
  if( ncol(betaValues) > nrow(betaValues)) betaValues <- t( betaValues )
  
  betaValues <- betaValues[order(rownames(betaValues)),]
  
  nbCpG <- length(intersect(rownames(betaValues), rownames(plCellCpGsThird)))
  # X remaining probes
  cat("Common reference probes:", nbCpG, "\n")

  CC.rpc <- epidish(beta.m = betaValues, 
                    ref.m = plCellCpGsThird, 
                    method = method) 
  
  CC.rpc <- CC.rpc$estF
  rownames(CC.rpc) <- colnames(betaValues)
  
  return( CC.rpc )
  
}

################################################################################
#------------------------------------------------------------------------------#
################################################################################
# 1. Processing Yuan's reference probes
#------------------------------------------------------------------------------#

data("plCellCpGsThird")

plCellCpGsThird <- plCellCpGsThird[order(rownames(plCellCpGsThird)),]
plCellCpGsThird <- plCellCpGsThird[,order(colnames(plCellCpGsThird))]

# Number of reference probes:
nbCpG <- nrow(plCellCpGsThird)
nbCpG
# 600

refCellTypes <- c("Endothelial", 
                  "Hofbauer", "nRBC", 
                  "Stromal", 
                  "Syncytiotrophoblast", "Trophoblasts")

################################################################################
# 2. Restrict to probes available in the dataset (after QC filtering)
# 3. Estimate CC using the RPC projection method
#------------------------------------------------------------------------------#

CCprop <- applyPlanet(methFile = methFile, plCellCpGsThird = plCellCpGsThird, betaValues = FALSE)
# Common reference probes: 510 
CCprop <- data.frame(Sample_Name = rownames(CCprop), CCprop)
rownames(CCprop) <- NULL

summary(CCprop)

################################################################################
# 4. Impute zero values in compositions
#------------------------------------------------------------------------------#

X <- CCprop[,refCellTypes]
X <- apply(X, 1:2, function(x) ifelse(x==0, NA, x))
X <- robCompositions::impCoda(x = X)$xImp
CCprop[,refCellTypes] <- X

apply(X,2,summary)

rm(X)

saveRDS(file = saveFile.prop, CCprop)

################################################################################
# Ilr transformation
# Save reference and estimates
#------------------------------------------------------------------------------#

CC.ilr <- compositions::ilr(CCprop[,refCellTypes])
CC.ilr <- apply(CC.ilr, 2, as.numeric)
colnames(CC.ilr) <- paste0("CC", 1:ncol(CC.ilr))

CC.ilr <- cbind.data.frame(Sample_Name = CCprop[, "Sample_Name"], CC.ilr)

saveRDS(file = saveFile.ilr, CC.ilr)

################################################################################
rm(list = ls())
################################################################################
