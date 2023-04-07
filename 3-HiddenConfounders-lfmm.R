#!/usr/bin/env Rscript
rm(list = ls()); gc()
cat("[AirPol-DNAme] Running Rscript 3-lfmm2...\n")
################################################################################
# Pool EPI+SEP+PEL - Technical Effects using lfmm2
################################################################################
# Author: Lucile
# Date of creation: 13/12/2022
# Check: 06/03/2023
################################################################################
# Effects:
# Model hidden confounders and technical confounders using lfmm
# Run Rscripts 0-Collect*, 1-PreprocessingDNAm and 2-CellComp-Planet before
# Run on a server
#------------------------------------------------------------------------------#
# Notes:
# lfmm: https://bcm-uga.github.io/lfmm/
################################################################################
# Paths and parameters
#------------------------------------------------------------------------------#

analysis <- "Main"

methFile <- "~/AirPollution-DNAme/MetaAnalysis/Data/EPI+SEP+PEL_adjFunnorm_Mval_noOutliers_autosomes.rds"

CCFile <- "~/AirPollution-DNAme/MetaAnalysis/Data/EPIC_CC.planet_rpc_ilr.rds"

metaFile <-  "~/AirPollution-DNAme/PooledAnalysis/Data/Pool_metaData_imp.rds"

expoFile <- "~/AirPollution-DNAme/PooledAnalysis/Data/Pool_envexposures.rds"

exposures <- c("NO2_p", "NO2_t1", "NO2_t2", "NO2_t3", 
               "PM10_p", "PM10_t1",  "PM10_t2", "PM10_t3", 
               "PM25_p", "PM25_t1", "PM25_t2", "PM25_t3")

################################################################################
# R packages
#------------------------------------------------------------------------------#

library(magrittr)

library(lfmm)

################################################################################
# R functions
#------------------------------------------------------------------------------#

# Transform categorical columns into numeric dummy covariate(s)
CatAsDummy <- function(df){
  
  for(i in 1:ncol(data.frame(df))){
    
    x <- fastDummies::dummy_cols(.data = data.frame(df[,i]), 
                                 remove_first_dummy = T, 
                                 remove_selected_columns = T)
    colnames(x) <- paste(colnames(df)[i], 1:ncol(x), sep = "_")
    if(i == 1){
      df.res <- x 
    }else{
      df.res <- cbind.data.frame(df.res, x)
    }
  }
  return( df.res )
}

for(exposure.name in exposures){
  
  saveFile <- paste0("~/AirPollution-DNAme/MetaAnalysis/Data/LF/", analysis, "/PoolEPIC-lfmm2-",exposure.name,".RData")
  
  ################################################################################
  # Load in metadata
  #------------------------------------------------------------------------------#
  
  metaData <- readRDS(file = metaFile)
  
  metaData$CentrexProject <- paste0(metaData$Centre, "x", metaData$Project)
  
  ################################################################################
  # Confounders
  #------------------------------------------------------------------------------#
  
  categorical_confounders <- c("Parity", 
                               "ChildSex",
                               "MaternalEducation", 
                               "Smoke",
                               "CentrexProject")
  
  continuous_confounders <- c("MaternalAge", "MaternalAge2", 
                              "MaternalBMI", "MaternalBMI2", 
                              "GestationalDuration","GestationalDuration2",
                              paste("CC", 1:5, sep = ""))
  
  PHENO <- cbind.data.frame(Sample_Name = metaData$Sample_Name, 
                            id = metaData$id,
                            Cohort = metaData$Cohort,
                            CatAsDummy(metaData[,categorical_confounders]),
                            GestationalDuration = metaData$GestationalDuration,
                            GestationalDuration2 = metaData$GestationalDuration^2,
                            MaternalAge = metaData$MaternalAge,
                            MaternalAge2 = metaData$MaternalAge^2,
                            MaternalBMI = metaData$MaternalBMI,
                            MaternalBMI2 = metaData$MaternalBMI^2)
  
  confounders <- c("Parity_1", 
                   "ChildSex_1", 
                   "MaternalEducation_1",
                   "MaternalEducation_2",
                   "Smoke_1",
                   "CentrexProject_1","CentrexProject_2","CentrexProject_3",
                   continuous_confounders)
  
  PHENO <- merge(PHENO, readRDS(CCFile), by = "Sample_Name")
  
  ################################################################################
  # Exposure
  #------------------------------------------------------------------------------#
  
  Expo <-  readRDS(expoFile)[, c("id", "Cohort", exposure.name)] %>% dplyr::distinct()
  PHENO <- merge(PHENO, Expo, by = c("id", "Cohort"))
  dim(PHENO)
  
  ################################################################################
  # Load methylation data
  #------------------------------------------------------------------------------#
  
  meth <- readRDS(methFile)
  dim(meth)
  #
  
  ################################################################################
  # Filter entries with missing values
  #------------------------------------------------------------------------------#
  
  meth <- na.omit(meth); dim(meth)
  meth <- t(meth)
  meth <- meth[order(rownames(meth)),]
  
  PHENO <- na.omit(PHENO); dim(PHENO)
  PHENO <- PHENO %>% dplyr::arrange(Sample_Name)
  
  meth <- meth[rownames(meth) %in% PHENO$Sample_Name,]
  PHENO <- PHENO[PHENO$Sample_Name %in% rownames(meth),]
  
  #------------------------------------------------------------------------------#
  stopifnot(identical(as.character(PHENO$Sample_Name), rownames(meth)))
  #------------------------------------------------------------------------------#
  
  ################################################################################
  # Estimate the number of latent factors
  #------------------------------------------------------------------------------#
  
  ## Compute residuals from adjusted linear models
  
  formula <- stats::as.formula(paste("Mvalue ~ ", exposure.name, " + ",
                                     paste(confounders, collapse = " + ")))
  print(formula)
  
  meth.resid <- apply(meth, 2,
                      function(y) resid(lm(formula = formula,
                                           data = cbind.data.frame(PHENO, Mvalue = y),
                      )))
  
  ## Compute PCA on lm residuals
  
  pc <- prcomp(meth.resid)
  
  par(mfrow = c(1,1))
  plot(pc$sdev[1:20]^2, xlab = 'PC', ylab = "Variance explained", main = exposure.name)
  points(6,pc$sdev[6]^2, type = "h", lwd = 3, col = "blue")
  
  ################################################################################
  # Save
  #------------------------------------------------------------------------------#
  
  save(pc, formula, file = saveFile)
  
  ################################################################################
  # Clean
  #------------------------------------------------------------------------------#
  
  rm(meth.resid); gc()  
  
  ################################################################################
  # Choose K based on the pca
  #------------------------------------------------------------------------------#
  
  K <- 8
  
  ################################################################################
  # Compute latent factors ajusted for exposure and confounders
  #------------------------------------------------------------------------------#
  
  X <- PHENO[,c(exposure.name,confounders)]
  X <- apply(X, 2, as.numeric)
  rownames(X) <- PHENO$Sample_Name
  
  RUN <- TRUE
  if( RUN ){
    
    coln <- paste("LF", 1:K, sep = "")
    
    lfmm_ridge.results <- lfmm_ridge(Y = meth, 
                                     X = X, 
                                     K = K, 
                                     lambda = 1e-05, 
                                     algorithm = "analytical",
                                     it.max = 100, 
                                     relative.err.min = 1e-06)
    
    rownames(lfmm_ridge.results$U) <- rownames(X)
    colnames(lfmm_ridge.results$U) <- coln
    rownames(lfmm_ridge.results$V) <- colnames(meth)
    colnames(lfmm_ridge.results$V) <- coln
    
  }
  
  ################################################################################
  # Save
  #------------------------------------------------------------------------------#
  
  save(pc, formula, K, lfmm_ridge.results, file = saveFile)
  
}

################################################################################
