#!/usr/bin/env Rscript
rm(list = ls()); gc()
################################################################################
# Pooled-DNAm: Main analysis
################################################################################
# Author: Lucile
# Date of creation: 13/12/2022
################################################################################
# Pre-processing features:
#
# - Normalisation: InterpolatedXY (Noob + Functional norm.)
# - Exposure assessment method: Hough's model
#
# Modelling:
#
# Site-wise robust Linear Models
# M-value ~ Exposure + 
#           ChildSex +
#           Smoke + 
#           Parity + 
#           MaternalEducation +
#           sp(GestationalDuration) +
#           sp(MaternalAge)  +
#           sp(MaternalBMI) + 
#           CentrexProject +
#           TechnicalEffects (Latent Factors) + 
#           RefCellTypes.ilr
#
################################################################################
# Notes:
# - Using newly curated EDEN database (Summer 2022; data provided by Emie)
################################################################################
# Paths and parameters
#------------------------------------------------------------------------------#

workDir <-  "~/AirPollution-DNAme/"

analysis <- "Main"

rm_samples <- c()

#------------------------------------------------------------------------------#
# Data (standardized)
#------------------------------------------------------------------------------#

# Methylation data
methFile <- "~/AirPollution-DNAme/MetaAnalysis/Data/EPI+SEP+PEL_adjFunnorm_Mval_noOutliers_autosomes.rds"
transformToMvalue <- FALSE

# Exposure data
expoFile <- "~/AirPollution-DNAme/PooledAnalysis/Data/Pool_envexposures.rds"

# Metadata
metaFile <- "~/AirPollution-DNAme/PooledAnalysis/Data/Pool_metaData_imp.rds"

# Cell-type composition (reference-based method)
CCFile <- "~/AirPollution-DNAme/MetaAnalysis/Data/EPIC_CC.planet_rpc_ilr.rds"
CCFile.prop <- "~/AirPollution-DNAme/MetaAnalysis/Data/EPIC_CC.planet_rpc.rds"

exposures <- c("NO2", "PM10", "PM25")
windows <- c("p", "t2", "t1", "t3")

################################################################################
# Set workspace
#------------------------------------------------------------------------------#

setwd(workDir)

resDir <- "Results/EWAS/"
if(!dir.exists(resDir)) dir.create(resDir)

resDir <- paste0(resDir, "/", analysis)
if(!dir.exists(resDir)) dir.create(resDir)

resDir <- paste0(resDir, "/DMP")
if(!dir.exists(resDir)) dir.create(resDir)

################################################################################
# R packages
#------------------------------------------------------------------------------#

library(magrittr)
library(rms)

#------------------------------------------------------------------------------#
# Source functions (copied from file)
#------------------------------------------------------------------------------#

source("~/AirPollution-DNAme/PooledAnalysis/Rscripts/Rfunctions.R")

################################################################################
# Association test : general parameters (to be set)
#------------------------------------------------------------------------------#

continuous_confounders <- c("MaternalAge_sp",
                            "GestationalDuration_sp",
                            "MaternalBMI_sp")

categorical_confounders <- c("ChildSex", 
                             "Parity", 
                             "MaternalEducation",
                             "Smoke",
                             "CentrexProject")

K <- 5 # how many latent factors?
hidden_confounders <- paste0("LF",1:K)

cellTypes <- paste0("CC", 1:5)  # set cellTypes <- 0 for not adjusting

mc.cores <- 3

################################################################################
#------------------------------------------------------------------------------#
################################################################################
# Load meta and exposure data
#------------------------------------------------------------------------------#

for(exposure in exposures){
  for(window in windows){
    
    # Latent factors (lfmm)
    LFFile <- paste0("~/AirPollution-DNAme/MetaAnalysis/Data/LF/Main/PoolEPIC-lfmm2-", exposure, "_", window,".RData")
    
    PHENO <- readRDS( metaFile )
    PHENO$CentrexProject <- paste0(PHENO$Centre, "x", PHENO$Project)
    dim(PHENO)
    # 1539 18
    
    PHENO[,categorical_confounders] <- apply(PHENO[,categorical_confounders], 2, as.factor)
    
    expo <- readRDS(expoFile) %>% dplyr::distinct()
    PHENO <- merge(PHENO, expo, by = c("Sample_Name", "Cohort", "Project"), all.y = F)
    PHENO <- PHENO %>% dplyr::arrange(Sample_Name)
    dim(PHENO)
    # 1539   39
    
    if( length(rm_samples)>0 ) PHENO <- PHENO %>% dplyr::filter(!(Sample_Name %in% rm_samples))
    
    # Reference Cell-type composition
    CC <- readRDS(CCFile)
    PHENO <- merge(PHENO, CC, by = "Sample_Name")
    dim(PHENO)
    # [1] 1539   44
    
    # Latent factors
    load(LFFile)
    
    length( rownames(lfmm_ridge.results$U) )
    
    PHENO <- merge(PHENO, 
                   data.frame(Sample_Name = rownames(lfmm_ridge.results$U),
                              lfmm_ridge.results$U[,hidden_confounders]),
                   by = "Sample_Name")
    dim(PHENO)
    #[1] 1538   52
    
    rm(lfmm_ridge.results, pc)
    
    ################################################################################
    # Filter out samples for which Syncytio < 25%
    #------------------------------------------------------------------------------#
    
    CC <- readRDS(CCFile.prop)
    
    rmSamples.cc <- CC$Sample_Name[which(CC[,"Syncytiotrophoblast"]<0.25)]
    
    if(length(rmSamples.cc)>0){
      cat(length(rmSamples.cc), "samples removed due to spurious cell-type composition estimates. \n")
      print(rmSamples.cc)
      PHENO <- PHENO %>% dplyr::filter(!(Sample_Name %in% rmSamples.cc))
    }
    
    dim(PHENO)
    # 391  23
    
    rm(CC) 
    
    #------------------------------------------------------------------------------#
    attach( PHENO )
    #------------------------------------------------------------------------------#
    
    ################################################################################
    # Load methylation data
    #------------------------------------------------------------------------------#
    
    meth <- readRDS(methFile)
    
    if( nrow(meth)<ncol(meth) ) meth <- t( meth )
    meth <- meth[,order(colnames(meth))]
    dim(meth)
    # [1] 346524   1539
    
    ################################################################################
    # Synchronize data sets
    #------------------------------------------------------------------------------#
    
    PHENO <- PHENO %>% dplyr::filter(Sample_Name %in% colnames(meth)) %>% dplyr::arrange(Sample_Name)
    meth <- meth[, colnames(meth) %in% PHENO$Sample_Name]
    
    stopifnot(identical(as.character(PHENO$Sample_Name), colnames(meth)))
    
    ################################################################################
    # Model continuous covariates as restricted cubic splines
    # - MaternalAge
    # - MaternalBMI
    # - GestationalDuration
    #------------------------------------------------------------------------------#
    
    MaternalAge_sp <- rms::rcs( PHENO[,"MaternalAge"] ) 
    
    MaternalBMI_sp <- rms::rcs( PHENO[,"MaternalBMI"] ) 
    
    GestationalDuration_sp <- rms::rcs( PHENO[,"GestationalDuration"] ) 
    
    ################################################################################
    # Run exome-wide association tests
    #------------------------------------------------------------------------------#
    
    confounders <- c(continuous_confounders, categorical_confounders, hidden_confounders)
    
    main_exposure <- paste0(exposure, "_", window)
    
    outFile <- paste0(resDir, "/", main_exposure, ".rds")
    print( Sys.time() )
    cat("Output file:", outFile, "\n")
    
    system.time({
      
      results <- pbmcapply::pbmcmapply(
        FUN = Lu_rlm_pollution, 
        y = rownames(meth),
        MoreArgs = list(meth = meth, 
                        main_exposure = main_exposure, 
                        xconf.names = confounders,
                        cellTypes = cellTypes,
                        method = "M"
        ),
        mc.cores = mc.cores
      )
      
    })
    
    results <- data.frame( t( results ) )
    results$window <- window
    results$exposure <- exposure
    
    saveRDS(results, file = outFile) 
    cat("[", Sys.time(), "]", "saved. \n")
    
  }
}

################################################################################
cat("Done. \n")
################################################################################
