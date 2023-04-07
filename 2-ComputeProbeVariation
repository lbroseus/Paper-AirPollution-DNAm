#!/usr/bin/env Rscript
rm(list = ls()); gc()
cat("[AirP-DNAme] Running Rscript X.R...\n")
################################################################################
# Pool EPIC - DNAme-technical variation 
################################################################################
# Author: Lucile
# Date of creation: 20/12/2022
# Note: computes intra-,inter-individual and batch variances
################################################################################
################################################################################
# Paths and parameters
#------------------------------------------------------------------------------#

#methFile <-  "~/AirPollution-DNAme/MetaAnalysis/Data/EPI+SEP+PEL_adjFunnorm_Mval.rds"
methFile <-  "~/AirPollution-DNAme/PooledAnalysis/Data/FDF+EPI+SEP+PEL_adjFunnorm_Mval.rds"

#outFile <- "~/AirPollution-DNAme/MetaAnalysis/Results/EPI+SEP+PEL_ProbeVariation_PRadj.rds"
outFile <- "~/AirPollution-DNAme/PooledAnalysis/Results/FDF+EPI+SEP+PEL_ProbeVariation_PRadj.rds"

metaFile <- "~/AirPollution-DNAme/PooledAnalysis/Data/Pool_metaData_imp.rds"

pooled.experiment <- TRUE

################################################################################
# R packages
#------------------------------------------------------------------------------#

library(magrittr)

################################################################################
# R functions
#------------------------------------------------------------------------------#

countReplicates <- function(ids, pooled.experiment){
  
  ids <- data.frame(id = ids)
  if(pooled.experiment){
    ids$sample_id <- tapply(ids$id, seq_along(ids$id), 
                            function(s) paste(stringr::str_split(string = s, pattern = "\\_", simplify = T)[c(1,2)], 
                                              collapse = "_"))
  }else{
    ids$sample_id <- tapply(ids$id, seq_along(ids$id), 
                            function(s) paste(stringr::str_split(string = s, pattern = "\\.", simplify = T)[c(1)]))
  }
  
  count <- table(ids$sample_id)
  ids <- merge(ids, data.frame(sample_id = names(count), Freq = as.vector(count)), by = "sample_id")
  
  return(ids)
}

estimateProbeVar <- function(meth, ids, grouping = NULL, verbose = TRUE){
  
  # Subset samples with at least 2 technical replicates
  ids <- ids[ids$Freq>=2,]
  ids <- ids %>% dplyr::arrange(id)
  meth <- meth[,order(colnames(meth))]
  
  index <- which(colnames(meth) %in% ids$id)
  meth <- meth[,index]
  if(!is.null(grouping)) grouping <- as.factor(grouping[index])
  colnames(meth) <- ids$sample_id
  ids <- ids %>% dplyr::distinct(sample_id, .keep_all = TRUE)
  if( verbose ) cat("There are", nrow(ids), "placental samples for which DNAm was measured at least two times.\n")
  
  gc()
  
  indiv <- as.factor(colnames(meth))
  
  #VCA::varPlot(CpG ~ indiv, Data = Data)
  if(is.null(grouping)){
    formula <- CpG ~ indiv
    Data <- data.frame(indiv = indiv,  CpG = NA)
    probeVar <- data.frame(CpG = rownames(meth), total = NA, indiv = NA, error = NA)
  }else{
    formula <- CpG ~ indiv + grouping
    Data <- data.frame(indiv = indiv, grouping = grouping, CpG = NA)
    probeVar <- data.frame(CpG = rownames(meth), total = NA, indiv = NA, grouping = NA, error = NA)
  }
  
  for(i in seq_along(rownames(meth))){
    Data$CpG <- meth[i,]
    vca.res <- try(VCA::fitVCA(form = formula, Data = Data, "reml"))
    if( "try-error" %in% class(vca.res) ){
      probeVar[i,2:ncol(probeVar)] <- c(NA,NA)
    }else{
      probeVar[i,2:ncol(probeVar)] <- as.vector(vca.res$aov.tab[,"VC"])
    }
  }
  
  return( probeVar )
}

################################################################################
################################################################################
# 1. Load methylation data
#------------------------------------------------------------------------------#

meth <- readRDS(methFile)
meth <- meth[,order(colnames(meth))]
cat("DNAm data set with technical replicates contains", ncol(meth), "samples. \n")

################################################################################
# 1'. Define grouping (experiment = project)
#------------------------------------------------------------------------------#

grouping <- tapply(colnames(meth), seq_along(colnames(meth)),
                   FUN = function(x) stringr::str_split(x, pattern = "_", simplify = T)[1])

################################################################################
# 2. 
# Detect replicated samples using ids
# Compute technical sd
#------------------------------------------------------------------------------#

ids <- countReplicates(ids = colnames(meth), pooled.experiment = pooled.experiment)

probeVar <- estimateProbeVar(ids, meth = meth, grouping = grouping)

################################################################################
# Save
#------------------------------------------------------------------------------#

saveRDS(probeVar, file = outFile)

################################################################################
cat("Done.\n")
################################################################################
