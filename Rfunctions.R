################################################################################
# Author: Lucile
# Date: July 2021 
################################################################################
# Functions for AirPollution-DNAm EWAS
################################################################################

################################################################################
#R packages
#------------------------------------------------------------------------------#

library(MASS) # rlm function for robust linear regression
library(sandwich) #Huber's estimation of the standard error
library(lmtest) # to use coeftest
library(stringr)

################################################################################
# Annotation of CpG sites 
#------------------------------------------------------------------------------#
#source("https://bioconductor.org/biocLite.R")
#biocLite("IlluminaHumanMethylation450k.db")

#library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
#data("IlluminaHumanMethylation450kanno.ilmn12.hg19")

#annotation.table <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

################################################################################
#Function removeOutliers()
#------------------------------------------------------------------------------#

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
  N_for_probe <- rowSums(!is.na(probes))
  Log <- data.frame(initial_NAs,removed_lower,removed_upper,N_for_probe)
  
  return( list(probes, Log) )
}

################################################################################
# Turn beta to Mvalue
#------------------------------------------------------------------------------#

logit2 <- function(p) ifelse(p*(1-p)!=0, return( log2(p/(1-p)) ), return( NA )) 

################################################################################
# Transform M-value to beta-value
#------------------------------------------------------------------------------#

m2beta <- function(m){return ( 2^m/(2^m + 1) ) }

################################################################################
# Function Lu_rlm_pollution()
# Performs robust linear regression based on M-values + intercept method 
#------------------------------------------------------------------------------#

# Data supposed to be M-values
Lu_rlm_pollution <- function(meth, y, 
                             main_exposure,
                             xconf.names, 
                             cellTypes = 0, 
                             method = "M"){ 
  
  #if( method %nin% c("M", "MM") ) stop("method argument must be either 'M' or 'MM', see documentation for MASS:rlm() \n")
  
  CpG <- meth[y,]
  median <- median(CpG, na.rm = T)
  IQR <- IQR(CpG, na.rm = T)
  
  ret <- cbind(median = median, IQR = IQR, MeanBeta = NA, Estimate = NA, SE = NA, zvalue = NA, raw_pvalue = NA)
  names(ret) <- c("median", "IQR", "MeanBeta", "Estimate","SE","zvalue","raw_pvalue")
  
  mod = try(
    rlm(formula = stats::as.formula(paste("CpG ~", main_exposure,
                                          "+",
                                          paste(xconf.names, collapse = "+"), 
                                          "+",
                                          paste(cellTypes, collapse = "+")
    )
    ), 
    maxit=400, method = method)
  )
  
  if( !("try-error" %in% class(mod)) ){
    cf <- try(
      coeftest(mod, vcov = vcovHC(mod, type = "HC0"))  
    )
    if( !("try-error" %in% class(cf)) ){
      ret <- c(median = median, IQR = IQR)
      #if( transformToMvalue ){
      Intercept <- mod$coefficients["(Intercept)"]
      Estimate <- mod$coefficients[main_exposure]
      MeanBeta <- m2beta(Intercept+Estimate) - m2beta(Intercept)
      ret <- c(ret, MeanBeta = MeanBeta)
      #}
      ret <- c(ret, cf[main_exposure, c(1,2,3,4)])
      names(ret) <- c("median", "IQR", "MeanBeta", "Estimate","SE","zvalue","raw_pvalue")
    }
  }
  
  return(ret)
  
}
################################################################################
