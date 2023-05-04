#!/usr/bin/env Rscript
cat("[AirP-DNAme] Running Rscript 4-DMRs.R...\n")
################################################################################
# AirP-DNAme -  Run ENMix::comb-p on results from the EWAS
################################################################################
# Author: Lucile
# Date of creation: September 2022
# Note: 
# Run Rscript 4-EWAS
################################################################################
################################################################################
# Paths and parameters
#------------------------------------------------------------------------------#

# Where EWAS results are to be found
resDir <- "~/Work/AirPollution-DNAme/Results/EWAS/Main/DMP/"
saveDir <- "~/Work/AirPollution-DNAme/Results/EWAS/Main/DMR/"
if(!dir.exists(saveDir)) dir.create(saveDir)

RUN.chrX <- FALSE

# Parameters for comb-p
seed <- 0.001
dist <- 500
minNprobes <- 2

ncores <- 3

################################################################################
# R packages
#------------------------------------------------------------------------------#

library(magrittr)

################################################################################
# R functions
#------------------------------------------------------------------------------#

simplifyGeneName <- function(name){ 
  name <- stringr::str_split(name, pattern = ";", simplify = T)
  name <- unique(c(name))
  name <- name[stringi::stri_order(name)]
  name <- paste0(name, collapse = ";")
  return( name )
}

requiredCol <- c("chr", "start", "end", "p", "probe")

runCombp <- function(data,
                     out.file,
                     seed = 0.001,
                     dist = 500,
                     region_plot = FALSE,
                     mht_plot = FALSE,
                     verbose = FALSE,
                     ncores = 1
){
  
  # Check
  if( !identical(intersect(requiredCol, colnames(data)),requiredCol) ) stop("Unvalid input data.frame")
  
  # Format chr names
  data <- data %>% dplyr::mutate(chr = stringr::str_remove(chr, pattern = "chr"))
  
  ENmix::combp(data = data,
               dist.cutoff = dist, 
               bin.size = 310,
               seed = seed,
               region_plot = region_plot, 
               mht_plot = mht_plot, 
               nCores = ncores, 
               verbose = verbose)
  
  #DMR results will be stored in a file with name resu_combp.csv
  system(command =  paste("cd", getwd()))
  system(command =  paste("mv resu_combp.csv", out.file))
  
}

library(GenomicFeatures)

seekNewAnnot <- function(dmrs, 
                         genes, 
                         gene.lab = "geneName", verbose = TRUE){
  
  dmrs.unknown <- dplyr::filter(dmrs, geneName == "")
  
  if( nrow(dmrs.unknown)>0 ){
    dmrs.gr <- GRanges(seqnames = dmrs.unknown$chr, 
                       ranges = IRanges(start = dmrs.unknown$start, end = dmrs.unknown$end),
                       strand = "*")
    
    seqlevelsStyle(genes) <- seqlevelsStyle(dmrs.gr) <- "UCSC"
    
    overlaps <- findOverlaps(query = dmrs.gr, subject = genes, ignore.strand = TRUE)
    
    if(length(overlaps)>0){
      
      if(verbose) cat("New annotation found!! \n") 
      
      gene.id <- genes$gene_id[subjectHits(overlaps)]
      gene.names <- genes$gene_name[subjectHits(overlaps)]
      for(i in seq_along(gene.names)) if(is.na(gene.names[i]) | gene.names[i]=="") gene.names[i] <- gene.id[i]
      
      dmrs.unknown[queryHits(overlaps),"geneName"] <- gene.names

      dmrs <- rbind.data.frame(dplyr::filter(dmrs, geneName != ""), dmrs.unknown)
   }
  }
  return( dmrs )
}
  
getLocations <- function(loc.vector, sep = ";"){
  
  locations <- tapply(loc.vector, seq_along(loc.vector),
                      FUN = function(s) stringr::str_split(s, pattern = sep, simplify = T))
  locations <- unlist(locations)
  locations <- unique(locations)
  locations <- locations[stringi::stri_order(locations)]
  locations <- paste0(locations, collapse = ";")
  
  return(locations)
}

################################################################################
# Load probe annotation
#------------------------------------------------------------------------------#

library( IlluminaHumanMethylationEPICanno.ilm10b4.hg19 )
annotation_object <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19::IlluminaHumanMethylationEPICanno.ilm10b4.hg19

FullAnnot <- minfi::getAnnotation(annotation_object)

FullAnnot <- data.frame(chr = FullAnnot$chr,
                        start = FullAnnot$pos,
                        end = FullAnnot$pos+1,
                        probe = FullAnnot$Name,
                        geneName = FullAnnot$UCSC_RefGene_Name,
                        location = FullAnnot$UCSC_RefGene_Group)
FullAnnot$geneName <- tapply(FullAnnot$geneName,
                              INDEX = seq_along(FullAnnot$geneName),
                              FUN = simplifyGeneName)

#------------------------------------------------------------------------------#
# Recently annotated genes
#------------------------------------------------------------------------------#

# To annotate some unknown regions
#library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#txdb <- makeTxDbFromGFF(file = "~/Data/Annotation/Homo_sapiens.GRCh38.104.gtf.gz",
#                        dbxrefTag = "gene_name")

library(rtracklayer)

genes <- import.gff("~/Data/Annotation/Homo_sapiens.GRCh38.104.gtf.gz")
genes <- genes[genes$type == "gene"]

path <- system.file(package = "liftOver", "extdata", "hg38ToHg19.over.chain")
chain <- import.chain(path)
seqlevelsStyle(genes) <- "UCSC"
genes <- liftOver(genes, chain) %>% unlist()

################################################################################
# RUN for all chromosomes
#------------------------------------------------------------------------------#

inFiles  <- list.files(path = resDir, pattern = ".rds", recursive = T, full.names = T)
outFiles <- stringr::str_remove_all(basename(inFiles), pattern = ".rds")

setwd(saveDir)

for(f in seq_along( inFiles )){
  
  print( inFiles[f] )
  res <- readRDS( inFiles[f] ) 
  res$probe <- rownames(res)

  colnames(res)[grep(pattern = "raw_p_value", x = colnames(res))] <- "raw_pvalue"
  
  #bc <- bacon::bacon(effectsizes = res$Estimate, standarderrors = res$SE)
  #res$pval_corrected <- bacon::pval(bc)
  
  res <- merge(res, FullAnnot, by = "probe")
  
  res <- res %>% dplyr::arrange(raw_pvalue)
  res <- res %>% dplyr::mutate(chr = stringr::str_remove(chr, pattern = "chr"))

  cat("Run comb-p...")
  outFile <- paste0(outFiles[f], ".combp.csv")
    
  bed <- res %>% dplyr::rename(p = raw_pvalue)
  #bed <- res %>% dplyr::rename(p = pval_corrected)
  
  bed <- bed[, c("chr", "start", "end", "p", "probe")]
  bed <- bed %>% dplyr::mutate(chr = stringr::str_remove(chr, pattern = "chr"))
    
  runCombp(data = bed, 
           out.file = outFile, 
           seed = seed, 
           dist = dist)
  cat("OK.\n")
    
  if( file.exists(outFile) ){
    
      bed <- read.csv(file = outFile)
      bed <- bed %>% dplyr::arrange(sidak)
      bed <- bed %>% dplyr::filter(nprobe>=minNprobes)
      
      if( nrow(bed)>0 ){
      
        bed$geneName <- NA
        bed$location <- NA
        bed$Estimate.mean <- NA 
        bed$MeanBeta.mean <- NA
        bed$pvalue.min <- NA
  
        for(d in seq_along(bed$start)){
          x <- res %>% dplyr::filter(chr == bed$chr[d] & start>=bed$start[d] & end<=bed$end[d])
          x <- x %>% dplyr::summarise(geneName = simplifyGeneName(geneName),
                                      location = getLocations(location),
                                      Estimate.mean = mean(Estimate),
                                      MeanBeta.mean = mean(MeanBeta),
                                      pvalue.min = min(raw_pvalue, na.rm = T)) %>%
                     data.frame()
          bed[d,colnames(x)] <- x
        }

        print(bed)
        bed <- seekNewAnnot(bed, genes = genes)
        print( bed )
        
        write.csv(bed, file = outFile, row.names = F)
        
        }else file.remove(outFile)
  }
}

################################################################################
# RUN for chromosome X
#------------------------------------------------------------------------------#

if(RUN.chrX){
  
  inFiles  <- list.files(path = resDir, pattern = ".rds", recursive = T, full.names = T)
  outFiles <- stringr::str_remove_all(basename(inFiles), pattern = "chrX.rds")
  
  setwd(saveDir)
  
  for(f in seq_along( inFiles )){
    
    print( inFiles[f] )
    res <- readRDS( inFiles[f] ) 
    res$probe <- rownames(res)
    
    colnames(res)[grep(pattern = "raw_p_value", x = colnames(res))] <- "raw_pvalue"
    
    #bc <- bacon::bacon(effectsizes = res$Estimate, standarderrors = res$SE)
    #res$pval_corrected <- bacon::pval(bc)
    
    res <- merge(res, FullAnnot, by = "probe")
    
    res <- res %>% dplyr::arrange(raw_pvalue)
    res <- res %>% dplyr::mutate(chr = stringr::str_remove(chr, pattern = "chr"))
    
    cat("Run comb-p...")
    outFile <- paste0(outFiles[f], ".combp.csv")
    
    bed <- res %>% dplyr::rename(p = raw_pvalue)
    #bed <- res %>% dplyr::rename(p = pval_corrected)
    
    bed <- bed[, c("chr", "start", "end", "p", "probe")]
    bed <- bed %>% dplyr::mutate(chr = stringr::str_remove(chr, pattern = "chr"))
    bed <- bed %>% dplyr::filter(chr == "X")
    
    runCombp(data = bed, 
             out.file = outFile, 
             seed = seed, 
             dist = dist)
    cat("OK.\n")
    
    if( file.exists(outFile) ){
      
      bed <- read.csv(file = outFile)
      bed <- bed %>% dplyr::arrange(sidak)
      bed <- bed %>% dplyr::filter(nprobe>=minNprobes)
      
      if( nrow(bed)>0 ){
        
        bed$geneName <- NA
        bed$location <- NA
        bed$Estimate.mean <- NA 
        bed$MeanBeta.mean <- NA
        bed$pvalue.min <- NA
        
        for(d in seq_along(bed$start)){
          x <- res %>% dplyr::filter(chr == bed$chr[d] & start>=bed$start[d] & end<=bed$end[d])
          x <- x %>% dplyr::summarise(geneName = simplifyGeneName(geneName),
                                      location = getLocations(location),
                                      Estimate.mean = mean(Estimate),
                                      MeanBeta.mean = mean(MeanBeta),
                                      pvalue.min = min(raw_pvalue, na.rm = T)) %>%
            data.frame()
          bed[d,colnames(x)] <- x
        }
        
        print(bed)
        bed <- seekNewAnnot(bed, genes = genes)
        print( bed )
        
        write.csv(bed, file = outFile, row.names = F)
        
      }else file.remove(outFile)
    }
  }
}

################################################################################
rm(list = ls())
################################################################################
