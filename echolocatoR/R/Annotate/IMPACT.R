
# IMPACT #
# https://github.com/immunogenomics/IMPACT


IMPACT.get_annotations <- function(chrom=NULL, 
                                subset_DT=NULL,
                                nThread=4){ 
  if(!is.null(subset_DT)){
    chrom <- subset_DT$CHR[1]
  }
  baseURL <- "https://github.com/immunogenomics/IMPACT/raw/master/IMPACT707/Annotations"
  URL <- file.path(baseURL, paste0("IMPACT707_EAS_chr",chrom,".annot.gz"))
  annot <- data.table::fread(URL, nThread = nThread)
  
  if(!is.null(subset_DT)){ 
    annot_merge <- data.table:::merge.data.table(data.table::data.table(subset_DT),
                                                 annot,
                                                 by.x = c("SNP","CHR","POS"), 
                                                 by.y = c("SNP","CHR","BP"),
                                                 all = F)  
    return(annot_merge)
  } else {
    return(annot)
  } 
}


IMPACT.get_ldscores <- function(chrom=NULL, 
                                subset_DT=NULL,
                                nThread=4){ 
  warning("LDSCores do not include any SNPs with MAF<0.5%, as they are restricted to HapMap3 SNPs. \
This may affect subsequent analyss (e.g. fine-mapping).")
  if(!is.null(subset_DT)){
    chrom <- subset_DT$CHR[1]
  }
  baseURL <- "https://github.com/immunogenomics/IMPACT/raw/master/IMPACT707/LDscores"
  URL <- file.path(baseURL, paste0("IMPACT707_EAS_chr",chrom,".l2.ldscore.gz"))
  ldscore <- data.table::fread(URL, nThread = nThread)
  
  if(!is.null(subset_DT)){ 
    ldscore_merge <- data.table:::merge.data.table(data.table::data.table(subset_DT),
                                  ldscore,
                                  by.x = c("SNP","CHR","POS"), 
                                  by.y = c("SNP","CHR","BP"),
                                  all = F)  
    return(ldscore_merge)
  } else {
    return(ldscore)
  } 
}

 