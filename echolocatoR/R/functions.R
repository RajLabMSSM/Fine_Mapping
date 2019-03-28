#  ((((((((((((((((()))))))))))))))))  # 
 # ((((((((((((((((())))))))))))))))) # 
 #    ------- echolocatoR -------    #
 # ((((((((((((((((())))))))))))))))) # 
#              ((((()))))              # 

# Author: Brian M. Schilder
  # Bioinformatician II
  # Icahn School of Medicine at Mount Sinai
  # New York City, New York, USA
  # https://bschilder.github.io/BMSchilder

     # =/\                  /\=
    #  / \'._   (\_/)   _.'/ \
   #  / .''._'--(o.o)--'_.''. \
  #  /.' _/ |`'=/ " \='`| \_ `.\
 #  /` .' `\;-,'\___/',-;/` '. '\
#  /.-'       `\(-V-)/`       `-.\
# `              "   "           `


# You can learn more about package authoring with RStudio at:
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'
 
# Load libraries
.libPaths()

library(readxl)
library(DT)
library(data.table)
library(dplyr)
library(ggplot2)
library(plotly)
library(cowplot)
library(ggrepel)
library(curl)
library(biomaRt)
# library(sqldf)
# Ensembl LD API
library(httr)
library(jsonlite)
library(xml2)
library(gaston)
library(RCurl)
library(tidyr)
library(biomaRt)

# *** susieR ****
# library(knitrBootstrap) #install_github('jimhester/knitrBootstrap')
library(susieR) # devtools::install_github("stephenslab/susieR")

# *** finemapr ****
## finemapr contains: finemap, CAVIAR, and PAINTOR
# library(finemapr) # devtools::install_github("variani/finemapr")
# library(roxygen2) #roxygenize()

# *** locuscomparer ****
# https://github.com/boxiangliu/locuscomparer
# library(locuscomparer); #devtools::install_github("boxiangliu/locuscomparer")

# thm <- knitr::knit_theme$get("bipolar")
# knitr::knit_theme$set(thm)



## General Functions
createDT <- function(DF, caption="", scrollY=400){
  data <- DT::datatable(DF, caption=caption,
                        extensions = 'Buttons',
                        options = list( dom = 'Bfrtip',
                                        buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                                        scrollY = scrollY, scrollX=T, scrollCollapse = T, paging = F,
                                        columnDefs = list(list(className = 'dt-center', targets = "_all"))
                        )
  )
  return(data)
}

createDT_html <- function(DF, caption="", scrollY=400){
  htmltools::tagList( createDT(DF, caption, scrollY))
} 


get_dataset_name <- function(file_path){
  dataset_name <- strsplit(dirname(file_path), "/")[[1]][3]
  return(dataset_name)
}


Data_dirs_to_table <- function(Data_dirs, writeCSV=F){
  df <- data.frame()
  for(n in names(Data_dirs)){
    newRow <- data.frame(Data_dirs[n], stringsAsFactors = F)
    colnames(newRow) <- c("Type","topSumStats", "fullSumStats","Reference")
    df <- rbind(df, cbind(Dataset=n,newRow))
  }  
  createDT(df)
  if(writeCSV!=F){ 
    data.table::fwrite(df, writeCSV, quote = F, sep = ",", row.names = F) 
  }
  return(df)
}

# Data Preprocessing
biomart_snps_to_geneInfo <- function(snp_list){
  # listMarts()
  snp_mart = useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp")
  # View(listFilters(snp_mart))
  # View(listAttributes(snp_mart))
  snp_results <- getBM(snp_mart, filters="snp_filter", values=snp_list,
                       attributes=c("refsnp_id","snp","chr_name", "chrom_start","chrom_end",
                                    "associated_gene","ensembl_gene_stable_id" ) )
  # # Split ensembl IDs
  gene_mart = useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
  gene_results <- getBM(mart = gene_mart, filters = "ensembl_gene_id",
                        # values = unlist(strsplit(snp_results$ensembl, ";")),
                        values = snp_results$ensembl_gene_stable_id,
                        attributes = c("hgnc_symbol","external_gene_name","ensembl_gene_id",
                                       "chromosome_name", "start_position", "end_position") )
  snp_results <-snp_results %>%
    mutate(ensembl = strsplit(as.character(ensembl_gene_stable_id), ";")) %>%
    tidyr::unnest(ensembl)
  merged_df <- data.table(gene_results, key = "ensembl_gene_id")[data.table(snp_results, key = "ensembl")]
  return(merged_df)
}
# biomart_snps_to_geneInfo(c("rs114360492"))

biomart_geneInfo <- function(geneList){
  # listDatasets(useMart("ENSEMBL_MART_ENSEMBL") )
  gene_mart = useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
  # View(listFilters(gene_mart))
  # View(listAttributes(gene_mart))
  gene_results <- getBM(mart = gene_mart, filters = "hgnc_symbol",
                        # values = unlist(strsplit(snp_results$ensembl, ";")),
                        values = geneList,
                        attributes = c("hgnc_symbol","external_gene_name","ensembl_gene_id",
                                      "chromosome_name", "start_position", "end_position") )
  return(gene_results)
}
# biomart_geneInfo(c("PTK2B","CLU","APOE"))



## TSS Data
# Use bioMart to get TSS positions for each gene
get_TSS_position <- function(gene){
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  # att <- listAttributes(mart)
  # grep("transcription", att$name, value=TRUE)
  TSS <- getBM(mart=mart,
               attributes=c("hgnc_symbol","transcription_start_site", "version"),
               filters="hgnc_symbol", values=gene)
}



## Import Sig GWAS/QTL Summary Statistics
# For each gene, get the position of top SNP (the one with the greatest effect size/Beta)
import_topSNPs <- function(file_path, caption="", sheet = 1,
                           chrom_col="CHR", position_col="POS", snp_col="SNP",
                           pval_col="P", effect_col="Effect", gene_col="Gene"){
  ## Only the significant subset of results
  if(endsWith(file_path, ".xlsx") | endsWith(file_path, ".xlsm")){
    top_SNPs <- read_excel(path = file_path, sheet = sheet)
  } else if (endsWith(file_path, ".csv")){
    top_SNPs <- fread(file=file_path, sep = ",", header = T, stringsAsFactors = F )
  }
  else if (endsWith(file_path, ".txt")){
    top_SNPs <- fread(file=file_path, sep = "\t", header = T, stringsAsFactors = F )
  } else {print("File type must be .xlsx, .cs, or tab-delimited .txt")}
  
  top_SNPs <- subset(top_SNPs, select=c(chrom_col, position_col, snp_col,
                                        pval_col, effect_col, gene_col ))
  # Standardize names
  colnames(top_SNPs) <- c("CHR","POS","SNP","P","Effect","Gene")
  # Get only the top SNP (sorted by lowest p-val, then highest Effect size) for each gene
  top_SNPs <- top_SNPs %>% arrange(P, desc(Effect)) %>% group_by(Gene) %>% slice(1)
  top_SNPs$CHR <- gsub("chr", "",top_SNPs$CHR)
  top_SNPs <- cbind(Coord= paste(top_SNPs$CHR, top_SNPs$POS, sep=":"),
                    top_SNPs)
  createDT(top_SNPs, caption)
  return(top_SNPs)
}
# top_SNPs <- import_topSNPs(
#   file_path = Data_dirs$Nalls_2019$topSS,
#   sheet="Data",
#   chrom_col = "CHR", position_col = "BP", snp_col="SNP",
#   pval_col="P, all studies", effect_col="Beta, all studies", gene_col="Nearest Gene",
#   caption= "Nalls et al. (2018) PD GWAS Summary Stats")


# import info from FUMA instead

import_FUMA <- function(topSS_path, geneList){
  # topSS_path = Data_dirs$Kunkle_2019$topSS
  # risk_loci <- fread(file.path(dirname(topSS_path), "FUMA/GenomicRiskLoci.txt"))
  # mapped_genes <- fread(file.path(dirname(topSS_path), "FUMA/magma.genes.out"))
  annovar <- fread(file.path(dirname(topSS_path), "FUMA/annov.txt")) 
  annovar_sub <- subset(annovar, symbol %in% geneList)
  
  mapped_genes_sub <- subset(mapped_genes, SYMBOL %in% geneList)
  
  # Get gene range directly from biomart, then query for snps in that range
  geneSNPs <- function(gene, file_path,
                       chrom_col="CHR", position_col="POS", snp_col="SNP",
                       pval_col="P", effect_col="Effect", stderr_col="StdErr", sep="\t"){ 
    # bm <- biomart_geneInfo(gene)
    gene_start <- min(annovar_sub$pos)
    gene_end <- max(annovar_sub$pos)
    gene_chr <- unique(annovar_sub$chr)[1]
    
    colDict <- column_dictionary(file_path)
    superpop <- translate_population(superpopulation) 
    dataset_name <- get_dataset_name(file_path)
    file_subset <- paste(dirname(file_path),"/",gene,"_",dataset_name,"_subset.txt",sep="") 
    awk_cmd <- paste("awk -F \"",sep,"\" 'NR==1{print $0}{ if(($",colDict[chrom_col]," == ",gene_chr,")",
                     " && ($",colDict[position_col]," >= ",gene_start," && $",colDict[position_col]," <= ",gene_end,")) { print } }' ",file_path,
                     " > ",file_subset,sep="")
    system(awk_cmd)
  } 
}



## Get Flanking SNPs

# Get all genes surrounding the index SNP (default is 500kb upstream + 500kb downstream)
# 1000000 bp

column_dictionary <- function(file_path){
  # Get the index of each column name
  cNames <- colnames(fread(file_path, nrows = 2))
  colDict <- setNames(1:length(cNames), cNames  )
  return(colDict)
}



get_flanking_SNPs <- function(gene, top_SNPs, bp_distance=500000, file_path,
                              chrom_col="CHR", position_col="POS", snp_col="SNP",
                              pval_col="P", effect_col="Effect", stderr_col="StdErr", superpopulation="",
                              minPos=NULL, maxPos=NULL, file_sep="\t"){
  colDict <- column_dictionary(file_path)
  if(gene %in% top_SNPs$Gene==F){
    cat("----- Could not find gene '",gene,"' in topSNPs dataframe. Try a different gene name. ----- \n" )
  }else{
    topSNP_sub <- top_SNPs[top_SNPs$Gene==gene & !is.na(top_SNPs$Gene),] #[1,]
    if(is.null(minPos)){minPos <- as.numeric(topSNP_sub$POS) - bp_distance}
    if(is.null(maxPos)){maxPos <- as.numeric(topSNP_sub$POS) + bp_distance} 
    cat("---Min snp position:",minPos, "---\n")
    cat("---Max snp position:",maxPos, "---\n")
    # Specify subset file name 
    dataset_name <- get_dataset_name(file_path)
    file_subset <- paste(dirname(file_path),"/",gene,"_",dataset_name,"_subset.txt",sep="") 
    if(file.exists(file_subset)){
      cat("Subset file already exists. Importing",file_subset,"...\n")
    } else {
      # Extract subset with awk
      cat("Extracting relevant variants from fullSS...")
      start <- Sys.time()
      
      awk_cmd <- paste("awk -F \"",file_sep,"\" 'NR==1{print $0}{ if(($",colDict[chrom_col]," == ",topSNP_sub$CHR,")",
                       " && ($",colDict[position_col]," >= ",minPos," && $",colDict[position_col]," <= ",maxPos,")) { print } }' ",file_path,
                       " > ",file_subset,sep="")
      cat("\n",awk_cmd)
      system(awk_cmd)
      end <- Sys.time()
      cat("Extraction completed in", round(end-start, 2),"seconds \n")
    }
    query <- fread(file_subset, header=T, stringsAsFactors = F, sep = file_sep)
    if(dim(query)[1]==0){ 
      file.remove(file_subset)
      stop("\n Could not find any rows in full data that matched query :(") 
    }else{  
      if(stderr_col=="calculate"){
        cat("Calculating Standard Error...\n")
        query$StdErr <- subset(query, select=effect_col) / subset(query, select="statistic")
        stderr_col="StdErr"
      }
      geneSubset <- subset(query, select=c(chrom_col,position_col, snp_col, pval_col, effect_col, stderr_col) )
      geneSubset <- geneSubset %>% dplyr::rename(CHR=chrom_col,POS=position_col, SNP=snp_col,
                                                 P=pval_col, Effect=effect_col, StdErr=stderr_col) %>%
        mutate(Location=paste(CHR,":",POS,sep=""))
      ## Remove SNPs with NAs in stats
      geneSubset[(geneSubset$P<=0)|(geneSubset$P>1),"P"] <- 1
      # Get just one SNP per location (just pick the first one)
      geneSubset <- geneSubset %>% group_by(Location) %>% slice(1)
      # Mark lead SNP
      geneSubset$leadSNP <- ifelse(geneSubset$SNP==topSNP_sub$SNP, T, F)
      # Only convert to numeric AFTER removing NAs (otherwise as.numeric will turn them into 0s)
      geneSubset <- geneSubset  %>%
        mutate(Effect=as.numeric(Effect), StdErr=as.numeric(StdErr), P=as.numeric(P))
      return(geneSubset)
    }
   
  }
} 
# flankingSNPs <- get_flanking_SNPs(gene="PTK2B", top_SNPs,
#                                   file_path = Data_dirs$MESA_CAU$topSS,
#                                   chrom_col = "chr",  pval_col = "pvalue", snp_col = "snps",
#                                   effect_col = "beta", position_col = "pos_snps",
#                                   stderr_col = "calculate")
#


# susieR

## Gaston LD
#
# * Nalls et al. 2018: Imputation Panel Notes
# + _"One of the limitations of this study is the use of multiple imputation panels, due to logistic constraints.
# Adding datasets from non-European populations would be helpful to further improve our granularity in association
# testing and ability to fine-map loci through integration of more variable LD signatures."_
# + _"Post-Chang 23andMe samples were imputed using a combination of Finch for phasing (an in-house developed fork of Beagle)
# and miniMac2 for imputation with all-ethnicity samples from the __September 2013 release of 1000 Genomes Phase1__
# as reference haplotypes."_
# + _"The Nalls et al . 2014 and Chang et al . 2017 samples were imputed with Minimac2 using
# __1000 Genomes phase 1 haplotypes__.
# All additional sample series except for the post-Chang et al . 2017 samples from 23andMe were imputed using the
# __Haplotype Reference Consortium (HRC)__  on the University of Michigan imputation server under default settings
# with Eagle v2.3 phasing based on reference panel HRC r1.1 2016"_

download_all_vcfs <- function(vcf_folder="../1000_Genomes_VCFs"){
  # PHASE 3 DATA
  path3 <- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
  for(chrom in c(1:22)){
    cat("\nDownloading Chromosome",chrom,"\n")
    URL <- paste("ALL.chr",chrom,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",sep = "")
    system(paste("wget -P",file.path(vcf_folder,"Phase3"), file.path(path3, URL) ))
  }
  X_chrom <-"ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz"
  system(paste("wget -P",file.path(vcf_folder,"Phase3"), file.path(path3, X_chrom)))
  Y_chrom <- "ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz"
  system(paste("wget -P",file.path(vcf_folder,"Phase3"), file.path(path3, Y_chrom) ))
  
  popDat_URL = file.path(path3, "integrated_call_samples_v3.20130502.ALL.panel")
  popDat <- read.delim(popDat_URL, header = F, row.names = NULL)
  write.table(popDat,file=file.path(vcf_folder,"Phase3","integrated_call_samples_v3.20130502.ALL.panel"), row.names = F, sep="\t", quote = F, col.names = F)
  
  # PHASE 1 DATA
  path1 <- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521"
  for(chrom in c(1:22)){
    cat("\nDownloading Chromosome",chrom,"\n")
    URL <- paste("ALL.chr",chrom, ".phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz", sep="")
    system(paste( "wget -P",file.path(vcf_folder,"Phase1"), file.path(path1, URL) ))
  }
  X_chrom <- "ALL.chrX.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz"
  system( paste("wget -P",file.path(vcf_folder,"Phase1"), file.path(path1, X_chrom)) )
  
  popDat_URL = file.path(path1, "phase1_integrated_calls.20101123.ALL.panel")
  popDat <- read.delim(popDat_URL, header = F, row.names = NULL)
  write.table(popDat,file=file.path(vcf_folder,"Phase1","phase1_integrated_calls.20101123.ALL.panel"),  row.names = F, sep="\t", quote = F, col.names = F)
}



gaston_LD <- function(flankingSNPs, gene, reference="1KG_Phase1", superpopulation="EUR", vcf_folder=F){
  # Download portion of vcf from 1KG website
  region <- paste(flankingSNPs$CHR[1],":",min(flankingSNPs$POS),"-",max(flankingSNPs$POS), sep="")
  chrom <- flankingSNPs$CHR[1]
  # PHASE 3 DATA
  if(reference=="1KG_Phase3"){
    cat("LD Reference Panel = 1KG_Phase3 \n")
    if(vcf_folder==F){## With internet
      vcf_URL <- paste("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr",chrom,
                       ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",sep="")
      popDat_URL = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
    }else{## WithOUT internet
      vcf_URL <- paste(vcf_folder, "/ALL.chr",chrom,
                       ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",sep="")
      popDat_URL = file.path(vcf_folder,"integrated_call_samples_v3.20130502.ALL.panel")
    }
    
    # PHASE 1 DATA
  } else if (reference=="1KG_Phase1") {
    cat("LD Reference Panel = 1KG_Phase1 \n")
    if(vcf_folder==F){## With internet
      vcf_URL <- paste("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr",chrom,
                       ".phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz", sep="")
      popDat_URL = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/phase1_integrated_calls.20101123.ALL.panel"
    }else{## WithOUT internet
      vcf_URL <- paste(vcf_folder,"/ALL.chr",chrom,
                       ".phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz", sep="")
      popDat_URL = file.path(vcf_folder, "phase1_integrated_calls.20101123.ALL.panel")
    }
  }
  
  phase <- gsub("1KG_","",reference) 
  popDat <- read.delim(popDat_URL, header = F, row.names = NULL)
  colnames(popDat) <- c("sample","population","superpop","gender")
  
  # library(Rsamtools); #BiocManager::install("Rsamtools")
  subset_vcf <- file.path("../1000_Genomes_VCFs", phase, paste(gene,"subset.vcf",sep="_"))
  # Create directory if it doesn't exist
  if(!dir.exists(dirname(dirname(subset_vcf))) ) {
    dir.create(path = dirname(subset_vcf),recursive =  T)
  }else(cat("Creating '../1000_Genomes_VCFs' directory.\n"))
  # Download and subset vcf if the subset doesn't exist already
  if(!file.exists(subset_vcf)){
    tabix_cmd <- paste("tabix -fh",vcf_URL, region, ">", subset_vcf)
    cat(tabix_cmd)
    system(tabix_cmd)
    vcf_name <- paste(basename(vcf_URL), ".tbi",sep="")
    # file.remove(vcf_name)
  }else{cat("Identified matching VCF subset file. Importing...", subset_vcf,"\n")} 
  # Import w/ gaston and further subset
  bed.file <- read.vcf(subset_vcf)
  ## Subset rsIDs
  bed <- select.snps(bed.file, id %in% flankingSNPs$SNP)
  # Subset Individuals
  selectedInds <- subset(popDat, superpop == superpopulation)
  bed <- select.inds(bed, id %in% selectedInds$sample)
  # Cleanup extra files
  remove(bed.file)
  # file.remove("subset.vcf")
  
  # Calculate pairwise LD for all SNP combinations
  #### "Caution that the LD matrix has to be correlation matrix" -SuSiER documentation
  ### https://stephenslab.github.io/susieR/articles/finemapping_summary_statistics.html
  ld.x <- gaston::LD(bed, lim = c(1,ncol(bed)), measure ="r" )
  LD_matrix <- ld.x
  LD_matrix[!is.finite(LD_matrix)] <- 0
  # LD plot
  try({
    leadSNP = subset(flankingSNPs, leadSNP==T)$SNP
    lead_index = match(leadSNP, rownames(LD_matrix))
    start = lead_index-10
    end = lead_index+10
    LD.plot( LD_matrix[start:end, start:end], snp.positions = bed@snps$pos[start:end] )
  })
  # Double check subsetting
  # LD_matrix <- ld.x[row.names(ld.x) %in% flankingSNPs$SNP, colnames(ld.x) %in% flankingSNPs$SNP]
  return(LD_matrix)
}
# LD_matrix <- gaston_LD(flankingSNPs)


## susieR Function
#
# * Notes on L parameter
# + L is the expected number of causal variants
# + Increasing L increases computational time
# + L=1: Gives a good amount of variation in PIP.
# + L=2: Warns "IBSS algorithm did not converge in 100 iterations!", but gives good variation in PIP.
# + L=3: Warns "IBSS algorithm did not converge in 100 iterations!". All PIPs 1s and 0s.
# + These results seem to be at least partially dependent on whether the ethnic composition of the LD matrix.
# * Notes on variance:
#   + If 'estimate_residual_variance' = TRUE _without_ providing 'var_y' _and_ L>1, susieR will throw error:
#   __"Estimating residual variance failed: the estimated value is negative"__
# + Running susieR with 'var_y = var(b)' provides _exactly_ the same results.
# * Statistical Terms:
#   + posterior inclusion probability (PIP)
# + coefficient estimate (Beta)
# + Effect allele frequency (EAF)
# + The I^2 statistic describes the percentage of variation across studies that seems not to be due to chance.
susie_on_gene <- function(gene, top_SNPs,
                          bp_distance=500000, file_path, num_causal=1,
                          chrom_col="CHR", position_col="POS", snp_col="SNP",
                          pval_col="P", effect_col="Effect", stderr_col="StdErr",
                          LD_reference="1KG_Phase1", superpopulation="EUR", vcf_folder=F,
                          minPos=NULL, maxPos=NULL, file_sep="\t"){ 
    cat("\n + Extracting SNPs flanking lead SNP... \n")
    flankingSNPs <- get_flanking_SNPs(gene, top_SNPs, bp_distance=bp_distance, file_path=file_path,
                                      chrom_col=chrom_col, position_col=position_col, snp_col=snp_col,
                                      pval_col=pval_col, effect_col=effect_col, stderr_col=stderr_col, 
                                      superpopulation=superpopulation, 
                                      minPos=minPos, maxPos=maxPos, file_sep=file_sep)
  
 
  ### Get LD matrix
  cat("\n + Creating LD matrix... \n")
  LD_matrix <- gaston_LD(flankingSNPs = flankingSNPs, gene=gene, reference = LD_reference, superpopulation = superpopulation, vcf_folder = vcf_folder)
  ## Turn LD matrix into positive semi-definite matrix
  # LD_matrix2 <- ifelse(matrixcalc::is.positive.semi.definite(LD_matrix),
  #        LDmatrix,
  #        Matrix::nearPD(LD_matrix)$mat %>% as.matrix() )
  
  ## Subset summary stats to only include SNPs found in query
  geneSubset <- flankingSNPs %>% subset(SNP %in% unique(row.names(LD_matrix), colnames(LD_matrix) ) )
  geneSubset <- geneSubset[complete.cases(geneSubset),] # Remove any NAs
  LD_matrix <- LD_matrix[geneSubset$SNP,  geneSubset$SNP]
  
  b <- geneSubset$Effect
  se <- geneSubset$StdErr
  
  # Run Susie
  cat("\n + Fine mapping with SusieR... \n")
  fitted_bhat <- susie_bhat(bhat = b, shat = se,
                            R = LD_matrix,
                            n = nrow(LD_matrix),
                            
                            L = num_causal, # we assume there are at *most* 'L' causal variables
                            # scaled_prior_variance = 0.1,
                            estimate_residual_variance = T, # TRUE
                            estimate_prior_variance = T, # FALSE
                            verbose = T,
                            
                            # var_y = var(b),
                            standardize = T
  )
  # Format results
  geneSubset$Coord <- paste(geneSubset$CHR, geneSubset$POS, sep=":")
  susieDF <- data.frame(SNP=names(fitted_bhat$X_column_scale_factors), PIP=fitted_bhat$pip) %>%
    base::merge(subset(geneSubset, select=c("CHR","POS","SNP","Effect","P","Coord","leadSNP")), by="SNP") %>%
    mutate(POS=as.numeric(POS))
  # Add credible set
  cat("\n Credible Set: \n")
  try({ 
    credible_set <- geneSubset[ as.numeric(strsplit( as.character(summary(fitted_bhat)$cs$variable) ,",")[[1]]), ] 
    cat("\n ******",length(credible_set),"SNPs included in Credible Set ******\n") 
  }) 
  
  if(!exists("credible_set")){
    cat("\n ****** Could NOT identify credible set. Default to SNPs with the top 5 PIPs ******\n") 
    CS <- susieDF %>% arrange(desc(PIP))
    credible_set <- CS$SNP[1:5]
  }
  susieDF$credible_set <- ifelse(susieDF$SNP %in% credible_set$SNP, T, F)
  # cat("\n Susie Plot: Credible Set")
  # susie_plot(fitted_bhat, y="PIP", b=b, add_bar = T, add_legend = T)
  return(susieDF)
}
# susieDF <- susie_on_gene("LRRK2", top_SNPs,
#                          file_path = "Data/Parkinsons/META.PD.NALLS2014.PRS.tsv",
#                           snp_col = "MarkerName", pval_col = "P.value")
# susieDF <- susie_on_gene(gene="CLU/PTK2B", top_SNPs,
#                          file_path="Data/Alzheimers/Posthuma/phase3.beta.se.hrc.txt",
#                          effect_col = "BETA", stderr_col = "SE", position_col = "BP")

## Before-After Plots
before_after_plots <- function(gene, susieDF, before_var="P"){
  roundBreaks <- seq(plyr::round_any(min(susieDF$POS),10000), max(susieDF$POS),250000)
  
  # Label top snps
  ## BEFORE fine-mapping
  if(before_var=="Effect"){
    leadSNP_before <- subset( susieDF %>% arrange(desc(abs(Effect))), leadSNP==T)[1,]
    yLimits1 <- c(min(susieDF[before_var]), max(susieDF[before_var])*1.1 )
    y_var = "Effect"
  }else{
    # Sort by pval and then absolute Effect size
    leadSNP_before <- subset( susieDF %>% arrange(P), leadSNP==T)[1,]
    yLimits1 <- c(min(-log10(susieDF[before_var])), max(-log10(leadSNP_before[before_var]))*1.1)
    y_var = "-log10(p-value)"
  }
  leadSNP_before$type <- "before"
  leadSNP_before$color <- "red"
  # AFTER fine-mapping
  leadSNP_after = subset(susieDF, credible_set==T) %>% arrange(desc(PIP))
  leadSNP_after$type <- "after"
  leadSNP_after$color <- "darkgreen"
  leadSNP_after[leadSNP_after$PIP==max(leadSNP_after$PIP), "color"] <- "green3"
  
  labelSNPs <- rbind(leadSNP_before, leadSNP_after)
  
  # -log10(eval(parse(text=before_var)))
  before_plot <- ggplot(susieDF, aes(x=POS, y=-log10(eval(parse(text=before_var))), label=SNP, color= -log10(P) )) +
    ylim(yLimits1) +
    geom_hline(yintercept=0, alpha=.5, linetype=1, size=.5) +
    geom_point(alpha=.5) +
    geom_point(data=labelSNPs[labelSNPs$type=="before",], pch=21, fill=NA, size=4, colour=labelSNPs[labelSNPs$type=="before","color"], stroke=1) +
    geom_point(data=labelSNPs[labelSNPs$type=="after",][1,], pch=21, fill=NA, size=4, colour=labelSNPs[labelSNPs$type=="after","color"][1], stroke=1) +
    geom_segment(aes(xend=POS, yend=yLimits1[1], color= -log10(P) ), alpha=.5) +
    geom_text_repel(data=labelSNPs, aes(label=SNP), color=labelSNPs$color, segment.alpha = .5, nudge_x = .5) +
    labs(title=paste(gene," (",length(susieDF$PIP)," variants)","\nBefore Fine Mapping",sep=""),
         y=y_var, x="Position", color="-log10(p-value)") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(breaks = roundBreaks)
  
  ## After fine-mapping
  yLimits2 <- c(min(susieDF$PIP), max(susieDF$PIP)*1.1)
  
  after_plot <- ggplot(susieDF, aes(x=POS, y=PIP, label=SNP, color= -log10(P) )) +
    # ylim(yLimits) +
    geom_hline(yintercept=0,alpha=.5, linetype=1, size=.5) +
    geom_point(alpha=.5) +
    geom_point(data=labelSNPs[labelSNPs$type=="before",], pch=21, fill=NA, size=4, colour=labelSNPs[labelSNPs$type=="before","color"], stroke=1) +
    geom_point(data=labelSNPs[labelSNPs$type=="after",][1,], pch=21, fill=NA, size=4, colour=labelSNPs[labelSNPs$type=="after","color"][1], stroke=1) +
    geom_segment(aes(xend=POS, yend=yLimits2[1], color= -log10(P)), alpha=.5) +
    geom_text_repel(data=labelSNPs, aes(label=SNP), color=labelSNPs$color, segment.alpha = .5, nudge_x = .5) +
    labs(title=paste(gene," (",length(susieDF$PIP)," variants)","\nAfter Fine Mapping",sep=""), y="PIP", x="Position",
         color="-log10(p-value)") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(breaks = roundBreaks)
  
  plot_grid(before_plot, after_plot, nrow = 2) %>% print()
  # susie_plot(fitted_bhat, y="PIP", b=b, add_bar = T)
  
  createDT_html(susieDF, paste("susieR Results: ", gene), scrollY = 200)
}
# before_after_plots(gene = "LRRK2", susieDF, topVariants = 3)
# before_after_plots(gene = "CLU/PTK2B", susieDF, topVariants = 3)

## Report SNP Overlap
before_after_consensus <- function(gene, top_SNPs, susieDF, max_SNPs=10){
  # Get top SNPs from Nalls et all fine mapping
  SS_dat  <- subset(top_SNPs, Gene==gene) %>%
    arrange(P) %>% mutate(SNP=as.character(SNP))
  SumStats_SNPs <-if(max_SNPs > length(SS_dat$SNP)){ SS_dat$SNP}else{SS_dat$SNP[1:max_SNPs]}
  # Get top SNPs from susieR fine mapping
  susieR_dat <- subset(susieDF, PIP!=0) %>%
    arrange(desc(PIP)) %>%  mutate(SNP=as.character(SNP))
  susieR_SNPs <-if(max_SNPs > length(susieR_dat$SNP)){ susieR_dat$SNP}else{susieR_dat$SNP[1:max_SNPs]}
  # Calculate percent
  overlap <-  length(intersect(SumStats_SNPs, susieR_SNPs))
  percentOverlap <- overlap / length(susieR_SNPs) * 100
  cat("\n",overlap," / ",length(susieR_SNPs), " (",round(percentOverlap,2),
      "%) of SNPs of the SNPs in the summary stats were confirmed after fine-mapping.","\n",sep="")
}
# before_after_consensus(top_SNPs, susieDF, max_SNPs=10)


# Fine Map Iteratively

#Nalls_SS %>% group_by(`Nearest Gene`) %>% tally() %>% subset(n>2)
finemap_geneList <- function(top_SNPs, geneList, file_path,
                             bp_distance=500000, num_causal=1,
                             chrom_col="CHR", position_col="POS", snp_col="SNP",
                             pval_col="P", effect_col="Effect", stderr_col="StdErr",
                             LD_reference="1KG_Phase1", superpopulation="EUR",
                             topVariants=3,vcf_folder=F, 
                             minPos=NULL, maxPos=NULL, file_sep="\t", force_new_subset=T){ 
  fineMapped_topSNPs <- data.table()
  fineMapped_allResults <- data.table()
  for (gene in geneList){
    cat('\n')
    cat("###", gene, "\n")
    # Force new file to be made
    if(force_new_subset==T){
      dataset_name <- get_dataset_name(file_path)
      file_subset <- paste(dirname(file_path),"/",gene,"_",dataset_name,"_subset.txt",sep="") 
      file.remove(file_subset)
    }
    susieDF <- susie_on_gene(gene=gene, top_SNPs=top_SNPs, num_causal = 1,
                             file_path=file_path, bp_distance=bp_distance,
                             chrom_col=chrom_col, position_col=position_col, snp_col=snp_col,
                             pval_col=pval_col, effect_col=effect_col, stderr_col=stderr_col,
                             LD_reference=LD_reference, superpopulation=superpopulation,
                             minPos=minPos, maxPos=maxPos, file_sep=file_sep)
    
    before_after_plots(gene, susieDF)
    # before_after_consensus(gene, top_SNPs, susieDF, max_SNPs=10)
    # Create summary table for all genes
    newEntry <- cbind(data.table(Gene=gene), susieDF) %>% as.data.table()
    fineMapped_topSNPs <- rbind(fineMapped_topSNPs, subset(newEntry, PIP==max(PIP)) )
    fineMapped_allResults <- rbind(fineMapped_allResults, newEntry)
    cat('\n')
    createDT_html(susieDF)
  }
  createDT_html(fineMapped_topSNPs, "Potential Causal SNPs Identified by susieR", scrollY = 200)
  return(fineMapped_allResults)
}



# *********** PREPROCESSING ***********
## 23andME
add_snpIDs_to_fullSS <-function(snpInfo_path, fullSS_path){
  snp_ids <- fread(snpInfo_path,
                   select=c("all.data.id", "assay.name", "scaffold", "position"),
                   key="all.data.id")
  snp_ids$scaffold <- gsub("chr", "",snp_ids$scaffold)
  fullSS <- fread(fullSS_path,  key="all.data.id")
  dim(fullSS)
  SS_merged <- snp_ids[fullSS] %>% rename(SNP="assay.name", CHR="scaffold")
  dim(SS_merged)
  fwrite(SS_merged, "Data/Parkinsons/23andMe/PD_all_post30APRIL2015_5.2_extended.txt",
         sep = "\t", quote = F, row.names = F)
  remove(snp_ids, fullSS, SS_merged)
}
# add_snpIDs_to_fullSS(snpInfo_path = "Data/Parkinsons/23andMe/all_snp_info-5.2.txt",
#                      fullSS_path = "Data/Parkinsons/23andMe/PD_all_post30APRIL2015_5.2.dat")



##### QTL #####
# Subset eQTL data to just the parts you need
subset_eQTL_SS <- function(fullSS_path, output_path, gene, gene_col="gene_name" ){
  if(file.exists(output_path)){
    cat("Subset file already exists. Importing",output_path,"...\n")
  } else {
    colDict <- column_dictionary(fullSS_path)
    
    # Extract subset with awk
    cat("Extracting relevant variants from fullSS...\n")
    start <- Sys.time()
    awk_cmd <- paste("awk -F \"\t\" 'NR==1{print $0} $",colDict[gene_col]," == \"",gene,"\"{print $0} ' ",fullSS_path,
                     " > ",output_path, sep="")
    cat("\n",awk_cmd)
    system(awk_cmd)
    end <- Sys.time()
    cat("Extraction completed in", round(end-start, 2),"seconds \n")
  }
  # query <- fread(output_path, header=T, stringsAsFactors = F, sep = "\t")
  # # query <- read.csv.sql(file = fullSS_path, sep="\t", header = T,
  # #                     sql = paste("select * from file where ",gene_col," = '",gene,"'",sep=""))
  # cat("File subset dimensions:\n")
  # dim(query)
}

translate_population <- function(superpopulation){
  pop_dict <- list("AFA"="AFR", "CAU"="EUR", "HIS"="AMR", 
                   "AFR"="AFR","EUR"="EUR", "AMR"="AMR")
  return(pop_dict[superpopulation][[1]])
}

finemap_eQTL <- function(superpopulation, gene, fullSS_path, num_causal=1,
                         chrom_col = "chr", position_col = "pos_snps", snp_col="snps",
                         pval_col="pvalue", effect_col="beta", gene_col="gene_name", stderr_col = "calculate"){ 
  superpop <- translate_population(superpopulation)
  subset_path <- paste("Data/eQTL/MESA/",gene,"_",superpop,"_subset.txt",sep="")
  subset_eQTL_SS(fullSS_path=fullSS_path,
                 output_path=subset_path,
                 gene=gene)
  # system("wc -l Data/eQTL/MESA/CAU_eQTL_PTK2B.txt", intern = T)# Get number of rows in file
  top_SNPs <- import_topSNPs(
    file_path = subset_path,
    chrom_col = chrom_col, position_col = position_col, snp_col=snp_col,
    pval_col=pval_col, effect_col=effect_col, gene_col=gene_col,
    caption= paste(population,": eQTL Summary Stats"))
  
 
  finemapped_eQTL <- finemap_geneList(top_SNPs, geneList=gene,
                                      file_path = subset_path,
                                      chrom_col = chrom_col,  pval_col = pval_col, snp_col = snp_col,
                                      effect_col = effect_col, position_col = position_col,
                                      stderr_col = stderr_col,superpopulation = superpop,
                                      vcf_folder = vcf_folder, num_causal = num_causal)
  return(finemapped_eQTL)
}


merge_finemapping_results <- function(named_results_list, credible_sets_only=T, csv_path=F){
  final_results <- data.table()
  for(n in names(named_results_list)){
    res <-  named_results_list[n][[1]] 
    res <- cbind(Dataset=n, res)
    if(credible_sets_only==T){res <- res %>% dplyr::filter(credible_set==T)}
    final_results <- rbind(final_results, res)
  }
  if(csv_path!=F){
    fwrite(final_results, file = csv_path, quote = F, sep = ",")
  }
  return(final_results)
}



# Preprocessing
# add_snp_info <-function(snpInfo_path, fullSS_path, newSS_path){
#   # sqldf::read.csv.sql(file=snpInfo_path, header = T, sep = "\t",
#   #                     sql = paste("select * from file where location IN (", paste("'",snp_locs,"'", collapse=",",sep=""),")",sep="")
#   #                     )
#   cat("\nLoading SNP Info file...\n")
#   snp_ids <- fread(snpInfo_path,header = T, stringsAsFactors = F, key="location", colClasses = rep("character",2))
#
#   cat("\nLoading full SS file...\n")
#   fullSS <- fread(fullSS_path, stringsAsFactors = F,
#                   colClasses = c(rep("character",3), rep("numeric",7), "character", rep("numeric",4)))
#   fullSS$MarkerName <- gsub("chr","",fullSS$MarkerName)
#   fullSS <- fullSS %>% rename(location="MarkerName")
#   fullSS <- data.table(fullSS, key="location")
#   dim(fullSS)
#   cat("\nMerging files...\n")
#   SS_merged <- snp_ids[fullSS]
#   dim(SS_merged)
#   cat("\nWriting new file...\n")
#   fwrite(SS_merged, newSS_path, sep = "\t", quote = F, row.names = F)
#   remove(snp_ids, fullSS, SS_merged)
# }
# add_snp_info(snpInfo_path = "Data/HRC.RSID.txt",
#              fullSS_path = "Data/Parkinsons/Nalls_2019/nallsEtal2019_no23andMe.tab.txt",
#              newSS_path="Data/Parkinsons/Nalls_2019/nallsEtal2019_no23andMe_extended.txt")

#
# split_location_col <- function(input_path, output_path){
#   awk_cmd <- paste("awk 'BEGIN {print \"CHR\tPOS\tRSID\tA1\tA2\tfreq\tbeta\tse\tp\tN_cases\tN_controls\"}; FNR>1 {split($1, c, \":\"); print c[1] \"\t\" c[2] \"\t\" $2 \"\t\" $4 \"\t\" $5 \"\t\" $6 \"\t\" $7 \"\t\" $8 \"\t\" $9 \"\t\" $10 \"\t\" $11}' "
#          input_path," > ", output_path, sep="")
#   cat(awk_cmd)
#   system(awk_cmd)
# }
