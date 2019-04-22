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
library(gaston)
library(tidyr)
library(biomaRt)
# library(bigsnpr) # BiocManager::install("bigsnpr")

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

source("./echolocatoR/R/directory_functions.R")
source("./echolocatoR/R/query_functions.R")
source("./echolocatoR/R/LD_functions.R")
source("./echolocatoR/R/plot_functions.R")

reload <- function(){
  source("echolocatoR/R/functions.R") 
}

get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}
# get_os() 

# quick_start
quickstart <- function(){ 
  # DO NOT include "./" in front of file paths. Confuses command line functions.
  # Assign global variables to test functions
  gene <<- "LRRK2"
  leadSNP <<- "rs76904798" 
  chrom_col <<- "CHR"
  position_col <<- "POS"
  snp_col <<- "RSID"
  pval_col <<- "p"
  effect_col <<- "beta"
  stderr_col <<- "se"
  vcf_folder <<- F#"../1000_Genomes_VCFs/Phase1"
  superpopulation <<- "EUR"
  force_new_subset <<- T 
  minPos <<- NULL
  maxPos <<- NULL
  file_sep <<- "\t"
  min_r2 <<- .2
  LD_block <<- F
  block_size <<- .7
  min_Dprime <<- 1
  plink_folder <<- "./plink_tmp" 
  reference <<- "1KG_Phase1"
  LD_reference <<- "1KG_Phase1"
  bp_distance <<- 500000
  num_causal <<- 1
  vcf_URL <<- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
  popDat_URL <<- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
  chr <<- 8
  vcf_name <<- "ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi"
  query_by <<- "coordinates"
  topSNP_sub <<- top_SNPs[top_SNPs$Gene==gene & !is.na(top_SNPs$Gene),]
  location_sep <<- ":"
  query_by <- ""
  
  fullSS_path <<- "Data/Parkinsons/Nalls_23andMe/nallsEtAl2019_allSamples_allVariants.mod.txt" 
  subset_path <<- get_subset_path(fullSS_path, gene, superpopulation)
  colDict <<- column_dictionary(fullSS_path)
  if(is.null(minPos)){minPos <<- topSNP_sub$POS - bp_distance}
  if(is.null(maxPos)){maxPos <<- topSNP_sub$POS + bp_distance} 
}
# quickstart()


quickstart_eQTL <- function(){ 
  # Assign global variables to test functions
  gene <<- "LRRK2"
  leadSNP <<- "rs76904798" 
  geneList <<- c("LRRK2")
  chrom_col <<- "chr"
  position_col <<- "pos_snps"
  snp_col <<- "snps"
  pval_col <<- "pvalue"
  effect_col <<- "beta"
  gene_col <<- "gene_name"
  stderr_col <<- "calculate"
  fullSS_path <<- "Data/eQTL/MESA/CAU_cis_eqtl_summary_statistics.txt"
  subset_path <<- "Data/eQTL/MESA/LRRK2_MESA_CAU_subset.txt"
  vcf_folder <<- F#"../1000_Genomes_VCFs/Phase1"
  superpopulation <<- "CAU"
  force_new_subset <<- F 
  minPos <<- NULL
  maxPos <<- NULL
  file_sep <<- "\t"
  min_r2 <<- .2
  LD_block <<- F
  block_size <<- .7
  min_Dprime <<- 1
  plink_folder <<- "./plink_tmp" 
  reference <<- "1KG_Phase1"
  LD_reference <<- "1KG_Phase1"
  bp_distance <<- 500000
  num_causal <<- 1
  vcf_URL <<- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
  popDat_URL <<- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
  chr <<- 8
  vcf_name <<- "ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi"
}
# quickstart_eQTL()

## ---------------- General Functions ----------------  ##

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
  dataset_name <- tail(strsplit(dirname(file_path), "/")[[1]],n = 1)
  return(dataset_name)
}

delete_subset <- function (force_new_subset, subset_path){
  # Force new file to be made
  if(force_new_subset==T){
    # dataset_name <- get_dataset_name(subset_path)
    # subset_path <- paste(dirname(subset_path),"/",gene,"_",superpopulation,"_",dataset_name,"_subset.txt",sep="") 
    file.remove(subset_path)
  }
}




## ---------------- Fine-mapping Functions ----------------  ##


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
susie_on_gene <- function(gene, fullSS_path, subset_path="auto", top_SNPs="auto", 
                          bp_distance=500000, num_causal=1,
                          chrom_col="CHR", position_col="POS", snp_col="SNP",
                          pval_col="P", effect_col="Effect", stderr_col="StdErr", 
                          tstat_col="t-stat", gene_col="Gene",
                          LD_reference="1KG_Phase1", superpopulation="EUR", vcf_folder=F,
                          minPos=NULL, maxPos=NULL, file_sep="\t", 
                          min_r2=0, LD_block=F, block_size=.7, min_Dprime=F, query_by="coordinates"){  
   cat("\n + Extracting SNPs flanking lead SNP... \n")
   flankingSNPs <- extract_SNP_subset(gene=gene, top_SNPs=top_SNPs, subset_path=subset_path, fullSS_path=fullSS_path,
                                      chrom_col=chrom_col, position_col=position_col, snp_col=snp_col,
                                      pval_col=pval_col, effect_col=effect_col, stderr_col=stderr_col,
                                      gene_col=gene_col, tstat_col=tstat_col, bp_distance=bp_distance,
                                      superpopulation=superpopulation,  
                                      minPos=minPos, maxPos=maxPos, file_sep=file_sep, 
                                      query_by=query_by)
  ### Get LD matrix
  cat("\n + Creating LD matrix... \n")
  LD_matrix <- gaston_LD(flankingSNPs = flankingSNPs, gene=gene, reference = LD_reference, 
                         superpopulation = superpopulation, vcf_folder = vcf_folder,
                         min_r2 = min_r2, LD_block=LD_block, block_size=block_size, min_Dprime=min_Dprime)
  ## Subset summary stats to only include SNPs found in query
  geneSubset <- flankingSNPs %>% subset(SNP %in% unique(row.names(LD_matrix), colnames(LD_matrix) ) )
  geneSubset <- geneSubset[complete.cases(geneSubset),] # Remove any NAs
  LD_matrix <- LD_matrix[row.names(LD_matrix) %in% geneSubset$SNP,  colnames(LD_matrix) %in% geneSubset$SNP]  
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
  ## *** IMPORTANT! ***: In case susieR cannot identify any credible set, 
  # take the snps with the top 5 PIPs and provide a warning message. Interpret these snps with caution.
  check_credible <- function(geneSubset, fitted_bhat){
    credible_set <- geneSubset[ as.numeric(strsplit( as.character(summary(fitted_bhat)$cs$variable) ,",")[[1]]), ]$SNP 
    cat("\n ******",length(credible_set),"SNPs included in Credible Set ******\n") 
    return(credible_set)
  }  
  error_handling <- function(code) {
    tryCatch(code,
             error = function(c) {
               cat("\n--- ERROR ---")
               cat("\n ****** Could NOT identify credible set. Default to SNPs with the top 5 PIPs ******\n") 
               CS <- susieDF %>% arrange(desc(PIP))
               credible_set <- as.character(CS[1:5,]$SNP)
               return(credible_set)
               },
             warning = function(c) "Warning",
             message = function(c) "Message"
    )
  }
  
  credible_set <- error_handling(check_credible(geneSubset, fitted_bhat))
  susieDF$credible_set <- ifelse(susieDF$SNP %in% credible_set, T, F)
  # cat("\n Susie Plot: Credible Set")
  # susie_plot(fitted_bhat, y="PIP", b=b, add_bar = T, add_legend = T)
  return(susieDF)
}
# susieDF <- susie_on_gene("LRRK2", top_SNPs,
#                          file_path = "Data/Parkinsons/META.PD.NALLS2014.PRS.tsv",
#                           snp_col = "MarkerName", pval_col = "P.value") 

# Fine Map Iteratively 
finemap_geneList <- function(geneList, fullSS_path, subset_path="auto", top_SNPs="auto",
                             bp_distance=500000, num_causal=1,
                             chrom_col="CHR", position_col="POS", snp_col="SNP",
                             pval_col="P", effect_col="Effect", stderr_col="StdErr",
                             tstat_col = "t-stat", gene_col="Gene",
                             LD_reference="1KG_Phase1", superpopulation="EUR",
                             topVariants=3, vcf_folder=F, 
                             minPos=NULL, maxPos=NULL, file_sep="\t", force_new_subset=T,
                             min_r2=0, LD_block=F, block_size=.7, min_Dprime=F, 
                             query_by="coordinates"){ 
  fineMapped_topSNPs <- data.table()
  fineMapped_allResults <- data.table()
  for (gene in geneList){
    start_gene <- Sys.time() 
    cat('\n')
    cat("###", gene, "\n")
    subset_path <- get_subset_path(fullSS_path=fullSS_path, gene=gene, subset_path = subset_path )
    delete_subset(force_new_subset, subset_path)
    susieDF <- susie_on_gene(gene=gene, top_SNPs=top_SNPs, subset_path=subset_path, fullSS_path=fullSS_path, 
                             num_causal=num_causal, bp_distance=bp_distance,
                             chrom_col=chrom_col, position_col=position_col, snp_col=snp_col,
                             pval_col=pval_col, effect_col=effect_col, stderr_col=stderr_col,
                             tstat_col=tstat_col, gene_col=gene_col,
                             LD_reference=LD_reference, superpopulation=superpopulation,
                             minPos=minPos, maxPos=maxPos, file_sep=file_sep, 
                             min_r2=min_r2,LD_block=LD_block, block_size=block_size, min_Dprime=min_Dprime, 
                             query_by=query_by)
    
    before_after_plots(gene, susieDF)
    # before_after_consensus(gene, top_SNPs, susieDF, max_SNPs=10)
    # Create summary table for all genes
    newEntry <- cbind(data.table(Gene=gene), susieDF) %>% as.data.table()
    fineMapped_topSNPs <- rbind(fineMapped_topSNPs, subset(newEntry, PIP==max(PIP)) )
    fineMapped_allResults <- rbind(fineMapped_allResults, newEntry)
    cat('\n')
    createDT_html(susieDF)
    end_gene <- Sys.time()
    cat("\n",gene,"fine-mapped in", round(end_gene-start_gene, 2),"seconds \n")
  }
  createDT_html(fineMapped_topSNPs, "Potential Causal SNPs Identified by susieR", scrollY = 200)
  return(fineMapped_allResults)
}



# 
# 
# ##### QTL #####
# finemap_eQTL <- function(superpopulation, geneList, fullSS_path, num_causal=1,
#                          chrom_col="CHR", position_col="POS", snp_col="SNP",
#                          pval_col="P", effect_col="Effect", stderr_col="StdErr",
#                          gene_col="gene_name", LD_reference="1KG_Phase1",
#                          force_new_subset=T, minPos=NULL, maxPos=NULL,
#                          min_r2=0, LD_block=F, block_size=.7, min_Dprime=F,
#                          bp_distance=500000){
#   # Setup variables
#   superpop <- translate_population(superpopulation)
#   dataset <- get_dataset_name(fullSS_path)
#   fineMapped_topSNPs <- data.table()
#   fineMapped_allResults <- data.table()
#   # Iterate over genes
#   for (gene in geneList){
#     start_gene <- Sys.time()
#     cat('\n')
#     cat("###", gene, "\n")
#     delete_subset(force_new_subset=force_new_subset, subset_path=subset_path)
#     # Subset data
#     subset_path <- paste("Data/eQTL/MESA/",gene,"_",dataset,"_",superpopulation,"_subset.txt",sep="")
#     top_SNPs <- subset_eQTL_SS(fullSS_path=fullSS_path,
#                    subset_path=subset_path,
#                    gene=gene, force_new_subset=force_new_subset,
#                    gene_col = gene_col)
# 
#     # Fine-map
#     susieDF <- susie_on_gene(gene=gene, subset_path=subset_path,
#                              top_SNPs=top_SNPs, num_causal = num_causal,
#                              subset_path=subset_path, bp_distance=bp_distance,
#                              chrom_col=chrom_col, position_col=position_col, snp_col=snp_col,
#                              pval_col=pval_col, effect_col=effect_col, stderr_col=stderr_col,
#                              LD_reference=LD_reference, superpopulation=superpopulation,
#                              minPos=minPos, maxPos=maxPos, file_sep=file_sep, min_r2=min_r2,
#                              LD_block=LD_block, block_size=block_size, min_Dprime=min_Dprime)
#     # Plot
#     before_after_plots(gene, susieDF)
#     # Create summary table for all genes
#     newEntry <- cbind(data.table(Gene=gene), susieDF) %>% as.data.table()
#     fineMapped_topSNPs <- rbind(fineMapped_topSNPs, subset(newEntry, PIP==max(PIP)) )
#     fineMapped_allResults <- rbind(fineMapped_allResults, newEntry)
#     cat('\n')
#     createDT_html(susieDF)
#     end_gene <- Sys.time()
#     cat("\n",gene,"fine-mapped in", round(end_gene-start_gene, 2),"seconds \n")
#   }
#   createDT_html(fineMapped_topSNPs, "Potential Causal SNPs Identified by susieR", scrollY = 200)
#   return(fineMapped_allResults)
# }
# 



merge_finemapping_results <- function(allResults, credible_sets_only=T, csv_path=F){
  final_results <- data.table()
  for(n in names(allResults)){ 
    try({ 
      res <- cbind(Dataset=n, allResults[[n]])
      if(credible_sets_only==T){res <- res %>% dplyr::filter(credible_set==T)}
      final_results <- rbind(final_results, res)
    })
  }
  if(csv_path!=F){
    fwrite(final_results, file = csv_path, quote = F, sep = ",")
  }
  return(final_results)
}


get_subset_path <- function(fullSS_path, gene, superpopulation="EUR", subset_path="auto"){
  # Specify subset file name 
  if(subset_path=="auto"){ 
    dataset_name <- get_dataset_name(fullSS_path)
    created_sub_path <- paste(dirname(fullSS_path),"/",gene,"_",superpopulation,"_",dataset_name,"_subset.txt",sep="") 
    return(created_sub_path)
  } else{return(subset_path)}
}
 
# run_finemapping_tools <- function(X, Y, Z=NULL){
#   library(abind)
#   mm_regression = function(X, Y, Z=NULL) {
#     if (!is.null(Z)) {
#       Z = as.matrix(Z)
#     }
#     reg = lapply(seq_len(ncol(Y)), function (i) simplify2array(susieR:::univariate_regression(X, Y[,i], Z)))
#     reg = do.call(abind, c(reg, list(along=0)))
#     # return array: out[1,,] is betahat, out[2,,] is shat
#     return(aperm(reg, c(3,2,1)))
#   }
#   sumstats = mm_regression(as.matrix(X), as.matrix(Y))
#   saveRDS(list(data=dat, sumstats=sumstats), 'N2.with_sumstats.rds')
#   
#   system(cat("export PATH=~/GIT/github/mvarbvs/dsc/modules/linux:$PATH
#   Rscript ~/GIT/susieR/inst/code/finemap.R input=\"N2.with_sumstats.rds\" output=\"N2finemapping.FINEMAP\" args=\"--n-causal-max\ 2\" &> /dev/null"))
#   
#   system("export PATH=~/GIT/github/mvarbvs/dsc/modules/linux:$PATH
#   Rscript ~/GIT/susieR/inst/code/caviar.R input=\"N2.with_sumstats.rds\" output=\"N2finemapping.CAVIAR\" args=\"-c\ 2\ -g\ 0.01\" &> /dev/null")
#   
#   system("export PATH=~/GIT/github/mvarbvs/dsc/modules/linux:$PATH
#   python ~/GIT/susieR/inst/code/dap-g.py N2.with_sumstats.rds N2finemapping.DAP -ld_control 0.20 --all &> /dev/null")
# 
# }


# Turn LD matrix into positive semi-definite matrix
# LD_matrix_check <- function(LD_matrix){
#   mat <- LD_matrix %>% data.table() %>% as.matrix()
#   mat <- Matrix::forceSymmetric(mat)
#   if(matrixcalc::is.positive.semi.definite(mat)){
#     return(LD_matrix)
#   } else{ 
#     # Makes matrix symmetric and positive definite 
#     ## (but super slow...)
#     mat <- Matrix::nearPD(mat)$mat %>% as.matrix() 
#   } 
# }  

