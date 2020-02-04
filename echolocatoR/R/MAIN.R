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
# `              v  v           `
#
#   message("      =/\                  /\=    ")
#   message("      / \'._   (\_/)   _.'/ \     ")
#   message("     / .''._'--(o.o)--'_.''. \    ")
#   message("    /.' _/ |`'=/ { \='`| \_ `.\   ")
#   message("   /` .' `\;-,'\___/',-;/` '. '\  ")
#   message("  /.-'       `\(-V-)/`       `-.\ ")
#   message(" `              v  v           `  ")
#


# You can learn more about package authoring with RStudio at:
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

# regex commands: http://www.endmemo.com/program/R/gsub.php

# Load libraries
# .libPaths()

# library(R.utils)
# library(devtools)
# library(readxl)
# library(DT)
library(data.table)
library(dplyr)
library(ggplot2)
library(plotly)
library(patchwork) #devtools::install_github("thomasp85/patchwork")
# library(cowplot)
# library(ggrepel)
# library(curl)
# library(gaston)
# library(tidyr)
# library(BiocManager)
# library(biomaRt) # BiocManager::install("biomaRt")
# library(refGenome)
# library(crayon)

# library(coloc)
# install.packages("haploR", dependencies = TRUE)
# library(haploR)
# library(GeneOverlap) #BiocManager::install("GeneOverlap")
# library(rtracklayer) #BiocManager::install("rtracklayer")
# BiocManager::install(c("supraHex","graph","Rgraphviz","dnet"))
# library(XGR)# install.packages("XGR")



# library(Rsamtools) # BiocManager::install("Rsamtools")

# install_fGWAS <- function(){
#   devtools::install_github("wzhy2000/fGWAS/pkg")
#   system("git clone https://github.com/wzhy2000/fGWAS.git & cd fGWAS & R CMD INSTALL pkg")
# }
# library(fGWAS)
# library(snpStats) #BiocManager::install("snpStats")
# library(bigsnpr) # BiocManager::install("bigsnpr")

# *** SUSIE ****
# library(knitrBootstrap) #install_github('jimhester/knitrBootstrap')
# library(susieR) # devtools::install_github("stephenslab/susieR")

# *** finemapr ****
## finemapr contains: finemap, CAVIAR, and PAINTOR
# library(finemapr) # devtools::install_github("variani/finemapr")
# library(roxygen2) #roxygenize()

# *** locuscomparer ****
# https://github.com/boxiangliu/locuscomparer
# library(locuscomparer); #devtools::install_github("boxiangliu/locuscomparer")

# thm <- knitr::knit_theme$get("bipolar")
# knitr::knit_theme$set(thm)

# General
source("./echolocatoR/R/directory.R")
# Download
source("./echolocatoR/R/Download/query.R")
source("./echolocatoR/R/Download/standardize.R")
source("./echolocatoR/R/Download/tabix.R")
source("./echolocatoR/R/Download/downloaders.R")
# LD
source("./echolocatoR/R/LD/LD.R")
source("./echolocatoR/R/LD/UKBiobank_LD.R")
# Fine-map
source("./echolocatoR/R/Finemap/multi_finemap.R")
source("./echolocatoR/R/Finemap/SUSIE.R")
source("./echolocatoR/R/Finemap/ABF.R")
source("./echolocatoR/R/Finemap/FINEMAP.R")
source("./echolocatoR/R/Finemap/PAINTOR.R")
source("./echolocatoR/R/Finemap/COJO.R")
source("./echolocatoR/R/Finemap/COLOC.R")
source("./echolocatoR/R/Finemap/fGWAS.R")
source("./echolocatoR/R/Finemap/POLYFUN.R")
# Plot
source("./echolocatoR/R/Plot/plot.R")
source("./echolocatoR/R/Plot/ggbio.R")
source("./echolocatoR/R/Plot/QTL_boxplots.R")
# Annotate
source("./echolocatoR/R/Annotate/annotate.R")
source("./echolocatoR/R/Annotate/GoShifter.R") # ***
source("./echolocatoR/R/Annotate/XGR.R")
source("./echolocatoR/R/Annotate/psychENCODE.R")
source("./echolocatoR/R/Annotate/mergeQTL.R")
source("./echolocatoR/R/Annotate/Nott_2019.R")

# When there's a ton of files, turn off indexing to speed up Rstudio:
# https://stackoverflow.com/questions/14599359/rsession-cpu-usage-when-idle


startup_image <- function(){
  try({
    col.text <- function(txt){
      library(dplyr)
      c(txt,"\n") %>%
        crayon::blurred() %>%
        crayon::bgBlack() %>%
        # crayon::col_align(align = "left") %>%
        crayon::cyan() %>%
        cat()
    }
    col.text("))))))))))>>))))))))))>  E c h o l o c a t o R  <((((((((((<<((((((((((")
    col.text("")
    col.text("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ V1.0 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    col.text("~~~~~~~~~~~~~~~~~~~~~~ © 2019 - Brian M. Schilder ~~~~~~~~~~~~~~~~~~~~~")
    col.text("~Department of Neuroscience, Department of Genetics & Genomic Sciences~")
    col.text("~~~~~~~~~~Icahn School of Medicine at Mount Sinai, NY, NYC, USA~~~~~~~~")
    col.text("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    grid::grid.newpage()
    img <- png::readPNG("./echolocatoR/images/echo_logo.png")
    grid::grid.raster(img)
  })
}
# startup_image()


reload <- function(){
  source("echolocatoR/R/MAIN.R")
}


rbind.file.list <- function(file.list, verbose=T){
  merged.dat <- lapply(file.list, function(x){
    printer(x, v = verbose)
    dat <- data.table::fread(x, nThread = 4)
    return(dat)
  }) %>% data.table::rbindlist(fill=T)
  return(merged.dat)
}

# data.table::fwrite( subset(subset_DT, SNP != "rs34637584"), "Data/GWAS/Nalls23andMe_2019/LRRK2/LocusZoomData_-rs34637584.txt", sep="\t")

quick_finemap <- function(locus="LRRK2", consensus_thresh = 3){
  gene <<- locus
  # locus <<- locus
  results_path <<- file.path("./Data/GWAS/Nalls23andMe_2019",gene)
  finemap_DT <<- data.table::fread(file.path(results_path, "Multi-finemap/Multi-finemap_results.txt"))
  finemap_DT <<- find_consensus_SNPs(finemap_DT, consensus_thresh = consensus_thresh)

  if(file.exists(file.path(results_path,"plink","UKB_LD.RDS"))){
    LD_matrix <<- readRDS(file.path(results_path,"plink","UKB_LD.RDS"))
  }
  if(file.exists(file.path(results_path,"plink","LD_matrix.RData")) ){
    file.path(results_path,"plink","LD_matrix.RData")
  }
  sub.out <- subset_common_snps(LD_matrix, finemap_DT)
  LD_matrix <<- sub.out$LD
  finemap_DT <<- sub.out$DT
  subset_DT <<- finemap_DT
}

effective_sample_size <- function(finemap_DT){
  finemap_DT$N <- (4.0 / (1.0/finemap_DT$N_cases + 1.0/finemap_DT$N_controls) )
  sample_size <- as.integer(median(finemap_DT$N))
  return(sample_size)
}

quick_finemap_soft <- function(locus="LRRK2"){
  gene <<- locus
  # locus <<- locus
  results_path <<- file.path("./Data/GWAS/Nalls23andMe_2019",gene)
  finemap_DT <- data.table::fread(file.path(results_path, "Multi-finemap/Multi-finemap_results.txt"))
  finemap_DT <- find_consensus_SNPs(finemap_DT)
  return(finemap_DT)
}


tryFunc <- function(input, func) {
  out <- tryCatch(
    {
      # Just to highlight: if you want to use more than one
      # R expression in the "try" part then you'll have to
      # use curly brackets.
      # 'tryCatch()' will return the last evaluated expression
      # in case the "try" part was completed successfully

      func(input)
      # The return value of `readLines()` is the actual value
      # that will be returned in case there is no condition
      # (e.g. warning or error).
      # You don't need to state the return value via `return()` as code
      # in the "try" part is not wrapped insided a function (unlike that
      # for the condition handlers for warnings and error below)
    },
    error=function(cond) {
      message(paste("URL does not seem to exist:", input))
      message("Here's the original error message:")
      message(cond)
      # Choose a return value in case of error
      return(NA)
    },
    warning=function(cond) {
      message(paste("URL caused a warning:", input))
      message("Here's the original warning message:")
      message(cond)
      # Choose a return value in case of warning
      return(NULL)
    },
    finally={
      # NOTE:
      # Here goes everything that should be executed at the end,
      # regardless of success or error.
      # If you want more than one expression to be executed, then you
      # need to wrap them in curly brackets ({...}); otherwise you could
      # just have written 'finally=<expression>'
      message(paste("Processed URL:", input))
      message("Some other message at the end")
    }
  )
  return(out)
}


# quick_start
quickstart <- function(){
  # reload()
  allResults <<- list()
  gene <<- "LRRK2"
  leadSNP <<- "rs76904798"
  chrom_col <<- "CHR"
  position_col <<- "POS"
  snp_col <<- "RSID"
  pval_col <<- "p"
  effect_col <<- "beta"
  stderr_col <<- "se"
  freq_col <<- "freq"
  MAF_col<<-"calculate"
  A1_col <<- "A1"
  A2_col <<- "A2"
  finemap_method_list <<- c("SUSIE","ABF","FINEMAP","COJO")
  finemap_methods <<- c("SUSIE","FINEMAP")
  method <<- finemap_methods[1]
  method_list <<- c("SUSIE","FINEMAP")
  force_new_LD <<- F
  before_var <<- "P"

  download_reference <<- T#"../1000_Genomes_VCFs/Phase1"
  superpopulation <<- "EUR"
  force_new_subset <<- F
  min_POS <<- NA
  max_POS <<- NA
  min_MAF <<- NA
  file_sep <<- "\t"
  min_r2 <<- F
  LD_block <<- F
  block_size <<- .7
  min_Dprime <<- F
  plink_folder <<- "./Data/GWAS/Nalls23andMe_2019/LRRK2/plink"
  reference <<- "1KG_Phase1"
  LD_reference <<- "1KG_Phase1"
  bp_distance <<- 100000
  n_causal <<- 10
  vcf_URL <<- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
  popDat_URL <<- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
  # chr <<- 12
  vcf_name <<- "ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi"
  vcf_folder <<- "./Data/Reference/1000_Genomes"
  query_by <<- "coordinates"
  remove_variants <<- "rs34637584"
  remove_correlates <<- "rs34637584"

  location_sep <<- ":"
  query_by <<- ""
  dataset_name <<- "Nalls23andMe_2019"
  dataset_type <<- "GWAS"
  N_cases_col <<- "N_cases"
  N_controls_col <<- "N_controls"
  proportion_cases <<- "calculate"


  top_SNPs <- Nalls_top_SNPs <- import_topSNPs(
    topSS_path = Directory_info(dataset_name, "topSS"),
    chrom_col = "CHR", position_col = "BP", snp_col="SNP",
    pval_col="P, all studies", effect_col="Beta, all studies", gene_col="Nearest Gene",
    caption= "Nalls et al. (2018) w/ 23andMe PD GWAS Summary Stats",
    group_by_locus = T,
    locus_col = "Locus Number",
    remove_variants = "rs34637584")
  topSNP_sub <<- top_SNPs[top_SNPs$Gene==gene & !is.na(top_SNPs$Gene),]


  results_path <<- make_results_path(dataset_name, dataset_type, gene)
  subset_path <<- get_subset_path(results_path, gene)
  subset_DT <<- finemap_DT <<- data.table::fread("Data/GWAS/Nalls23andMe_2019/LRRK2/Multi-finemap/Multi-finemap_results.txt", sep="\t")
  multi <<- T
  subtitle=""
  # subset_path <<- 'Data/GWAS/Nalls23andMe_2019/LRRK2/LRRK2_combined_meta_subset.txt'
  # subset_DT <<- data.table::fread(subset_path, sep="\t")
  LD_path <<- file.path(results_path, "plink/LD_matrix.RData")

  fullSS_path <<- Directory_info(dataset_name, "fullSS.local")
  show_plot<<-T
  subtitle<<-NA
  multi<<-T
  LD_SNP<<- NA


  colDict <<- column_dictionary(fullSS_path)
  if(is.na(min_POS)){min_POS <<- topSNP_sub$POS - bp_distance}
  if(is.na(max_POS)){max_POS <<- topSNP_sub$POS + bp_distance}

  # COJO
  conditioned_snps <<- "rs76904798"#"rs34637584"
  excluded_snps <<- "rs34637584"
  min_MAF <<- 0
  GCTA_path <<- "echolocatoR/tools/gcta_1.92.1beta5_mac/bin/gcta64"
  bfiles <<- "plink_tmp/plink"
  # load("Data/GWAS/Nalls23andMe_2019/LRRK2/plink/LD_matrix.RData")
  LD_matrix <<- readRDS(file.path(results_path,"plink/UKB_LD.RDS"))
  # subset_DT <- subset(subset_DT, SNP %in% unique(row.names(LD_matrix), colnames(LD_matrix) ) )
  # subset_DT <- subset_DT[complete.cases(subset_DT),] # Remove any NAs
  # LD_matrix <- LD_matrix[row.names(LD_matrix) %in% subset_DT$SNP,  colnames(LD_matrix) %in% subset_DT$SNP]

  scaled_prior_variance <<- 0.1
  sample_size <<- NA
  freq_cutoff <<- 0.1

  finemap_method_list <<- c("SUSIE","FINEMAP")#  "ABF", "COJO"

  consensus <<- T

  dataset1_path <<- "./Data/GWAS/Nalls23andMe_2019/LRRK2/LRRK2_Nalls23andMe_2019_subset_500kb.txt"
  dataset2_path <<- "./Data/eQTL/MESA_CAU/LRRK2/LRRK2_MESA_CAU_subset.txt"
  # shared_MAF <<- data.table::fread("Data/GWAS/Nalls23andMe_2019/LRRK2/LRRK2_Nalls23andMe_2019_subset_500kb.txt.gz", sep="\t")$MAF
  plot_subtitle <<- "Fairfax (2014) + CD14 eQTL"
  dataset2_proportion_cases <<- 5e-324
  PP_threshold <<- 0.8
  force_new_finemap <<- F
  dataset1_type <<- "cc"
  dataset2_type <<- "quant"
  shared_MAF <<- 1
  which_merge <<-1
  show_plot <<- T

  goshifter_path <<- "./echolocatoR/tools/goshifter"
  permutations <<- 1000
  remove_tmps <<- T
  chromatin_states <<- c("TssA","EnhA1","EnhA2")
  diff_freq <<- 0.1

  paintor_path <<- "./echolocatoR/tools/PAINTOR_V3.0"
  locus_name<<- NULL
  GWAS_datasets <<- dataset_name
  QTL_datasets <<- NULL
  populations <<- "EUR"
  use_annotations <<- F
}


quickstart_AD <- function(locus="PTK2B"){
  loci <<- "PTK2B"
  trim_gene_limits <<- F
  dataset_name <<- dataset_name
  dataset_type <<- "GWAS"
  query_by <<-"tabix"
  finemap_method <<- c("ABF","SUSIE","POLYFUN_SUSIE","FINEMAP")
  force_new_subset <<- T
  force_new_LD <<- F
  force_new_finemap <<- T
  # file_sep <<- " "
  fullSS_path <<- Directory_info(dataset_name, "fullSS.local")
  chrom_col <<- "Chromosome"
  position_col <<- "Position"
  snp_col <<- "MarkerName"
  pval_col <<- "Pvalue"
  effect_col <<- "Beta"
  stderr_col <<- "SE"
  A1_col <<- "Effect_allele"
  A2_col <<- "Non_Effect_allele"
  N_cases <<- 21982
  N_controls <<- 41944
  proportion_cases <<- "calculate"
  MAF_col <<- "MAF"
  freq_col <<- "freq"
  N_cases_col<<-"N_cases"
  N_controls_col<<-"N_controls"
  tstat_col <<-"calculate"


  bp_distance <<- 500000
  download_reference <<- T
  LD_reference <<- "UKB"
  superpopulation <<- "EUR"
  LD_block <<- F
  min_MAF <<- 0.001
  PP_threshold <<- .95
  n_causal <<- 5
  remove_tmps <<- F
  # server <<- F
  gene <<- locus
  results_path <<- file.path("Data/GWAS/Kunkle_2019",gene)
  subset_path <<- file.path(results_path, paste0(gene,"_Kunkle_2019_subset.tsv.gz"))
  # finemap_DT <<- data.table::fread(subset_path)
  finemap_DT <<- data.table::data.table(standardize_subset(subset_path = subset_path, gene=gene))
  subset_DT <<- finemap_DT
  LD_matrix <<- readRDS(file.path(results_path,"plink/UKB_LD.RDS"))
  sub.out <- subset_common_snps(LD_matrix=LD_matrix,
                                finemap_DT=subset_DT)
  LD_matrix <- sub.out$LD
  subset_DT <- sub.out$DT

  polyfun<<-"./echolocatoR/tools/polyfun"
  force_new_priors <<- F
}





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


dt.replace <- function(DT, target, replacement){
  # Annoyingly, there is no native function for this in data.table
  for(col in names(DT)) set(DT, i=which(DT[[col]]==target), j=col, value=replacement)
  return(DT)
}


subset_common_snps <- function(LD_matrix, finemap_DT, verbose=F){
  printer("+ Subsetting LD matrix and finemap_DT to common SNPs...", v=verbose)
  # Remove duplicate SNPs
  dups <- which(!duplicated(LD_matrix))
  LD_matrix <- LD_matrix[dups,dups]
  ld.snps <- row.names(LD_matrix)

  # Remove duplicate SNPs
  finemap_DT <- finemap_DT[which(!duplicated(finemap_DT$SNP)),]
  fm.snps <- finemap_DT$SNP
  common.snps <- base::intersect(ld.snps, fm.snps)
  printer("+ LD_matrix =",length(ld.snps),"SNPs.", v=verbose)
  printer("+ finemap_DT =",length(fm.snps),"SNPs.", v=verbose)
  printer("+",length(common.snps),"SNPs in common.", v=verbose)
  # Subset/order LD matrix
  new_LD <- LD_matrix[common.snps, common.snps]
  new_LD[is.na(new_LD)] <- 0
  # Subset/order finemap_DT
  finemap_DT <- data.frame(finemap_DT)
  row.names(finemap_DT) <- finemap_DT$SNP
  new_DT <- data.table::as.data.table(finemap_DT[common.snps, ])
  new_DT <- unique(new_DT)
  # Reassign the lead SNP if it's missing
  if(sum(new_DT$leadSNP)==0){
    printer("+ leadSNP missing. Assigning new one by min p-value.", v=verbose)
    top.snp <- head(arrange(new_DT, P, desc(Effect)))[1,]$SNP
    new_DT$leadSNP <- ifelse(new_DT$SNP==top.snp,T,F)
  }
  # Check dimensions are correct
  if(nrow(new_DT)!=nrow(new_LD)){
    warning("+ LD_matrix and finemap_DT do NOT have the same number of SNPs.",v=verbose)
    warning("+ LD_matrix SNPs = ",nrow(new_LD),"; finemap_DT = ",nrow(finemap_DT), v=verbose)
  }
  printer("++ Subsetting complete.",v=verbose)
  return(list(LD=new_LD,
              DT=new_DT))
}

gene_trimmer <- function(subset_DT,
                         gene,
                         min_POS=NULL,
                         max_POS=NULL){
    printer("BIOMART:: Trimming data to only include SNPs within gene coordinates.")
    gene_info <- biomart_geneInfo(gene)
    gene_info_sub <- subset(gene_info, hgnc_symbol==gene)
    # Take most limiting min position
    min_POS <- max(min_POS, gene_info_sub$start_position, na.rm = T)
    # Take most limiting max position
    max_POS <- min(max_POS, gene_info_sub$end_position, na.rm = T)
    subset_DT <- subset(subset_DT, CHR==gene_info$chromosome_name[1] & POS>=min_POS & POS<=max_POS)
    printer("BIOMART::",nrow(subset_DT),"SNPs left after trimming.")
    return(subset_DT)
}

limit_SNPs <- function(max_snps=500, subset_DT){
  printer("echolocator:: Limiting to only",max_snps,"SNPs.")
  if(nrow(subset_DT)>max_snps){
    orig_n <- nrow(subset_DT)
    lead.index <- which(subset_DT$leadSNP==T)
    i=1
    tmp.sub<-data.frame()
    while(nrow(tmp.sub)<max_snps-1){
      # print(i)
      snp.start <- max(1, lead.index-i)
      snp.end <- min(nrow(subset_DT), lead.index+i)
      tmp.sub <- subset_DT[snp.start:snp.end]
      i=i+1
    }
    # Need to add that last row on only one end
    snp.start <- max(1, lead.index-i+1) # +1 keep it the same index as before
    snp.end <- min(nrow(subset_DT), lead.index+i)
    tmp.sub <- subset_DT[snp.start:snp.end]

  printer("+ Reduced number of SNPs:",orig_n,"==>",nrow(tmp.sub))
  return(tmp.sub)
  } else {
    printer("+ Data already contains less SNPs than limit (",nrow(subset_DT),"<",max_snps,")")
    return(subset_DT)
  }
}

printer <- function(..., v=T){if(v){print(paste(...))}}


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

echoR.filter_snps <- function(subset_DT,
                              bp_distance,
                              remove_variants,
                              gene,
                              verbose=T,
                              min_POS=NULL,
                              max_POS=NULL,
                              max_snps=NULL,
                              min_MAF=NULL,
                              trim_gene_limits=F){
  if(remove_variants!=F){
    printer("Removing specified variants:",paste(remove_variants, collapse=','), v=verbose)
    try({subset_DT <- subset(subset_DT, !(SNP %in% remove_variants) )})
  }
  # Trim subset according to annotations of where the gene's limit are
  if(trim_gene_limits){
    subset_DT <- gene_trimmer(subset_DT=subset_DT,
                              gene=gene,
                              min_POS=min_POS,
                              max_POS=min_POS)
  }
  if(!is.null(max_snps)){
    subset_DT <- limit_SNPs(max_snps = max_snps, subset_DT = subset_DT)
  }
  if(!is.null(min_MAF) & length(min_MAF>0)>0){
    printer("echolocatoR:: Removing SNPs w/ MAF <",min_MAF)
    subset_DT <- subset(subset_DT, MAF>=min_MAF)
  }
  # Limit range
  if(!is.null(bp_distance)){
    lead.snp <- subset(subset_DT, leadSNP)
    subset_DT <- subset(subset_DT,
                        POS >= lead.snp$POS - bp_distance &
                        POS <= lead.snp$POS + bp_distance)
  }
  if(!is.na(min_POS)){subset_DT <- subset(subset_DT, POS>=min_POS)}
  if(!is.na(max_POS)){subset_DT <- subset(subset_DT, POS<=max_POS)}
  printer("++ Post-filtered data:",paste(dim(subset_DT), collapse=" x "))
  return(subset_DT)
}


finemap_pipeline <- function(gene,
                             fullSS_path,
                             dataset_name, dataset_type="general",
                             top_SNPs="auto",
                             force_new_subset=F,
                             force_new_LD=F,
                             force_new_finemap=T,
                             finemap_methods=c("SUSIE","FINEMAP"),
                             bp_distance=500000,
                             n_causal=5,
                             sample_size=NA,
                             chrom_col="CHR",
                             position_col="POS",
                             snp_col="SNP",
                             pval_col="P",
                             effect_col="Effect",
                             stderr_col="StdErr",
                             tstat_col="t-stat",
                             gene_col="Gene",
                             freq_col="Freq",
                             MAF_col="MAF",
                             A1_col = "A1",
                             A2_col = "A2",
                             N_cases_col="N_cases",
                             N_controls_col="N_controls",
                             N_cases=NULL,
                             N_controls=NULL,
                             proportion_cases="calculate",

                             LD_reference="1KG_Phase1",
                             superpopulation="EUR",
                             download_reference=T,
                             min_POS=NA,
                             max_POS=NA,
                             min_MAF=NA,
                             trim_gene_limits=F,
                             max_snps=NULL,

                             file_sep="\t",
                             min_r2=0,
                             LD_block=F,
                             block_size=.7,
                             min_Dprime=F,
                             query_by="coordinates",
                             remove_variants=F,
                             remove_correlates=F,
                             probe_path = "./Data/eQTL/gene.ILMN.map",
                             conditioned_snps,
                             plot_LD = F,
                             verbose=T,
                             remove_tmps=T,
                             plot_types=c("simple"),
                             PAINTOR_QTL_datasets=NULL,
                             server=F,
                             PP_threshold=.95){
   # Create paths
   results_path <- make_results_path(dataset_name, dataset_type, gene)
   # Extract subset
   subset_DT <- extract_SNP_subset(gene = gene,
                                    top_SNPs = top_SNPs,
                                    fullSS_path = fullSS_path,
                                    results_path  =  results_path,
                                    force_new_subset = force_new_subset,

                                    chrom_col = chrom_col,
                                    position_col = position_col,
                                    snp_col = snp_col,
                                    pval_col = pval_col,
                                    effect_col = effect_col,
                                    stderr_col = stderr_col,
                                    gene_col = gene_col,
                                    tstat_col = tstat_col,
                                    MAF_col = MAF_col,
                                    freq_col = freq_col,
                                    A1_col = A1_col,
                                    A2_col = A2_col,

                                    N_cases_col = N_cases_col,
                                    N_controls_col = N_controls_col,
                                    N_cases = N_cases,
                                    N_controls = N_controls,
                                    proportion_cases = proportion_cases,

                                    bp_distance = bp_distance,
                                    superpopulation = superpopulation,
                                    min_POS = min_POS,
                                    max_POS = max_POS,
                                    file_sep = file_sep,
                                    query_by = query_by,
                                    probe_path = probe_path,
                                    remove_tmps = remove_tmps)

  ### Compute LD matrix
  message("--- Step 2: Calculate Linkage Disequilibrium ---")
  LD_matrix <- LD.load_or_create(results_path=results_path,
                                 subset_DT=subset_DT,
                                 gene=gene,
                                 force_new_LD=force_new_LD,
                                 LD_reference=LD_reference,
                                 superpopulation=superpopulation,
                                 download_reference=download_reference,
                                 min_r2=min_r2,
                                 LD_block=LD_block,
                                 block_size=block_size,
                                 min_Dprime=min_Dprime,
                                 remove_correlates=remove_correlates,
                                 verbose=verbose,
                                 server=server,
                                 remove_tmps=remove_tmps)

  #### ***** SNP Filters ***** ###
  # Remove pre-specified SNPs
  ## Do this step AFTER saving the LD to disk so that it's easier to re-subset in different ways later without having to redownload LD.
  message("-------------- Step 3: Filter SNPs -------------")
  subset_DT <- echoR.filter_snps(subset_DT=subset_DT,
                                  bp_distance=bp_distance,
                                  remove_variants=remove_variants,
                                  gene=gene,
                                  verbose=verbose,
                                  min_POS=min_POS,
                                  max_POS=max_POS,
                                  max_snps=max_snps,
                                  trim_gene_limits=trim_gene_limits)
  # Subset LD and df to only overlapping SNPs
  sub.out <- subset_common_snps(LD_matrix, subset_DT)
  LD_matrix <- sub.out$LD
  subset_DT <- sub.out$DT
  # finemap
  finemap_DT <- finemap_handler(results_path = results_path,
                                fullSS_path = fullSS_path,
                                finemap_methods = finemap_methods,
                                force_new_finemap = force_new_finemap,
                                dataset_type = dataset_type,
                                subset_DT = subset_DT,
                                LD_matrix = LD_matrix,
                                n_causal = n_causal,
                                sample_size = sample_size,
                                conditioned_snps = conditioned_snps,

                                snp_col = snp_col,
                                freq_col = freq_col,
                                effect_col = effect_col,
                                stderr_col = stderr_col,
                                pval_col = pval_col,
                                N_cases_col = N_cases_col,
                                N_controls_col = N_controls_col,
                                A1_col = A1_col,
                                A2_col = A2_col,
                                PAINTOR_QTL_datasets = PAINTOR_QTL_datasets,
                                PP_threshold = PP_threshold)
  finemap_DT <- find_consensus_SNPs(finemap_DT, credset_thresh = PP_threshold)
  # Step 6: COLOCALIZE
  # Step 7: Functionally Fine-map

  # Plot
  message("--------------- Step 7: Visualize --------------")
  if("simple" %in% plot_types){
    try({
      # mf_plot <- multi_finemap_plot(finemap_DT = finemap_DT,
      #                               LD_matrix = LD_matrix,
      #                               results_path = results_path,
      #                               finemap_method_list = finemap_methods,
      #                               conditioned_snps = conditioned_snps,
      #                               gene = gene,
      #                               original = T,
      #                               save = T)
      mf_plot <- GGBIO.plot(finemap_DT=finemap_DT,
                            LD_matrix=LD_matrix,
                            gene=gene,
                            results_path=results_path,
                            method_list=finemap_methods,
                            Nott_sn_epigenome = T,
                            XGR_libnames = NULL,
                            save_plot = T,
                            show_plot = T)
    })
  }
  if("fancy" %in% plot_types){
    try({
      trx <- GGBIO.plot(finemap_DT = finemap_DT,
                        gene = gene,
                        LD_matrix = LD_matrix,
                        results_path = results_path,
                        method_list = finemap_methods, #c("SUSIE","FINEMAP","PAINTOR","PAINTOR_Fairfax")
                        XGR_libnames = c("ENCODE_TFBS_ClusteredV3_CellTypes",
                                         "ENCODE_DNaseI_ClusteredV3_CellTypes"))
    })
  }

  # Plot LD
  if(plot_LD){
    try({
      LD_plot(LD_matrix=LD_matrix, subset_DT=subset_DT, span=10)
    })
  }




  # Show results table
  printer("+",gene,"Credible Set SNPs", v=verbose)
  print(createDT_html( subset(finemap_DT, Support >0) ))


  # Cleanup:
  if(remove_tmps){
    tmp_files <- file.path(results_path,"plink",
                           c("plink.bed",
                             "plink.bim",
                             "plink.fam",
                             "plink.ld",
                             "plink.ld.bin",
                             "plink.log",
                             "plink.nosex",
                             "SNPs.txt") )
    suppressWarnings(file.remove(tmp_files))
  }
  return(finemap_DT)
}

arg_list_handler <- function(arg, i){
  output <- if(length(arg)>1){arg[i]}else{arg}
  return(output)
}

snps_to_condition <- function(conditioned_snps, top_SNPs, loci){
  if(conditioned_snps=="auto"){
    lead_SNPs_DT <- subset(top_SNPs, Gene %in% loci)
    # Reorder
    lead_SNPs_DT[order(factor(lead_SNPs_DT$Gene,levels= loci)),]
    return(lead_SNPs_DT$SNP)
  } else {return(conditioned_snps)}
}

# Fine-ap iteratively over genes/loci
finemap_loci <- function(loci,
                         fullSS_path,
                         dataset_name,
                         dataset_type="general",
                         force_new_subset=F,
                         force_new_LD=F,
                         force_new_finemap=T,
                         subset_path="auto",
                         top_SNPs="auto",
                         finemap_methods=c("SUSIE","FINEMAP"),
                         bp_distance=500000,
                         n_causal=5,
                         sample_size=NA,
                         chrom_col="CHR",
                         position_col="POS",
                         snp_col="SNP",
                         pval_col="P",
                         effect_col="Effect",
                         stderr_col="StdErr",
                         tstat_col = "t-stat",
                         MAF_col="MAF",
                         gene_col="Gene",
                         freq_col="Freq",
                         A1_col = "A1",
                         A2_col = "A2",
                         N_cases_col="N_cases",
                         N_controls_col="N_controls",
                         N_cases=NULL,
                         N_controls=NULL,
                         proportion_cases="calculate",

                         LD_reference="1KG_Phase1",
                         superpopulation="EUR",
                         topVariants=3,
                         download_reference=T,
                         min_POS=NA,
                         max_POS=NA,
                         min_MAF=NA,
                         trim_gene_limits=F,
                         max_snps=NULL,
                         file_sep="\t",
                         min_r2=0, LD_block=F, block_size=.7, min_Dprime=F,
                         query_by="coordinates",
                         remove_variants=F,
                         remove_correlates=F,
                         probe_path = "./Data/eQTL/gene.ILMN.map",
                         conditioned_snps="auto",
                         plot_LD=F,
                         verbose=T,
                         remove_tmps=T,
                         plot_types = c("simple"),
                         PAINTOR_QTL_datasets=NULL,
                         server=F,
                         PP_threshold=.95){
  fineMapped_topSNPs <- data.table()
  fineMapped_allResults <- data.table()
  lead_SNPs <- snps_to_condition(conditioned_snps, top_SNPs, loci)

  for (i in 1:length(loci)){
    try({
      gene <- loci[i]
      message("🦇 🦇 🦇 ",gene," 🦇 🦇 🦇 ")
      lead_SNP <- arg_list_handler(lead_SNPs, i)
      gene_limits <- arg_list_handler(trim_gene_limits, i)
      conditioned_snp <- arg_list_handler(conditioned_snps, i)
      start_gene <- Sys.time()
      message("^^^^^^^^^ Running echolocatoR on: ",gene," ^^^^^^^^^")
      cat('  \n###', gene, '  \n')
      # Delete the old subset if force_new_subset == T
      finemap_DT <- finemap_pipeline(gene=gene,
                                     top_SNPs=top_SNPs,
                                     fullSS_path=fullSS_path,
                                     finemap_methods=finemap_methods,
                                     force_new_subset=force_new_subset,
                                     force_new_LD=force_new_LD,
                                     force_new_finemap=force_new_finemap,
                                     dataset_name=dataset_name,
                                     dataset_type=dataset_type,
                                     n_causal=n_causal,
                                     bp_distance=bp_distance,

                                     chrom_col=chrom_col,
                                     position_col=position_col,
                                     snp_col=snp_col,
                                     pval_col=pval_col,
                                     effect_col=effect_col,
                                     stderr_col=stderr_col,
                                     tstat_col=tstat_col,
                                     gene_col=gene_col,
                                     MAF_col=MAF_col,
                                     freq_col=freq_col,
                                     A1_col = A1_col,
                                     A2_col = A2_col,
                                     N_cases_col = N_cases_col,
                                     N_controls_col = N_controls_col,
                                     N_cases = N_cases,
                                     N_controls = N_controls,
                                     proportion_cases = proportion_cases,

                                     LD_reference=LD_reference,
                                     superpopulation=superpopulation,
                                     min_POS=min_POS,
                                     max_POS=max_POS,
                                     min_MAF=min_MAF,

                                     trim_gene_limits=gene_limits,
                                     max_snps=max_snps,
                                     file_sep=file_sep,
                                     min_r2=min_r2,
                                     LD_block=LD_block,
                                     block_size=block_size,
                                     min_Dprime=min_Dprime,
                                     query_by=query_by,
                                     remove_variants=remove_variants,
                                     remove_correlates=remove_correlates,
                                     probe_path=probe_path,
                                     conditioned_snps=lead_SNP,
                                     plot_LD=plot_LD,
                                     remove_tmps=remove_tmps,
                                     plot_types=plot_types,
                                     PAINTOR_QTL_datasets=PAINTOR_QTL_datasets,
                                     server=server,
                                     PP_threshold=PP_threshold)

      # Create summary table for all genes
      printer("Generating summary table...", v=verbose)
      newEntry <- cbind(data.table(Gene=gene), finemap_DT) %>% as.data.table()
      fineMapped_allResults <- rbind(fineMapped_allResults, newEntry, fill=T)
      end_gene <- Sys.time()
      printer(gene,"fine-mapped in", round(end_gene-start_gene, 2),"seconds", v=verbose)
      cat('  \n')
    })

  }
  return(fineMapped_allResults)
}




