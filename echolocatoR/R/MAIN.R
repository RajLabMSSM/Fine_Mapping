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

library(R.utils)
library(devtools)
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
library(BiocManager)
library(biomaRt) # BiocManager::install("biomaRt")
library(snpStats)  #BiocManager::install("snpStats") 
library(coloc)

# library(bigsnpr) # BiocManager::install("bigsnpr")
# install.packages("haploR", dependencies = TRUE)
library(haploR)
library(GeneOverlap) #BiocManager::install("GeneOverlap")


# *** SUSIE ****
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

source("./echolocatoR/R/finemap.R")
source("./echolocatoR/R/directory.R")
source("./echolocatoR/R/query.R")
source("./echolocatoR/R/LD.R")
source("./echolocatoR/R/plot.R")
source("./echolocatoR/R/conditional.R")
source("./echolocatoR/R/colocalization.R")
source("./echolocatoR/R/annotate.R")
source("./echolocatoR/R/eQTL_boxplots.R")




reload <- function(){
  source("echolocatoR/R/MAIN.R") 
}

# data.table::fwrite( subset(subset_DT, SNP != "rs34637584"), "Data/GWAS/Nalls23andMe_2019/LRRK2/LocusZoomData_-rs34637584.txt", sep="\t")

# quick_start
quickstart <- function(){ 
  reload()
  
  Data_dirs <<- read.csv("./Data/directories_table.csv")
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
  method <<- finemap_method <<- "SUSIE"
  method_list <<- finemap_method
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
  bp_distance <<- 500000
  n_causal <<- 5
  vcf_URL <<- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
  popDat_URL <<- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
  chr <<- 8
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
    topSS_path = Directory_info(Data_dirs, dataset_name, "topSumStats"),
    chrom_col = "CHR", position_col = "BP", snp_col="SNP",
    pval_col="P, all studies", effect_col="Beta, all studies", gene_col="Nearest Gene",
    caption= "Nalls et al. (2018) w/ 23andMe PD GWAS Summary Stats",
    group_by_locus = T, 
    locus_col = "Locus Number",
    remove_variants = "rs34637584")
  topSNP_sub <<- top_SNPs[top_SNPs$Gene==gene & !is.na(top_SNPs$Gene),]
  
  
  results_path <<- make_results_path(dataset_name, dataset_type, gene)
  subset_path <<- get_subset_path(results_path, gene)
  subset_DT <<- data.table::fread("Data/GWAS/Nalls23andMe_2019/LRRK2/Multi-finemap/Multi-finemap_results.txt", sep="\t")
  # subset_path <<- 'Data/GWAS/Nalls23andMe_2019/LRRK2/LRRK2_combined_meta_subset.txt'
  # subset_DT <<- data.table::fread(subset_path, sep="\t")
  LD_path <<- file.path(results_path, "plink/LD_matrix.RData")
  
  fullSS_path <<- file.path("Data",dataset_type,dataset_name,"nallsEtAl2019_allSamples_allVariants.mod.txt")
  
  

  colDict <<- column_dictionary(fullSS_path)
  if(is.na(min_POS)){min_POS <<- topSNP_sub$POS - bp_distance}
  if(is.na(max_POS)){max_POS <<- topSNP_sub$POS + bp_distance} 
  
  # COJO
  conditioned_snps <<- "rs76904798"#"rs34637584"
  excluded_snps <<- "rs34637584"
  min_MAF <<- 0
  GCTA_path <<- "echolocatoR/tools/gcta_1.92.1beta5_mac/bin/gcta64"
  bfiles <<- "plink_tmp/plink"
  load("Data/GWAS/Nalls23andMe_2019/LRRK2/plink/LD_matrix.RData")
  LD_matrix <<- LD_matrix
  # subset_DT <- subset(subset_DT, SNP %in% unique(row.names(LD_matrix), colnames(LD_matrix) ) )
  # subset_DT <- subset_DT[complete.cases(subset_DT),] # Remove any NAs
  # LD_matrix <- LD_matrix[row.names(LD_matrix) %in% subset_DT$SNP,  colnames(LD_matrix) %in% subset_DT$SNP]
  
  scaled_prior_variance <<- 0.1
  sample_size <<- NA 
  freq_cutoff <<- 0.1
  
  finemap_method_list <<- c("SUSIE", "ABF", "FINEMAP", "COJO")
  consensus <<- T
  
  dataset1_path <<- "./Data/GWAS/Nalls23andMe_2019/LRRK2/LRRK2_Nalls23andMe_2019_subset_500kb.txt"
  dataset2_path <<- "./Data/eQTL/MESA_CAU/LRRK2/LRRK2_MESA_CAU_subset.txt"
  shared_MAF <<- data.table::fread("Data/GWAS/Nalls23andMe_2019/LRRK2/LRRK2_Nalls23andMe_2019_subset_500kb.txt", sep="\t")$MAF
  plot_subtitle <<- "Fairfax (2014) + CD14 eQTL"
  dataset2_proportion_cases <<- 5e-324
  PP_threshold <<- 0.8
  force_new_finemap <<- F
  dataset1_type <<- "cc"
  dataset2_type <<- "quant"
  shared_MAF <<- 1
  which_merge <<-1
  show_plot <<- T
}
# quickstart()

# 
# quickstart_eQTL <- function(){
#   # Assign global variables to test functions
#   gene <<- "LRRK2"
#   leadSNP <<- "rs76904798"
#   gene_list <<- c("LRRK2")
#   chrom_col <<- "chr"
#   position_col <<- "pos_snps"
#   snp_col <<- "snps"
#   pval_col <<- "pvalue"
#   effect_col <<- "beta"
#   gene_col <<- "gene_name"
#   stderr_col <<- "calculate"
#   freq_col <<- "freq"
# 
#   fullSS_path <<- "Data/eQTL/MESA/CAU_cis_eqtl_summary_statistics.txt"
#   subset_path <<- "Data/eQTL/MESA/LRRK2_MESA_CAU_subset.txt"
#   vcf_folder <<- F#"../1000_Genomes_VCFs/Phase1"
#   superpopulation <<- "CAU"
#   force_new_subset <<- F
#   min_POS <<- NA
#   max_POS <<- NA
#   file_sep <<- "\t"
#   min_r2 <<- .2
#   LD_block <<- F
#   block_size <<- .7
#   min_Dprime <<- 1
#   plink_folder <<- "./plink_tmp"
#   reference <<- "1KG_Phase1"
#   LD_reference <<- "1KG_Phase1"
#   bp_distance <<- 500000
#   n_causal <<- 1
#   vcf_URL <<- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
#   popDat_URL <<- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
#   chr <<- 8
#   vcf_name <<- "ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi"
# }
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


gene_trimmer <- function(subset_DT, trim_gene_limits, gene, min_POS, max_POS){
  if(trim_gene_limits){
    gene_info <- biomart_geneInfo(gene)
    gene_info_sub <- subset(gene_info, hgnc_symbol==gene)
    # Take most limiting min position
    min_POS <- max(min_POS, gene_info_sub$start_position, na.rm = T) 
    # Take most limiting max position
    max_POS <- min(max_POS, gene_info_sub$end_position, na.rm = T) 
    return(subset(subset_DT, POS>=min_POS & POS<=max_POS))
  } else{return(subset_DT)} 
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
finemap_pipeline <- function(gene, 
                             fullSS_path, 
                             dataset_name, dataset_type="general", 
                             top_SNPs="auto", 
                             force_new_subset=F, 
                             force_new_LD=F,
                             force_new_finemap=T,
                             finemap_method="SUSIE",
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
                             N_cases_col="N_cases",
                             N_controls_col="N_controls",
                             A1_col = "A1",
                             A2_col = "A2",
                             LD_reference="1KG_Phase1", 
                             superpopulation="EUR", 
                             download_reference=T,
                             min_POS=NA, 
                             max_POS=NA,
                             min_MAF=NA,
                             trim_gene_limits=T,
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
                             remove_tmps=T
                          ){
   # Create paths 
   results_path <- make_results_path(dataset_name, dataset_type, gene)
   subset_path <- get_subset_path(results_path=results_path, gene=gene, subset_path="auto")
   # delete_subset(force_new_subset, subset_path) 
   
   # Extract subset 
   subset_DT <- extract_SNP_subset(gene = gene, 
                      top_SNPs = top_SNPs, 
                      fullSS_path = fullSS_path,
                      subset_path  =  subset_path,
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
                      
                      bp_distance = bp_distance,
                      superpopulation = superpopulation,  
                      min_POS = min_POS, 
                      max_POS = max_POS, 
                      file_sep = file_sep, 
                      query_by = query_by,
                      probe_path = probe_path,
                      remove_tmps = remove_tmps
                      ) 
   # Remove pre-specified SNPs
   if(remove_variants!=F){
     printer("Removing specified variants:",paste(remove_variants, collapse=','), v=verbose)
     subset_DT <- subset(subset_DT, !SNP %in% remove_variants )
   }
   # Filter by MAF
   if(!is.na(min_MAF)){ 
     subset_DT <- subset(subset_DT, MAF >= min_MAF) 
   } 
   # Trim subset according to annotations of where the gene's limit are 
   subset_DT <- gene_trimmer(subset_DT, trim_gene_limits, gene, min_POS, max_POS) 
   
  ### Compute LD matrix  
  printer("",v=verbose)
  printer("--- Step 2: Calculate Linkage Disequilibrium ---", v=verbose)
  
  LD_path <- file.path(results_path,"plink/LD_matrix.RData")
  if(!file.exists(LD_path) | force_new_LD==T){
    printer("+ Computing LD matrix... \n", verbose) 
    LD_matrix <- compute_LD_matrix(results_path = results_path, 
                                   subset_DT = subset_DT, 
                                   gene = gene,
                                   reference = LD_reference,
                                   superpopulation = superpopulation, 
                                   download_reference = download_reference,
                                   
                                   min_r2 = min_r2,
                                   LD_block = LD_block,
                                   block_size = block_size,
                                   min_Dprime = min_Dprime,
                                   remove_correlates = remove_correlates) 
    # Save LD matrix 
    # data.table::fwrite(LD_matrix, LD_path, sep="\t") 
    printer("+ Saving LD matrix to:",LD_path, v=verbose) 
    save(LD_matrix, file = LD_path) 
    # write.table(LD_matrix, LD_path, sep="\t", quote = F) 
  } else { 
    printer("+ Previously computed LD matrix detected. Importing...",LD_path, v=verbose) 
    # LD_matrix <- data.table::fread(LD_path, sep="\t", stringsAsFactors = F)  
    load(LD_path)
  }
  # Plot LD 
  if(plot_LD){
    try({ 
      LD_plot(LD_matrix=LD_matrix, subset_DT=subset_DT, span=10)
    })
  }
  
  
  # Final filtering  
  printer("",v=verbose)
  printer("-------------- Step 3: Filter SNPs -------------",v=verbose) 
  ## Subset summary stats to only include SNPs found in query
  subset_DT <- subset(subset_DT, SNP %in% unique(row.names(LD_matrix), colnames(LD_matrix) ) )
  subset_DT <- subset_DT[complete.cases(subset_DT),] # Remove any NAs
  LD_matrix <- LD_matrix[row.names(LD_matrix) %in% subset_DT$SNP,  colnames(LD_matrix) %in% subset_DT$SNP]
  
  # finemap 
  finemap_DT <- finemap_handler(results_path = results_path, 
                                fullSS_path = fullSS_path,
                                finemap_method = finemap_method, 
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
                                A2_col = A2_col)  
  # Step 6: COLOCALIZE
  # Step 7: Functionally Fine-map
  
  # Plot   
  printer("",v=verbose)
  printer("--------------- Step 7: Visualize --------------", v=verbose) 
  mf_plot <- multi_finemap_plot(finemap_DT = finemap_DT,
                     results_path = results_path,
                     finemap_method_list = finemap_method, 
                     conditioned_snps = conditioned_snps,
                     gene = gene, 
                     original = T, 
                     save = T)
  print(mf_plot)  
  
  # Show results table 
  printer("+",gene,"Credible Set SNPs", v=verbose)
  createDT( subset(finemap_DT, Support >0) ) %>% print() 
  
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
snps_to_condition <- function(conditioned_snps, top_SNPs, gene_list){ 
  if(conditioned_snps=="auto"){
    lead_SNPs_DT <- subset(top_SNPs, Gene %in% gene_list)
    # Reorder
    lead_SNPs_DT[order(factor(lead_SNPs_DT$Gene,levels= gene_list)),] 
    return(lead_SNPs_DT$SNP)
  } else {return(conditioned_snps)}
}

# Fine-ap iteratively over genes/loci
finemap_gene_list <- function(gene_list, fullSS_path, 
                             dataset_name,
                             dataset_type="general",
                             force_new_subset=F, 
                             force_new_LD=F,
                             force_new_finemap=T,
                             subset_path="auto", 
                             top_SNPs="auto", 
                             finemap_method="SUSIE",
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
                             LD_reference="1KG_Phase1", 
                             superpopulation="EUR",
                             topVariants=3, 
                             download_reference=T, 
                             min_POS=NA, 
                             max_POS=NA,
                             min_MAF=NA,
                             trim_gene_limits=T,
                             file_sep="\t",
                             min_r2=0, LD_block=F, block_size=.7, min_Dprime=F, 
                             query_by="coordinates", 
                             remove_variants=F,
                             remove_correlates=F,
                             probe_path = "./Data/eQTL/gene.ILMN.map",
                             conditioned_snps="auto",
                             plot_LD=F,
                             verbose=T,
                             remove_tmps=T
                             ){ 
  fineMapped_topSNPs <- data.table()
  fineMapped_allResults <- data.table()
  lead_SNPs <- snps_to_condition(conditioned_snps, top_SNPs, gene_list)
  
  for (i in 1:length(gene_list)){
    try({ 
      gene <- gene_list[i]
      lead_SNP <- arg_list_handler(lead_SNPs, i) 
      gene_limits <- arg_list_handler(trim_gene_limits, i) 
      conditioned_snp <- arg_list_handler(conditioned_snps, i)  
      start_gene <- Sys.time()  
      cat('  \n###', gene, '  \n')  
      # Delete the old subset if force_new_subset == T  
      finemap_DT <- finemap_pipeline(gene=gene, 
                                     top_SNPs=top_SNPs, 
                                     fullSS_path=fullSS_path,
                                     finemap_method=finemap_method,
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
                                     
                                     LD_reference=LD_reference, 
                                     superpopulation=superpopulation,
                                     min_POS=min_POS, 
                                     max_POS=max_POS,
                                     min_MAF=min_MAF,
                                     
                                     trim_gene_limits=gene_limits,
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
                                     remove_tmps=remove_tmps)  
      
      # Create summary table for all genes
      printer("Generating summary table...", v=verbose)
      newEntry <- cbind(data.table(Gene=gene), finemap_DT) %>% as.data.table() 
      fineMapped_allResults <- rbind(fineMapped_allResults, newEntry) 
      end_gene <- Sys.time()
      printer(gene,"fine-mapped in", round(end_gene-start_gene, 2),"seconds", v=verbose)
    }) 
    cat('  \n') 
  }  
  return(fineMapped_allResults)
}
 



