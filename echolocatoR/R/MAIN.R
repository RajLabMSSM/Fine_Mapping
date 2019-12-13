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
 
# Load libraries
# .libPaths()  

library(R.utils)
library(devtools)
library(readxl)
# library(DT)
library(data.table)
library(dplyr)
library(ggplot2)
library(plotly)
library(patchwork) #devtools::install_github("thomasp85/patchwork")
# library(cowplot)
library(ggrepel)
library(curl) 
library(gaston)
library(tidyr)
# library(BiocManager)
library(biomaRt) # BiocManager::install("biomaRt") 
# library(refGenome)

library(crayon) 

# install_fGWAS <- function(){
#   devtools::install_github("wzhy2000/fGWAS/pkg")
#   system("git clone https://github.com/wzhy2000/fGWAS.git & cd fGWAS & R CMD INSTALL pkg")
# }
# library(fGWAS) 
# library(snpStats) #BiocManager::install("snpStats") 
# library(coloc)


# library(bigsnpr) # BiocManager::install("bigsnpr")
# install.packages("haploR", dependencies = TRUE)
# library(haploR)
# library(GeneOverlap) #BiocManager::install("GeneOverlap")
# library(rtracklayer) #BiocManager::install("rtracklayer")
# BiocManager::install(c("supraHex","graph","Rgraphviz","dnet")) 
# library(XGR)# install.packages("XGR")



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
source("./echolocatoR/R/query.R")
source("./echolocatoR/R/LD.R")
# Fine-mapping
source("./echolocatoR/R/Finemapping/multi_finemap.R")
source("./echolocatoR/R/Finemapping/SUSIE.R")
source("./echolocatoR/R/Finemapping/ABF.R")
source("./echolocatoR/R/Finemapping/FINEMAP.R") 
source("./echolocatoR/R/Finemapping/PAINTOR.R")
source("./echolocatoR/R/Finemapping/COJO.R")
source("./echolocatoR/R/Finemapping/COLOC.R")
source("./echolocatoR/R/Finemapping/fGWAS.R")
source("./echolocatoR/R/Finemapping/POLYFUN.R")

# Plotting
source("./echolocatoR/R/Plotting/plot.R")
source("./echolocatoR/R/Plotting/ggbio.R")
source("./echolocatoR/R/Plotting/QTL_boxplots.R")
# Annotation
source("./echolocatoR/R/Annotation/annotate.R")
source("./echolocatoR/R/Annotation/GoShifter.R") # ***
source("./echolocatoR/R/Annotation/XGR.R")
source("./echolocatoR/R/Annotation/psychENCODE.R")
source("./echolocatoR/R/Annotation/mergeQTL.R")



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
startup_image()


reload <- function(){
  source("echolocatoR/R/MAIN.R") 
}


rbind.file.list <- function(file.list, verbose=T){
  merged.dat <- lapply(file.list, function(x){
    printer(x, v = verbose)
    dat <- data.table::fread(x, nThread = 4) 
    return(dat)
  }) %>% data.table::rbindlist()
  return(merged.dat)
}

# data.table::fwrite( subset(subset_DT, SNP != "rs34637584"), "Data/GWAS/Nalls23andMe_2019/LRRK2/LocusZoomData_-rs34637584.txt", sep="\t")

quick_finemap <- function(){
  gene <- "LRRK2"
  results_path <<- file.path("./Data/GWAS/Nalls23andMe_2019",gene)
  finemap_DT <<- data.table::fread(file.path(results_path, "Multi-finemap/Multi-finemap_results.txt"))
  # return(finemap_DT)
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
# quickstart()

# 
# quickstart_eQTL <- function(){
#   # Assign global variables to test functions
#   gene <<- "LRRK2"
#   leadSNP <<- "rs76904798"
#   loci <<- c("LRRK2")
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


dt.replace <- function(DT, target, replacement){
  # Annoyingly, there is no native function for this in data.table
  for(col in names(DT)) set(DT, i=which(DT[[col]]==target), j=col, value=replacement)
  return(DT)
}


subset_common_snps <- function(LD_matrix, finemap_DT){
  ld.snps <- unique(row.names(LD_matrix))
  fm.snps <- unique(finemap_DT$SNP)
  common.snps <- intersect(ld.snps, fm.snps)
  printer("+ LD_matrix =",length(ld.snps),"SNPs.")
  printer("+ finemap_DT =",length(fm.snps),"SNPs.")
  printer("+",length(common.snps),"SNPs in common.")
  LD_matrix <- LD_matrix[common.snps, common.snps] 
  finemap_DT <- subset(finemap_DT, SNP %in% common.snps)
  return(list(LD_matrix=LD_matrix,
              finemap_DT=finemap_DT))
}

gene_trimmer <- function(subset_DT, 
                         gene, 
                         trim_gene_limits=F, 
                         min_POS=NULL, 
                         max_POS=NULL){
  if(trim_gene_limits){
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
   
  } else{return(subset_DT)} 
}

limit_SNPs <- function(max_snps=500, subset_DT){
  printer("echolocator:: Limiting to only",max_snps,"SNPs.")
  if(nrow(subset_DT)>max_snps){
    orig_n <- nrow(subset_DT) 
    lead.index <- which(subset_DT$leadSNP==T) 
    i=1
    tmp.sub<-data.frame() 
    while(nrow(tmp.sub)<max_snps){
      # print(i)
      snp.start <- max(1, lead.index-i)
      snp.end <- min(nrow(subset_DT), lead.index+i)
      tmp.sub <- subset_DT[snp.start:snp.end]
      i=i+1
    }
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
                             plot_types=c("simple","fancy"),
                             PAINTOR_QTL_datasets=NULL,
                             server=F){
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
                      remove_tmps = remove_tmps) 
   # Remove pre-specified SNPs
   if(remove_variants!=F){
     printer("Removing specified variants:",paste(remove_variants, collapse=','), v=verbose)
     try({subset_DT <- subset(subset_DT, !(SNP %in% remove_variants) )}) 
   }
   # Filter by MAF
   if(!is.na(min_MAF)){ 
     try({subset_DT <- subset(subset_DT, MAF >= min_MAF) })
   } 
   # Trim subset according to annotations of where the gene's limit are 
   subset_DT <- gene_trimmer(subset_DT, gene=gene,
                             trim_gene_limits=trim_gene_limits,
                             min_POS=min_POS,
                             max_POS=min_POS) 
   
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
                                 server=server)
  # Subset LD and df to only overlapping SNPs 
  sub.out <- subset_common_snps(LD_matrix, subset_DT)
  LD_matrix <- sub.out$LD_matrix
  subset_DT <- sub.out$finemap_DT 
  
  if(!is.null(max_snps)){
    subset_DT <- limit_SNPs(max_snps = max_snps, subset_DT = subset_DT)
    sub.out <- subset_common_snps(LD_matrix, subset_DT)
    LD_matrix <- sub.out$LD_matrix
    subset_DT <- sub.out$finemap_DT
  }
   
  # Plot LD 
  if(plot_LD){
    try({ 
      LD_plot(LD_matrix=LD_matrix, subset_DT=subset_DT, span=10)
    })
  }
  
  
  # Final filtering   
  message("-------------- Step 3: Filter SNPs -------------") 
  ## Subset summary stats to only include SNPs found in query
  # subset_DT <- subset(subset_DT, SNP %in% unique(row.names(LD_matrix), colnames(LD_matrix) ) )
  # subset_DT <- subset_DT[complete.cases(subset_DT),] # Remove any NAs
  # LD_matrix <- LD_matrix[row.names(LD_matrix) %in% subset_DT$SNP,  colnames(LD_matrix) %in% subset_DT$SNP]
  
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
                                PAINTOR_QTL_datasets = PAINTOR_QTL_datasets)  
  finemap_DT <- find_consensus_SNPs(finemap_DT)
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
      mf_plot <- ggbio_plot(finemap_DT=finemap_DT, 
                            LD_matrix=LD_matrix, 
                            gene=gene, 
                            results_path=results_path, 
                            method_list=finemap_methods,
                            XGR_libnames = NULL)
      print(mf_plot)
    }) 
  }
  if("fancy" %in% plot_types){
    try({
      trx <- ggbio_plot(finemap_DT = finemap_DT,
                        gene = gene,
                        LD_matrix = LD_matrix,
                        results_path = results_path,
                        method_list = finemap_methods, #c("SUSIE","FINEMAP","PAINTOR","PAINTOR_Fairfax")
                        XGR_libnames = c("ENCODE_TFBS_ClusteredV3_CellTypes",
                                         "ENCODE_DNaseI_ClusteredV3_CellTypes")) 
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
finemap_loci <- function(loci, fullSS_path, 
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
                             plot_types = c("simple","fancy"),
                             PAINTOR_QTL_datasets=NULL,
                             server=F){ 
  fineMapped_topSNPs <- data.table()
  fineMapped_allResults <- data.table()
  lead_SNPs <- snps_to_condition(conditioned_snps, top_SNPs, loci)
  
  for (i in 1:length(loci)){
    try({ 
      gene <- loci[i]
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
                                     server=server)  
      
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
 



