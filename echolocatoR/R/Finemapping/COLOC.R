#     %%%%%%%%%%%%%%%%%      #
####### Colocalization ####### 
#     %%%%%%%%%%%%%%%%%     #


# Jansen et al. 2017: https://eqtl.onderzoek.io/index.php?page=gene_cis_details&gene=BST1

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

COLOC.construct_dataset <- function(subset_DT, 
                                    sample_size=NA, 
                                    proportion_cases=5e-324, # Doesn't allow actual 0, so use smallest number R allows
                                    MAF=NA,
                                    type="cc"){ 
  sample_size <- get_sample_size(subset_DT, sample_size)
  if("proportion_cases" %in% colnames(subset_DT)){
    printer("++ Extracting Proportion of Cases...")
    proportion_cases <- getmode(subset_DT$proportion_cases)
  }
  if(length(MAF)==1){
    printer("++ Extracting MAF...")
    if(is.na(MAF) & "MAF" %in% colnames(subset_DT)){
      MAF <- subset_DT$MAF
    }  
  }
  # List form
  dataset <- list(pvalues = subset_DT$P, 
                  beta = subset_DT$Effect,
                  varbeta = subset_DT$StdErr^2, # MUST be squared
                  snp = subset_DT$SNP,
                  
                  N = sample_size, # [optional]
                  s = proportion_cases, # use overall proportions
                  MAF = MAF, # [required]
                  type = type)
  return(dataset)
}

 

COLOC <- function(gene,
                  # dataset1_path, 
                  # dataset2_path,
                  subset_DT1,
                  subset_DT2,
                  shared_MAF=1, 
                  dataset2_proportion_cases=5e-324, 
                  PP_threshold=0.8,
                  save_results=T,
                  show_plot=T){
  printer("\n******** Step 5: COLOCALIZE ********\n")  
  # The Approximate Bayes Factor colocalisation analysis described in the next section 
  ## essentially works by fine mapping each trait under a single causal variant assumption 
  ##and then integrating over those two posterior distributions to calculate probabilities that 
  ## those variants are shared.
  # https://cran.r-project.org/web/packages/coloc/vignettes/vignette.html
  # subset_DT1 <- data.table::fread(dataset1_path, stringsAsFactors = F, sep="\t") %>% data.frame()
  # subset_DT2 <- data.table::fread(dataset2_path, stringsAsFactors = F, sep="\t") %>% data.frame()
  common_SNPs <- intersect(subset_DT1$SNP, subset_DT2$SNP)
  subset_DT1 <- subset(subset_DT1, SNP %in% common_SNPs) %>% group_by(SNP) %>% slice(1) %>% data.frame()
  subset_DT2 <- subset(subset_DT2, SNP %in% common_SNPs) %>% group_by(SNP) %>% slice(1) %>% data.frame()

  dataset1 <- COLOC.construct_dataset(subset_DT1,
                                      type = dataset1_type)
  dataset2 <- COLOC.construct_dataset(subset_DT2, 
                                      proportion_cases = dataset2_proportion_cases,
                                      type = dataset2_type)
  if(shared_MAF==1){
    dataset2$MAF <- dataset1$MAF
  } else if(shared_MAF==2){
    dataset1$MAF <- dataset2$MAF
  }
  ## NOTES: MESA and Fairfax: No sample size (SNP-level), proportion of cases, or freq/MAF info available?   
  printer("\n\n")
  coloc.res <- coloc::coloc.abf(dataset1 = dataset1,
                                    dataset2 = dataset2)
                                    # MAF = dataset1$MAF) 
 hypothesis_key <- setNames(
   c("Neither trait has a genetic association in the region.",
     "Only trait 1 has a genetic association in the region.",
     "Only trait 2 has a genetic association in the region.",
     "Both traits are associated, but with different causal variants (one in each dataset).",
     "Both traits are associated and share a single causal variant.") ,
   c("PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf")) 
 # Report hypothess results
  printer("\n Hypothesis Results @ PP_threshold =",PP_threshold,":")
  true_hyp <-""
  for(h in names(hypothesis_key)){
    if(coloc.res$summary[h]>=PP_threshold){
      hyp <- hypothesis_key[h]
      printer("\n    ",h,"== TRUE: **",hyp )
      true_hyp <- paste0(names(hyp),": ", hyp)
    } else{
      printer("\n    ",h,"== FALSE: ")
    } 
  } 
  
  # Save raw results   
  coloc_DT <- coloc.res$results
  results_path <- dirname(dirname(dataset1_path))
  data.table::fwrite(coloc_DT, file.path(results_path, "COLOC/COLOC_raw.txt"), sep="\t")
  
  # Find the causal SNP that coloc.abf identified in each dataset via finemap.abf
  # DT1
  causal_DT1 <- coloc_DT %>% arrange(desc(lABF.df1))
  causal_DT1 <- causal_DT1$snp[1]
  # DT2
  causal_DT2 <- coloc_DT %>% arrange(desc(lABF.df2))
  causal_DT2 <- causal_DT2$snp[1]
  
  # Process results  
  coloc_DT$Colocalized <- ifelse(coloc_DT$SNP.PP.H4 >= PP_threshold, T, F)
  colocalized_snps <- subset(coloc_DT, Colocalized==T)$snp# subset(coloc_DT, Colocalized==1)$SNP
  coloc_datasets <- coloc_plot_data(coloc.res, subset_DT1, subset_DT2)
  subtitle2 <- paste0("Colocalized SNPs: ", paste(colocalized_snps,sep=", "))
  if((coloc.res$summary["PP.H3.abf"] + coloc.res$summary["PP.H4.abf"] >= PP_threshold) & 
     (coloc.res$summary["PP.H4.abf"]/coloc.res$summary["PP.H3.abf"] >= 2)){
    # "We called the signals colocalized when (coloc H3+H4 ≥ 0.8 and H4∕H3 ≥ 2)" -Yi et al. (2019)
    report <- paste("Datasets colocalized")  
  } else {
    report <- paste("Datasets NOT colocalized") 
  }   
  printer("\n++",report,"at PP.H3 + PP.H4 >=",PP_threshold," and PP.H3 / PP.H4 >= 2.","\n") 
 
  
  # Plot 
  if(show_plot){
    title1 <- paste(gene,":",strsplit(dataset1_path,"/")[[1]][4],strsplit(dataset1_path,"/")[[1]][3])
    title2 <- paste(gene,":",strsplit(dataset2_path,"/")[[1]][4],strsplit(dataset2_path,"/")[[1]][3]) 
    
    COLOC.plot(coloc_DT1 = coloc_datasets$coloc_DT1,
               coloc_DT2 = coloc_datasets$coloc_DT2, 
               title1 = title1,
               subtitle1 = report,
               title2 = title2, 
               subtitle2 = subtitle2,
               SNP_list = c("rs34637584","rs76904798","rs117073808"),
               alt_color_SNPs = colocalized_snps, 
               show_plot = T)
  } 
  # Save
  if(save_results){
    coloc_path <- file.path(dirname(dataset1_path),"COLOC")
    dir.create(coloc_path, recursive = T, showWarnings = F)
    data.table::fwrite(coloc_DT, file.path(coloc_path,"COLOC_results.txt"), sep="\t")
  }
  return(coloc_DT)
} 



COLOC.plot_data <- function(coloc.res, 
                            subset_DT1, 
                            subset_DT2){
  coloc_DT <- coloc.res$results
  # Merge Dataset 1
  coloc_DT1 <- coloc_DT %>% dplyr::select(SNP="snp", 
                             V = "V.df1", 
                             Z= "z.df1", 
                             lABF = "lABF.df1", 
                             H4.Probability = "SNP.PP.H4")
  coloc_DT1 <- data.table:::merge.data.table(subset_DT1, coloc_DT1, on = "SNP", all = F)
  # Merge Dataset 2
  coloc_DT2 <- coloc_DT %>% dplyr::select(SNP="snp", 
                                          V = "V.df2", 
                                          Z= "z.df2", 
                                          lABF = "lABF.df2", 
                                          H4.Probability = "SNP.PP.H4")
  coloc_DT2 <- data.table:::merge.data.table(subset_DT2, coloc_DT2, on = "SNP", all = F) 
  return(list(coloc_DT1=coloc_DT1, coloc_DT2=coloc_DT2))
}

COLOC.plot <- function(gene, 
                       coloc_DT1, 
                       coloc_DT2, 
                       title1, 
                       title2, 
                       subtitle1,
                       subtitle2,
                       SNP_list=c(), 
                       alt_color_SNPs=c(), 
                       show_plot=T){ 
  mp1 <- manhattan_plot(subset_DT = coloc_DT1, 
                 gene = gene, 
                 SNP_list = SNP_list,
                 alt_color_SNPs = alt_color_SNPs,
                 title = title1,
                 subtitle = subtitle1, 
                 show_plot = F)
  mp2 <- manhattan_plot(subset_DT = coloc_DT2, 
                       gene = gene, 
                       SNP_list = SNP_list,
                       alt_color_SNPs = alt_color_SNPs,
                       title = title2,
                       subtitle = subtitle2, 
                       show_plot = F)
  cp <- cowplot::plot_grid(mp1, mp2, ncol = 1)
  if(show_plot){print(cp)}else{return(cp)}
}














COLOC.report_summary <- function(coloc.res, PP_threshold=.8){ 
  # MAF = dataset1$MAF) 
  hypothesis_key <- setNames(
    c("Neither trait has a genetic association in the region.",
      "Only trait 1 has a genetic association in the region.",
      "Only trait 2 has a genetic association in the region.",
      "Both traits are associated, but with different causal variants (one in each dataset).",
      "Both traits are associated and share a single causal variant.") ,
    c("PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf")) 
  # Report hypothess results
  printer("\n Hypothesis Results @ PP_threshold =",PP_threshold,":")
  true_hyp <-""
  for(h in names(hypothesis_key)){
    if(!is.na(coloc.res$summary[h])){
      if(coloc.res$summary[h]>=PP_threshold){
        hyp <- hypothesis_key[h]
        printer("    ",h,"== TRUE: **",hyp )
        true_hyp <- paste0(names(hyp),": ", hyp)
      }
    } else{
      printer("    ",h,"== FALSE: ")
    } 
  } 
  
  # Save raw results   
  coloc_DT <- coloc.res$results
  # Process results  
  coloc_DT$Colocalized <- ifelse(coloc_DT$SNP.PP.H4 >= PP_threshold, T, F)
  colocalized_snps <- subset(coloc_DT, Colocalized==T)$snp# subset(coloc_DT, Colocalized==1)$SNP
  subtitle2 <- paste0("Colocalized SNPs: ", paste(colocalized_snps,sep=", "))
  if(!is.na(coloc.res$summary)["PP.H4.abf"] ){
    if((coloc.res$summary["PP.H3.abf"] + coloc.res$summary["PP.H4.abf"] >= PP_threshold) & 
       (coloc.res$summary["PP.H4.abf"]/coloc.res$summary["PP.H3.abf"] >= 2)){
      # "We called the signals colocalized when (coloc H3+H4 ≥ 0.8 and H4∕H3 ≥ 2)" -Yi et al. (2019)
      report <- paste("Datasets colocalized")  
    } 
  } else {
    report <- paste("Datasets NOT colocalized") 
  }   
  printer("+ COLOC::",report,"at: PP.H3 + PP.H4 >=",PP_threshold," and PP.H3 / PP.H4 >= 2.") 
  return(coloc_DT)
}


COLOC.plot_results <- function(){
  # Find the causal SNP that coloc.abf identified in each dataset via finemap.abf
  # DT1
  causal_DT1 <- coloc_DT %>% arrange(desc(lABF.df1))
  causal_DT1 <- causal_DT1$snp[1]
  # DT2
  causal_DT2 <- coloc_DT %>% arrange(desc(lABF.df2))
  causal_DT2 <- causal_DT2$snp[1]
  coloc_datasets <- coloc_plot_data(coloc.res, subset_DT1, subset_DT2)
  # Plot 
  title1 <- paste(gene,":",strsplit(dataset1_path,"/")[[1]][4],strsplit(dataset1_path,"/")[[1]][3])
  title2 <- paste(gene,":",strsplit(dataset2_path,"/")[[1]][4],strsplit(dataset2_path,"/")[[1]][3]) 
  
  COLOC.plot(coloc_DT1 = coloc_datasets$coloc_DT1,
             coloc_DT2 = coloc_datasets$coloc_DT2, 
             title1 = title1,
             subtitle1 = report,
             title2 = title2, 
             subtitle2 = subtitle2,
             SNP_list = c("rs34637584","rs76904798","rs117073808"),
             alt_color_SNPs = colocalized_snps, 
             show_plot = T) 
}


COLOC.PP4_plot <- function(COLOC_DT){
  COLOC_DT$coloc <- (COLOC_DT$PP.H3.abf + COLOC_DT$PP.H4.abf >= PP_threshold) &  (COLOC_DT$PP.H4.abf/COLOC_DT$PP.H3.abf >= 2)
  COLOC_DT$Tissue <- gsub(paste0(GTEx_version,"_"),"",COLOC_DT$Dataset2)
  
  cp <- ggplot(subset(COLOC_DT, coloc==T), aes(x=Tissue, y=PP.H4.abf, fill=Tissue)) + 
    facet_grid(~Locus) + 
    geom_col(show.legend = F) + 
    coord_flip() + 
    theme_bw()  
  cp
}



COLOC.iterate_GTEx <- function(GTEx_version="GTEx_V7"){
  FM_all <- merge_finemapping_results(minimum_support = 0, 
                                      include_leadSNPs = T, 
                                      dataset = "./Data/GWAS/Nalls23andMe_2019")
  QTL_files <- list.files(path = "./Data/QTL", pattern = "*.finemap.txt.gz", recursive = T)
  QTL_files <- grep(paste(c("GTEx*","MESA*"),collapse="|"), QTL_files, value = T)
 
  # qtl_file = QTL_files[28] # Need allele/MAF info: 1:4 (Brain_xQTL_Serve), 5:6 (Cardiogenics), 7:10 (Fairfax), 28:31 (psychENCODE) ##### (GTEx and MESA are good (tho MESA needs sample size))
  COLOC_DT <- lapply(QTL_files, function(qtl_file){
    dataset_name <- gsub(".finemap.txt.gz","",basename(qtl_file))
    FM_merge <- mergeQTL.merge_handler(FM_all = FM_all, qtl_file = qtl_file)
    FM_merge$mean.PP <- rowMeans(subset(FM_merge, select=grep(".Probability",colnames(FM_merge))))
    FM_merge$Adjusted.Effect <- FM_merge$Effect * FM_merge$mean.PP
    
    coloc_dt <- lapply(unique(FM_merge$Gene), function(gene){
          printer("+COLOC::",gene) 
          FM_gene <- subset(FM_merge, Gene==gene) 
          FM_gene <- FM_gene[!is.na(FM_gene$QTL.Effect),]  
          default_results <- data.table::data.table(Locus=gene,
                                                    Dataset1=dataset, 
                                                    Dataset2=dataset_name, 
                                                    nSNPs_locus=nrow(FM_merge),
                                                    nSNPs_overlap=nrow(FM_gene),
                                                    PP.H0.abf=NA,
                                                    PP.H1.abf=NA,
                                                    PP.H2.abf=NA,
                                                    PP.H3.abf=NA,
                                                    PP.H4.abf=NA)
          if(nrow(FM_gene)>0){
            # GWAS
            dataset.gwas <- list(pvalues = FM_gene$P, 
                                 beta = FM_gene$mean.PP,
                                 varbeta = FM_gene$StdErr^2, # MUST be squared
                                 snp = FM_gene$SNP,
                                 
                                 N =  max(FM_gene$N_cases, na.rm = T) + max(FM_gene$N_controls, na.rm = T), # [optional]
                                 s = getmode(FM_gene$proportion_cases), # use overall proportions
                                 MAF = FM_gene$MAF, # [required]
                                 type = "cc")
            # QTL
            if(all(is.na(FM_gene$QTL.MAF))){FM_gene$QTL.MAF <- FM_gene$MAF; printer("+ COLOC:: No QTL.MAF given. Using GWAS.MAF instead.")} 
            dataset.qtl <- list(pvalues = FM_gene$QTL.P, 
                                beta = FM_gene$QTL.Effect,
                                varbeta = FM_gene$QTL.StdErr^2, # MUST be squared
                                snp = FM_gene$SNP,
                                
                                N = max(FM_gene$QTL.SampleSize), # [optional]
                                # s = NA, # use overall proportions
                                MAF = FM_gene$QTL.MAF, # [required]
                                type = "quant")
            
            coloc.res <- coloc::coloc.abf(dataset1 = dataset.gwas,
                                          dataset2 = dataset.qtl)
            COLOC.report_summary(coloc.res, PP_threshold = .8) 
            # dat <- data.table::data.table(coloc.res$results)
            dat <- data.table::data.table(t(coloc.res$summary))
            dat <- cbind(Locus=gene, 
                         Dataset1=dataset, 
                         Dataset2=dataset_name, 
                         nSNPs_locus=n_SNPs,
                         dat) %>% dplyr::rename(nSNPs_overlap=nsnps)
            return(dat)
          } else{
            printer("COLOC:: No overlapping SNPs.")
            return(default_results)
          }  
    }) %>% data.table::rbindlist(fill=T) 
    return(coloc_dt)
  }) %>% data.table::rbindlist(fill=T)
 
  # Save
  if(save_results){
    coloc_path <- file.path("./Data/GWAS/",dataset,"_genome_wide","COLOC")
    dir.create(coloc_path, recursive = T, showWarnings = F)
    data.table::fwrite(COLOC_DT, file.path(coloc_path,"COLOC_results.txt"), sep="\t")
  } 
  COLOC.PP4_plot(COLOC_DT)
  return(COLOC_DT)
}

  


# MASHR <- function(){ 
#   devtools::install_github("stephenslab/mashr@v0.2-11")
#   library(ashr)
#   library(mashr)
#   set.seed(1)
#   simdata = simple_sims(500,5,1)
# }
