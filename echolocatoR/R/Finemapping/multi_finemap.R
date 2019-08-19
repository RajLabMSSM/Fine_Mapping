#///////////////#
# Multi-finemap # 
#///////////////#

find_consensus_SNPs <- function(finemap_DT, verbose=T){
  printer("+++ Identifying Consensus SNPs...",v=verbose)
  # Find SNPs that are in the credible set for all fine-mapping tools
  CS_cols <- colnames(finemap_DT)[endsWith(colnames(finemap_DT),".Credible_Set")] 
  consensus_SNPs <- finemap_DT %>%  
    as.data.frame() %>%
    dplyr::filter_at(.vars = CS_cols, all_vars(.>0))
  finemap_DT$Consensus_SNP <- finemap_DT$SNP %in% consensus_SNPs$SNP 
  finemap_DT <- finemap_DT %>% arrange(desc(Consensus_SNP))
  # Number of tools supporting
  CS_DT <- subset(finemap_DT, select = CS_cols) 
  finemap_DT$Support <- rowSums(CS_DT > 0)
  return(finemap_DT)
}

# results_path <- "Data/GWAS/Kunkle_2019/PTK2B"
# subset_DT <- fread( file.path(results_path,"PTK2B_Kunkle_2019_subset.txt"))
# dataset_type <- "GWAS"
# load(file.path(results_path, "plink/LD_matrix.RData"))
# 


# Fine-map using multiple fine-mapping tools
multi_finemap <- function(results_path,
                          fullSS_path,
                          finemap_method_list, 
                          subset_DT, 
                          dataset_type,
                          LD_matrix=NULL, 
                          n_causal=5, 
                          sample_size,
                          snp_col="SNP",
                          freq_col="Freq",
                          effect_col="Effect",
                          stderr_col="StdErr",
                          pval_col="P",
                          N_cases_col="N_cases",
                          N_controls_col="N_controls",
                          A1_col="A1",
                          A2_col="A2"){
  printer("++ Fine-mapping using multiple tools:", paste(finemap_method_list, collapse=", "))
  # finemap_method_list <- finemap_method_list[finemap_method_list!="COJO"] 
 
  select_cols <- colnames(subset_DT)[!grepl(colnames(subset_DT), 
                                            pattern = paste(c(finemap_method_list,"Support","Consensus_SNP"),collapse = "|"))]
  merged_DT <- subset(subset_DT, select = select_cols)
  
  
  for(i in 1:length(finemap_method_list)){  
    m <- finemap_method_list[i];
    message("Multi-finemap:: ",m)
    DT <- null_DT <- data.table::data.table(SNP = merged_DT$SNP, Credible_Set = NA, Probability = NA);
    # DT <- tryCatch({
      # EXPRESSION
     try({
      DT <- finemap_method_handler(fullSS_path = fullSS_path,
                                   results_path = results_path,
                                   finemap_method = m, 
                                   subset_DT = subset_DT, 
                                   dataset_type = dataset_type,
                                   LD_matrix = LD_matrix, 
                                   n_causal = n_causal, 
                                   sample_size = sample_size,
                                   snp_col = snp_col,
                                   freq_col = freq_col,
                                   effect_col = effect_col,
                                   stderr_col = stderr_col,
                                   pval_col = pval_col,
                                   N_cases_col = N_cases_col,
                                   N_controls_col = N_controls_col,
                                   A1_col = A1_col,
                                   A2_col = A2_col)
     }) 
      # },
      # # WARNING
      # warning = function(w){printer("WARNING")},
      # # ERROR
      # error = function(e){
      #   message("SUSIE::Error:: Could not identify Credible Set."); 
      #   return(null_DT)
      #   } 
      # ) ## End tryCatch
      try({printer("++ SUSIE: Credible Set SNPS identified =",sum(!is.na(DT$Credible_Set)) )})
      # Add results to method-specific columns
      printer("+++ Merging",m,"results with multi-finemap data.");
      value_var <- if(m=="COJO"){"Conditioned_Effect"}else{"Probability"};
      DT_select <- subset(DT, select = c("SNP","Credible_Set",value_var) );
      # Rename columns according to method name
      cols <- colnames(DT_select);
      colnames(DT_select) <- c("SNP", paste(m, cols[2:length(cols)], sep="." ));
      
      # Merge new columns into DT
      merged_DT <- data.table:::merge.data.table(data.table::as.data.table(merged_DT),
                                                 data.table::as.data.table(DT_select),
                                                 by="SNP", all = T);

  }  
  finemap_DT <- find_consensus_SNPs(merged_DT)  
  return(finemap_DT)
}
 

####### Fine-mapping Handler ####### 
create_method_dir <- function(results_path, finemap_method){
  method_dir <- file.path(results_path, finemap_method)
  # Make finemapping results folder 
  dir.create(method_dir, recursive = T, showWarnings = F)
  # Return results file name
  file_dir <- file.path(method_dir,paste0(finemap_method,"_results.txt"))
  return(file_dir)
}

save_finemap_results <- function(finemap_DT, file_dir){ 
  data.table::fwrite(finemap_DT, file_dir, sep = "\t", na = NA, quote = F)
}


finemap_method_handler <- function(results_path,
                                   fullSS_path,
                                    finemap_method="SUSIE", 
                                    subset_DT, 
                                    dataset_type="GWAS",
                                    force_new_finemap=T,
                                    LD_matrix=NULL, 
                                    n_causal=5, 
                                    sample_size,
                                    conditioned_snps,
                                    snp_col="SNP",
                                    freq_col="Freq",
                                    effect_col="Effect",
                                    stderr_col="StdErr",
                                    pval_col="P",
                                    N_cases_col="N_cases",
                                    N_controls_col="N_controls",
                                    A1_col="A1",
                                    A2_col="A2"){
  printer("\n",finemap_method)  
  # INITIATE FINE-MAPPING 
  if(finemap_method=="SUSIE"){ 
    # SUSIE
    finemap_DT <- SUSIE(subset_DT = subset_DT, 
                        dataset_type = dataset_type,
                        LD_matrix = LD_matrix,
                        n_causal = n_causal, 
                        sample_size = sample_size)
    
    
  } else if(finemap_method=="ABF"){
    # coloc - finemap.abf
    finemap_DT <- ABF(subset_DT = subset_DT,
                      PP_threshold = .5)
    
    
  } else if(finemap_method=="FINEMAP"){
    # FINEMAP
    finemap_DT <- FINEMAP(subset_DT = subset_DT,
                          results_path = results_path, 
                          n_samples = sample_size,   
                          n_causal = n_causal)
    
    
  } else if("COJO" %in% finemap_method){
    #COJO
    conditioned_snps <- subset(subset_DT, leadSNP==T)$SNP
    finemap_DT <- COJO(subset_DT = subset_DT,
                       results_path = results_path, 
                       fullSS_path = fullSS_path, 
                       conditioned_snps = conditioned_snps, 
                       conditional_analysis = T, 
                       stepwise_procedure = F,
                       
                       snp_col = snp_col,
                       freq_col = freq_col,
                       effect_col = effect_col,
                       stderr_col = stderr_col,
                       pval_col = pval_col,
                       N_cases_col = N_cases_col,
                       N_controls_col = N_controls_col,
                       A1_col = A1_col,
                       A2_col = A2_col) 
    
  } else {
    stop("[::ERROR::] Enter valid finemap_method: 'SUSIE', 'ABF', 'FINEMAP', 'COJO', and 'PAINTOR' are currently available.")
  }
  return(finemap_DT)
}


finemap_handler <- function(results_path,
                            fullSS_path,
                            finemap_methods=c("SUSIE","FINEMAP"), 
                            subset_DT, 
                            dataset_type="GWAS",
                            force_new_finemap=T,
                            LD_matrix=NULL, 
                            n_causal=5, 
                            sample_size,
                            conditioned_snps,
                            snp_col="SNP",
                            freq_col="Freq",
                            effect_col="Effect",
                            stderr_col="StdErr",
                            pval_col="P",
                            N_cases_col="N_cases",
                            N_controls_col="N_controls",
                            A1_col="A1",
                            A2_col="A2"){
  message("-------- Step 4: Statistically Fine-map --------")
  start_FM <- Sys.time()
  set.seed(1)  
  # First, check if there's more than one fin-mapping method given. If so, switch to multi-finemap function
  # if(length(finemap_methods)>1){ 
    ## Next, see if fine-mapping has previously been done (with multi-finemap)
    file_dir <- create_method_dir(results_path, "Multi-finemap")
    ### If so, import the previous results
    if(file.exists(file_dir) & force_new_finemap==F){
      printer("++ Previously multi-fine-mapped results identified. Importing...")
      finemap_DT <- data.table::fread(file_dir, sep="\t") %>% data.table::data.table() 
    } else {
      ### If not, or if forcing new fine-mapping is set to TRUE, fine-map using multiple tools
      finemap_DT <- multi_finemap(results_path = results_path, 
                                 fullSS_path = fullSS_path,
                                 finemap_method_list = finemap_methods,
                                 subset_DT = subset_DT, 
                                 dataset_type = dataset_type,
                                 LD_matrix = LD_matrix,
                                 n_causal = n_causal, 
                                 sample_size = sample_size,
                                 
                                 snp_col = snp_col,
                                 freq_col = freq_col,
                                 effect_col = effect_col,
                                 stderr_col = stderr_col,
                                 pval_col = pval_col,
                                 N_cases_col = N_cases_col,
                                 N_controls_col = N_controls_col,
                                 A1_col = A1_col,
                                 A2_col = A2_col)
      save_finemap_results(finemap_DT, file_dir)
    } 
  # } else {
  #   # If only one tool is given, fine-map using only that tool
  #   ## Next, see if fine-mapping has previously been done (with this tool).
  #   file_dir <- create_method_dir(results_path, finemap_methods)
  #   ### If so, import the previous results
  #   if(file.exists(file_dir) & force_new_finemap==F){
  #     printer("\n +++ Previously fine-mapped", finemap_methods, "results identified. Importing...") 
  #     finemap_DT <- data.table::fread(file_dir, sep="\t") %>% data.table::data.table() 
  #   } else{ 
  #     ### If not, or if forcing new fine-mapping is set to TRUE, fine-map with that given tool
  #     finemap_DT <- finemap_method_handler(results_path = results_path, 
  #                                  finemap_methods = finemap_methods,
  #                                  subset_DT = subset_DT, 
  #                                  dataset_type = dataset_type,
  #                                  LD_matrix = LD_matrix,
  #                                  n_causal = n_causal, 
  #                                  sample_size = sample_size) 
  #     save_finemap_results(finemap_DT, file_dir) 
  #   }
  # } 
  end_FM <- Sys.time()
  printer("++ Fine-mapping with '", paste0(finemap_methods, collapse=", "),"' completed in ",round(end_FM-start_FM,2)," seconds.", sep="")
  return(finemap_DT)
}

 
# 
# run_finemapping_tools <- function(X, Y, Z=NA){
#   # https://stephenslab.github.io/susie-paper/manuscript_results/pedagogical_example.html
#   data("N2finemapping")
#   attach(N2finemapping)
#   dat = N2finemapping
#   names(dat)
#   
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
#   printer("export PATH=~/GIT/github/mvarbvs/dsc/modules/linux:$PATH",
#       "Rscript ~/GIT/susieR/inst/code/finemap.R input=\"N2.with_sumstats.rds\" output=\"N2finemapping.FINEMAP\" args=\"--n-causal-max\ 2\" &> /dev/null")
# 
#   system("export PATH=~/GIT/github/mvarbvs/dsc/modules/linux:$PATH
#   Rscript ~/GIT/susieR/inst/code/caviar.R input=\"N2.with_sumstats.rds\" output=\"N2finemapping.CAVIAR\" args=\"-c\ 2\ -g\ 0.01\" &> /dev/null")
# 
#   system("export PATH=~/GIT/github/mvarbvs/dsc/modules/linux:$PATH
#   python ~/GIT/susieR/inst/code/dap-g.py N2.with_sumstats.rds N2finemapping.DAP -ld_control 0.20 --all &> /dev/null")
# 
# }

 


