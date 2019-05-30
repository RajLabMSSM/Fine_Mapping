#     %%%%%%%%%%%%%%%%%      #
####### Fine-mapping ####### 
#     %%%%%%%%%%%%%%%%%     #


get_sample_size <- function(subset_DT, sample_size=NA){
  if(is.na(sample_size)){
    if("N_cases" %in% colnames(subset_DT) & "N_controls" %in% colnames(subset_DT)){
      sample_size <- max(subset_DT$N_cases) + max(subset_DT$N_controls)
    } else {
      sample_size <- 10000
    }
  }
  return(sample_size)
}



get_var_y <- function(subset_DT, dataset_type){ 
  if(dataset_type=="GWAS" & "N_cases" %in% colnames(subset_DT) & "N_controls" %in% colnames(subset_DT)){
    cat("\n Computing phenotype variance...") 
    phenotype_variance <- var(c(rep(0, max(subset_DT$N_cases)),
                                rep(1, max(subset_DT$N_controls)))
                              ) 
  } else if(dataset_type=="eQTL" & "Expression" %in% colnames(subset_DT)){ 
    phenotype_variance <- var(subset_DT$Expression)
  } else {
    cat("\n Estimating prior ") 
    phenotype_variance <- NA
  }
  return(list(phenotype_variance=phenotype_variance))
}



check_credible <- function(subset_DT, fitted_bhat){
  ## *** IMPORTANT! ***: In case susieR cannot identify any credible set, 
  # take the snps with the top 5 PIPs and provide a warning message. Interpret these snps with caution.
  
  # Credible_Set <- subset_DT[ as.numeric(strsplit( as.character(summary(fitted_bhat)$cs$variable) ,",")), ]$SNP  
  cat("\n ******",length(Credible_Set),"SNPs included in Credible Set ******\n") 
  return(Credible_Set)
}  
error_handling <- function(code) {
  tryCatch(code,
           error = function(c) {
             cat("\n--- ERROR ---")
             cat("\n ****** Could NOT identify credible set. Default to SNPs with the top 5 PIPs ******\n") 
             CS <- finemap_DT %>% arrange(desc(PIP))
             Credible_Set <- as.character(CS[1:5,]$SNP)
             return(Credible_Set)
           },
           warning = function(c) "Warning",
           message = function(c) "Message"
  )
}

SUSIE <- function(subset_DT, 
                  LD_matrix, 
                   dataset_type="GWAS",
                   n_causal=5,
                   sample_size=NA, 
                   var_y="estimate"){
  # Sum of Single Effects (SuSiE): Iterative Bayesian Step-wise Selection
  # https://stephenslab.github.io/susieR/
  vars <- get_var_y(subset_DT, dataset_type)
  sample_size <- get_sample_size(subset_DT, sample_size)
  
  cat("\n + Fine-mapping with SusieR... \n")
  fitted_bhat <- susieR::susie_bhat(bhat = subset_DT$Effect,
                            shat = subset_DT$StdErr,
                            R = LD_matrix,
                            n = sample_size, # Number of samples/individuals in the dataset
                            L = n_causal, # we assume there are at *most* 'L' causal variables
                            ## NOTE: setting L == 1 has a strong tendency to simply return the SNP with the largest effect size.
                            
                            estimate_prior_variance = TRUE, # default = FALSE
                            scaled_prior_variance = 0.1, # 0.1: Equates to "proportion of variance explained"
                            
                            standardize = TRUE,
                            estimate_residual_variance = TRUE, # TRUE
                            var_y = vars$phenotype_variance, # Variance of the phenotype (e.g. gene expression, or disease status)
                            
                            verbose = FALSE
  ) 
  cat("\n ++ Extracting Credible Sets...") 
  susie_snps <- names(fitted_bhat$X_column_scale_factors)
  CS_indices <- susieR::susie_get_cs(fitted_bhat)$cs
  Credible_Sets <- lapply(CS_indices, function(x){susie_snps[x]}) #%>% unlist()
  
  cat("\n ++ Merging susieR results with original data")
  # Create coordinates col
  # subset_DT$Coord <- paste(subset_DT$CHR, subset_DT$POS, sep=":")  
  res_DT <- data.table::data.table(SNP = names(fitted_bhat$X_column_scale_factors), 
                                   Probability = fitted_bhat$pip )
  finemap_DT <- data.table:::merge.data.table(x = data.table::data.table(subset_DT, key="SNP"), 
                                              y = data.table::data.table(res_DT, key="SNP")) 
  finemap_DT <- finemap_DT %>% arrange(desc(Probability))
  # Assign credible set #, or 0 to denote that it's not part of any credible set
  ## NOTE: if a SNP is part of more than one list, the top-ranked group to which is belong is used
  CS_dict <- list() 
  for(i in 1:length(Credible_Sets)){
    for(s in Credible_Sets[[i]]){
      CS_dict <- append(CS_dict, setNames(i,s))
    }  
  }
  finemap_DT$Credible_Set <- lapply(finemap_DT$SNP, function(x){ if(x %in% names(CS_dict)){ CS_dict[[x]] } else{0}}) %>% unlist()   
  return(finemap_DT)
}


# ABF 
ABF <- function(subset_DT, PP_threshold=.5){
  cat("\n Fine-mapping with ABF... \n") 
  #data.table::fread("Data/GWAS/Nalls23andMe_2019/LRRK2/LRRK2_Nalls23andMe_2019_subset.txt") 
  finemap_DT <- coloc::finemap.abf(dataset = list(beta = subset_DT$Effect,
                                           varbeta = subset_DT$StdErr^2, # MUST be squared
                                           N = length(subset_DT$Effect),
                                           s = subset_DT$proportion_cases,
                                           snp = subset_DT$SNP,
                                           MAF = subset_DT$MAF, 
                                           type="cc")) 
  finemap_DT <- subset(finemap_DT, snp!="null") %>% rename(SNP=snp, Probability=SNP.PP) %>%
    arrange(desc(Probability))
  # Arbitarily assign SNPs with the top N probability as the top candidate causal SNPs
  # finemap_DT$Credible_Set <- c(rep(1,n_causal), rep(0,dim(finemap_DT)[1]-n_causal))
  # Any SNPs with a PP greater than the set threshold get included in the credible set
  finemap_DT$Credible_Set <- ifelse(finemap_DT$Probability > PP_threshold, 1, 0)
  finemap_DT <- data.table:::merge.data.table(x=subset_DT, 
                                              y=subset(finemap_DT, select=c("SNP","Probability","Credible_Set")), 
                                              on="SNP")
  finemap_DT <- finemap_DT %>% arrange(desc(Probability))
  return(finemap_DT)
} 





#### FINEMAP ####   
construct_FINEMAP_data <- function(results_path, subset_DT){
  ####### data.z ####### 
  cat("\n++ Formatting data.z file for FINEMAP")
  data.z <- subset_DT %>% dplyr::select(rsid=SNP, 
                                        chromosome=CHR, 
                                        position=POS,
                                        allele1=A1, 
                                        allele2=A2, 
                                        maf=MAF, 
                                        beta=Effect, # *required
                                        se=StdErr # *required
                                        )
  data.z$flip <- 0 # [optional] - flip==1, don't flip==0
  
  ####### data.ld #######
  cat("\n++ Formatting LD Matrix for FINEMAP")
  ## The order of the SNPs in the dataset.ld must correspond to the order of variants in dataset.z.
  load(file.path(results_path,"plink","LD_matrix.RData")) 
  
  # Filter 
  data.z <- subset(data.z, rsid %in% rownames(LD_matrix))
  ## This filters AND sorts LD_matrix by the order of rsids in data.z
  LD_filt <- LD_matrix[rownames(LD_matrix) %in% data.z$rsid, colnames(LD_matrix) %in% data.z$rsid] 
  
  # Write files
  ## MUST be space-delimited
  cat("\n+++ Writing FINEMAP z and ld files...")
  if( dim(data.z)[1]==dim(LD_filt)[1] ){
    # data.z
    data.z_path <- file.path(results_path,"FINEMAP","data.z")
    data.table::fwrite(data.z, data.z_path, sep = " ") 
    # Sys.chmod(data.z_path, "777", use_umask = FALSE)
    # data.ld
    data.ld_path <- file.path(results_path,"FINEMAP","data.ld")
    data.table::fwrite(data.table:::as.data.table.matrix(LD_filt),
                       data.ld_path, sep=" ", quote = F, col.names = F)
    # Sys.chmod(data.ld_path, "777", use_umask = FALSE)
  } else {warning("\n --- FINEMAP: Summary statistics file (data.z) and LD matrix (data.ld) must contain the same number of SNPs.---")}
}

construct_FINEMAP_master <- function(results_path,   
                                     n_samples,
                                     dataset_number=1,
                                     file.k=NA){ # [optional input]){
  cat("\n ++ Constructing FINEMAP master file.")
  # For full list of parameters: http://www.christianbenner.com 
  header <- "z;ld;snp;config;cred;log;n_samples"
  # pathList <-  paste(c(file.z, file.ld, file.snp, file.config, file.log, n_samples), collapse=";")
  files <- c("data.z",  # [required input]
             "data.ld", # [required input]
             "data.snp", # [output]
             "data.config", # [optional output]
             "data.cred", # [optional output]
             "data.log"# [optional output]
             )
  if(!is.na(file.k)){ pathList <- append(pathList, file.k) } 
  paths_list <- paste(c(file.path("FINEMAP",files),n_samples), collapse = ";")  
  # Write master file
  dir.create(file.path(results_path, "FINEMAP"), recursive = T, showWarnings = F)
  data.table::fwrite(list(header,paths_list), file.path(results_path,"FINEMAP","master"), quote=F, sep="\n")
}


process_FINEMAP_results <- function(results_path, subset_DT){
  # Import credible sets
  top_config <- data.table::fread(file.path(results_path,"FINEMAP/data.config"), sep=" ", nrows = 1) 
  Credible_Set <- strsplit(top_config$config, ",")[[1]]
  # Import snp-level results
  snp_level <- data.table::fread(file.path(results_path,"FINEMAP/data.snp"), sep=" ")
  # Merge with original data 
  subset_DT$Credible_Set <- ifelse(subset_DT$SNP %in% Credible_Set, 1, 0)
  subset_DT <- data.table:::merge.data.table(subset_DT, 
                                             subset(snp_level, select=c("rsid","prob")), by.x = "SNP", by.y="rsid")
  subset_DT <- subset_DT %>% rename(Probability=prob) %>% arrange(desc(Credible_Set))
  return(subset_DT)
}



FINEMAP <- function(subset_DT,
                    results_path,
                    FINEMAP_path="./echolocatoR/tools/FINEMAP/finemap_v1.3_MacOSX",
                    n_samples=NA,
                    n_causal=5,# Max number of allowed causal SNPs
                    model="cond" # cond (stepwise conditional search) vs. sss (stochastic shotgun search)
                    ){ 
  # The stepwise conditional search starts with a causal configuration containing the 
  ## SNP with the lowest P-value alone and then iteratively adds to the causal configuration 
  ## the SNP given the highest posterior model probability until no further SNP yields
  ## a higher posterior model probability.
  if(is.na(n_samples) & "N_cases" %in% colnames(subset_DT) & "N_controls" %in% colnames(subset_DT)){
    n_samples <- max(subset_DT$N_cases) + max(subset_DT$N_controls)
    }
  # Setup files
  construct_FINEMAP_master(results_path = results_path, n_samples = n_samples)
  construct_FINEMAP_data(results_path = results_path, subset_DT = subset_DT) 
  # Command line
  ## Example: 
  ## cmd <- paste(FINEMAP_path," --sss --in-files",file.path(dirname(FINEMAP_path),"example","master"), "--dataset 1 --n-causal-snps 5") 
  file.copy(from=FINEMAP_path, to=file.path(results_path))
  cmd <- paste("cd",results_path,"&&",
               "./finemap_v1.3_MacOSX",
               paste0("--",model),
               "--in-files",file.path("FINEMAP/master"),
               "--log",
               "--n-causal-snps",n_causal)
  cat(cmd)
  system(cmd) 
  file.remove(file.path(results_path,"finemap_v1.3_MacOSX"))
  # Process results
  finemap_DT <- process_FINEMAP_results(results_path, subset_DT)
  return(finemap_DT)
}
 


find_consensus_SNPs <- function(finemap_DT){
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

# Fine-map using multiple fine-mapping tools
multi_finemap <- function(results_path, 
                          finemap_method_list, 
                          subset_DT, 
                          dataset_type,
                          LD_matrix=NULL, 
                          n_causal=5, 
                          sample_size){
  cat("\n Fine-mapping using multiple tools:\n  ", paste(finemap_method_list, collapse=", "),"\n")
  # finemap_method_list <- finemap_method_list[finemap_method_list!="COJO"]
  
  merged_DT <- subset_DT 
  for(i in 1:length(finemap_method_list)){  
    m <- finemap_method_list[i] 
    DT <- finemap_methods(results_path = results_path,
                          finemap_method = m, 
                          subset_DT = subset_DT, 
                          dataset_type = dataset_type,
                          LD_matrix = LD_matrix, 
                          n_causal = n_causal, 
                          sample_size = sample_size) 
    # Add results to method-specific columns
    cat("\n --- Merging",m,"results with multi-finemap data.")
    value_var <- if(m=="COJO"){"Conditioned_Effect"}else{"Probability"}
    DT_sub <- subset(DT, select = c("SNP","Credible_Set",value_var) )
    cols <- colnames(DT_sub)
    colnames(DT_sub) <- c("SNP", paste(m, cols[2:length(cols)], sep="." ))
    merged_DT <- data.table:::merge.data.table(merged_DT, DT_sub, by="SNP") 
    # Arrange by newest probability measure
    # merged_DT <- merged_DT %>% arrange(desc(eval(parse(text=paste0(m,".Probability") ))))
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
  data.table::fwrite(finemap_DT, file_dir, sep = "\t")
}


finemap_methods <- function(results_path,
                            finemap_method="SUSIE", 
                            subset_DT, 
                            dataset_type="GWAS",
                            force_new_finemap=T,
                            LD_matrix=NULL, 
                            n_causal=5, 
                            sample_size,
                            conditioned_snps){
  cat("\n",finemap_method) 
  finemap_DT <- data.table::data.table() #Initialize DT
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
    finemap_DT <- COJO(results_path = results_path, 
                       subset_DT = subset_DT, 
                       conditioned_snps = conditioned_snps, 
                       conditional_analysis = T, 
                       stepwise_procedure = T) 
    
  } else {
    cat("\n [::ERROR::] Enter valid finemap_method: 'SUSIE', 'ABF', 'FINEMAP', and 'COJO' are currently available. \n")
  }
}


finemap_handler <- function(results_path,
                            finemap_method="SUSIE", 
                            subset_DT, 
                            dataset_type="GWAS",
                            force_new_finemap=T,
                            LD_matrix=NULL, 
                            n_causal=5, 
                            sample_size,
                            conditioned_snps){
  cat("\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n")
  cat("  -------- Step 4: Statistically Fine-map --------  ")
  cat("\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n")
  start_FM <- Sys.time()
  set.seed(1)  
  # First, check if there's more than one fin-mapping method given. If so, switch to multi-finemap function
  if(length(finemap_method)>1){ 
    ## Next, see if fine-mapping has previously been done (with multi-finemap)
    file_dir <- create_method_dir(results_path, "Multi-finemap")
    ### If so, import the previous results
    if(file.exists(file_dir) & force_new_finemap==F){
      cat("\n +++ Previously multi-fine-mapped results identified. Importing...")
      finemap_DT <- data.table::fread(file_dir, sep="\t") %>% data.table::data.table() 
    } else {
      ### If not, or if forcing new fine-mapping is set to TRUE, fine-map using multiple tools
      finemap_DT <- multi_finemap(results_path = results_path, 
                                 finemap_method_list = finemap_method,
                                 subset_DT = subset_DT, 
                                 dataset_type = dataset_type,
                                 LD_matrix = LD_matrix,
                                 n_causal = n_causal, 
                                 sample_size = sample_size)
      save_finemap_results(finemap_DT, file_dir)
    } 
  } else {
    # If only one tool is given, fine-map using only that tool
    ## Next, see if fine-mapping has previously been done (with this tool).
    file_dir <- create_method_dir(results_path, finemap_method)
    ### If so, import the previous results
    if(file.exists(file_dir) & force_new_finemap==F){
      cat("\n +++ Previously fine-mapped", finemap_method, "results identified. Importing...") 
      finemap_DT <- data.table::fread(file_dir, sep="\t") %>% data.table::data.table() 
    } else{ 
      ### If not, or if forcing new fine-mapping is set to TRUE, fine-map with that given tool
      finemap_DT <- finemap_method(results_path = results_path, 
                                   finemap_method_list = finemap_method,
                                   subset_DT = subset_DT, 
                                   dataset_type = dataset_type,
                                   LD_matrix = LD_matrix,
                                   n_causal = n_causal, 
                                   sample_size = sample_size) 
      save_finemap_results(finemap_DT, file_dir) 
    }
  } 
  end_FM <- Sys.time()
  cat("\n + Fine-mapping with '", paste0(finemap_method, collapse=", "),"' completed in ",round(end_FM-start_FM,2)," seconds.", sep="")
  return(finemap_DT)
}




multi_finemap_results_table <- function(results_path, finemap_method_list, fancy_table=T){
  finemap_DT <- data.table::fread(file.path(results_path,"Multi-finemap/Multi-finemap_results.txt"),stringsAsFactors = F)
  CS_cols <- paste0(finemap_method_list,".Credible_Set")
  support_DT <- subset(finemap_DT, Support>0, select=c("SNP","CHR","POS","P",CS_cols,"Support")) %>% 
    arrange(desc(Support))
 
  
  # Plot table 
  if(fancy_table){
    customGreen0 = "#DeF7E9" 
    customGreen = "#71CA97" 
    customRed = "#ff7f7f" 
  CS_formatter <- 
    formattable::formatter("span", 
                style = x ~ style( 
                  color = ifelse(x > 0, customGreen, ifelse(x == 0, "black", "black")))) 
  formattable::formattable(support_DT, 
                align =c("l","c","c","c","c", "c", "c", "c", "r"), 
                list( P = formattable::color_tile(customGreen, customGreen0),
                      SUSIE.Credible_Set = CS_formatter,
                      ABF.Credible_Set = CS_formatter,
                      FINEMAP.Credible_Set = CS_formatter,
                      COJO.Credible_Set = CS_formatter,
                      Support = formattable::color_tile("white", "green")) 
                )
   
   } else {
    print(support_DT)
  }
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
#   cat("export PATH=~/GIT/github/mvarbvs/dsc/modules/linux:$PATH",
#       "Rscript ~/GIT/susieR/inst/code/finemap.R input=\"N2.with_sumstats.rds\" output=\"N2finemapping.FINEMAP\" args=\"--n-causal-max\ 2\" &> /dev/null")
# 
#   system("export PATH=~/GIT/github/mvarbvs/dsc/modules/linux:$PATH
#   Rscript ~/GIT/susieR/inst/code/caviar.R input=\"N2.with_sumstats.rds\" output=\"N2finemapping.CAVIAR\" args=\"-c\ 2\ -g\ 0.01\" &> /dev/null")
# 
#   system("export PATH=~/GIT/github/mvarbvs/dsc/modules/linux:$PATH
#   python ~/GIT/susieR/inst/code/dap-g.py N2.with_sumstats.rds N2finemapping.DAP -ld_control 0.20 --all &> /dev/null")
# 
# }


