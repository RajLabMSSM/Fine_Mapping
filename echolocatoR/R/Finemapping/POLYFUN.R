# $$$$$$$$$$$$$$$$ $$$$$$$ $$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$ PolyFun $$$$$$$$$$$$$$$$ 
# $$$$$$$$$$$$$$$$ $$$$$$$ $$$$$$$$$$$$$$$$
# https://github.com/omerwe/polyfun

#####-------------------------------------------------------------
# There are three ways to run PolyFun:
  
# 1. Using precomputed prior causal probabilities of 19 million imputed UK Biobank SNPs with MAF>0.1%, based on a meta-analysis of 15 UK Biobank traits. This is the simplest approach, but it may not include all your SNPs of interest (especially when analyzing non-European populations) and the prior causal probabilities may not be optimal for some traits.

# 2. Computing prior causal probabilities via an L2-regularized extension of stratified LD-score regression (S-LDSC). This is a relatively simple approach, but the prior causal probabilities may not be robust to modeling misspecification.

# 3. Computing prior causal probabilities non-parametrically. This is the most robust approach, but it is computationally intensive and requires access to individual-level genotypic data from a large reference panel (optimally >10,000 population-matched individuals). 
#####-------------------------------------------------------------



# NOTES for Omer:
     
# 6. polyfun.py script reduced the number of SNPs from 7764212 to 34851, printing "[WARNING]  number of SNPs is smaller than 200k; this is almost always bad."  Is all the filtering just because the example reference set is small?
# 9. If a variant has been filtered out by PolyFun at some point, is it appropriate to assign that variant a prior prob value of 0? I'm looking for ways to retain as many SNPs as possible for the fine-mapping analysis so I can compare functional fine-mapping (e.g. PolyFun+SUSIE) vs. statistical fine-mapping (e.g. SUSIE).
 
# GitHub Notes:
## How to pull changes from original repo into the forked repo
# https://digitaldrummerj.me/git-syncing-fork-with-original-repo/


POLYFUN.help <- function(){
  cmd <- paste("python3",file.path(polyfun,"polyfun.py"),
               "--help")
  system(cmd)
}


# %%%%%%%%%%%%%%%% PolyFun approach 1 %%%%%%%%%%%%%%%% 
## Using precomputed prior causal probabilities based on a meta-analysis of 15 UK Biobank traits

POLYFUN.prepare_snp_input <- function(PF.output.path,
                                       results_path="./Data/GWAS/Nalls23andMe_2019/LRRK2",
                                       finemap_DT=NULL){
  # PolyFun requires a space-delimited (gzipped or not) file with these columns:
  ## CHR BP A1 A2
  printer("PolyFun:: Preparing SNP input file...")  
  PF.dat <- dplyr::select(finemap_DT, SNP, CHR, BP=POS, A1, A2)
  printer("+ PolyFun::",nrow(PF.dat),"SNPs identified.")
  snp.path <- file.path(PF.output.path,"snps_to_finemap.txt.gz")
  printer("+ PolyFun:: Writing SNP file ==>",snp.path)
  data.table::fwrite(PF.dat, file = snp.path, 
                     nThread = 4, sep = " ")
  return(snp.path)
}


POLYFUN.initialize <- function(results_path, 
                               finemap_DT=NULL, 
                               dataset="Nalls23andMe_2019", 
                               locus="_genome_wide"){
  # Create path
  PF.output.path <- file.path(results_path, "PolyFun")
  dir.create(PF.output.path, showWarnings = F, recursive = T)
  # Import SNPs
  if(is.null(finemap_DT)){
    if(is.null(locus)){
      printer("POLYFUN:: Importing summary stats: genome-wide")
      finemap_DT <- data.table::fread(Directory_info(dataset, "fullSS.local"), nThread = 4)
    }
    printer("POLYFUN:: Importing summary stats:",locus)
    finemap_DT <- data.table::fread(file.path(dirname(Directory_info(dataset, "fullSS.local")),locus,
                                            "Multi-finemap/Multi-finemap_results.txt"), nThread = 4) 
  } 
  return(finemap_DT)
}

POLYFUN.get_precomputed_priors <- function(polyfun="./echolocatoR/tools/polyfun", 
                                           results_path="./Data/GWAS/Nalls23andMe_2019/LRRK2",
                                           dataset="Nalls23andMe_2019",
                                           finemap_DT=NULL, 
                                           locus="LRRK2"){
  PF.output.path <- file.path(results_path, "PolyFun")
  finemap_DT <- POLYFUN.initialize(finemap_DT=finemap_DT, 
                                    results_path=results_path, 
                                    locus=locus)
  # Prepare input
  snp.path <- POLYFUN.prepare_snp_input(PF.output.path=PF.output.path, 
                                        results_path=results_path,
                                        finemap_DT=finemap_DT)
  PF.output.file <- file.path(PF.output.path,"snps_with_priors.snpvar.gz")
  # Retrieve priors
  cmd <- paste("python3",file.path(polyfun,"extract_snpvar.py"),
        "--snps",snp.path,
        "--out",file.path(PF.output.path,"snps_with_priors"))
  print(cmd)
  system(cmd) 
  # Import results
  priors <- data.table::fread(PF.output.file, nThread = 4) 
  return(priors)
}


 

# %%%%%%%%%%%%%%%% PolyFun approaches 2 & 3 %%%%%%%%%%%%%%%% 
## 2. Computing prior causal probabilities via an L2-regularized extension of S-LDSC
## 3. Computing prior causal probabilities non-parametrically
 


POLYFUN.munge_summ_stats <- function(polyfun="./echolocatoR/tools/polyfun",  
                                     dataset="Nalls23andMe_2019",
                                     sample.size=1474097,
                                     min_INFO=0,
                                     min_MAF=0){
  results_path <- file.path(dirname(Directory_info(dataset_name = dataset, "fullSS.local")), "_genome_wide")
  PF.output.path <- file.path(results_path, "PolyFun")
  munged.path <- file.path(PF.output.path,"sumstats_munged.parquet") 
  # Requires space-delimited file with the following columns (munging can recognize several variations of these names):
  ## SNP CHR BP ....and....
  ## either a p-value, an effect size estimate and its standard error, a Z-score or a p-value 
  if(!file.exists(munged.path)){ 
    printer("+ PolyFun:: Initiating data munging pipeline...")
    cmd <- paste("python3",file.path(polyfun,"munge_polyfun_sumstats.py"),
                 "--sumstats",Directory_info(dataset_name = dataset, "fullSS.local"),
                 "--n",sample.size, # Study sample size
                 "--out",munged.path,
                 "--min-info",min_INFO,
                 "--min-maf", min_MAF)
    print(cmd)
    system(cmd)
  } else {printer("+ PolyFun:: Existing munged summary stats files detected.")} 
}



POLYFUN.gather_ldscores <- function(output_prefix){
  ldscore.files <-  list.files(dirname(output_prefix), pattern = ".l2.ldscore.parquet", full.names = T)
  parquor <- POLYFUN.read_parquet(ldscore.files[1])
}



POLYFUN.read_parquet <- function(parquet_path){
  # Converts parquet to data.table
  SparkR::sparkR.session()
  parquor <- SparkR::read.parquet(parquet_path)
  parquor <- SparkR::as.data.frame(parquor) %>% 
    data.table::data.table() 
  return(parquor)
}


POLYFUN.compute_priors <- function(polyfun="./echolocatoR/tools/polyfun",
                                    PF.output.path,
                                    munged.path, 
                                    min_INFO = 0.6,
                                    min_MAF = 0.05,
                                    annotations.path=file.path(polyfun,"example_data/annotations."),
                                    weights.path=file.path(polyfun,"example_data/weights."), 
                                    prefix="PD_GWAS",
                                    chrom="all",
                                    compute_ldscores=F, 
                                    allow_missing_SNPs=T){
  # Quickstart:
  # polyfun="./echolocatoR/tools/polyfun"; parametric=T;  weights.path=file.path(polyfun,"example_data/weights."); annotations.path=file.path(polyfun,"example_data/annotations."); munged.path= "./Data/GWAS/Nalls23andMe_2019/_genome_wide/PolyFun/sumstats_munged.parquet"; parametric=T; dataset="Nalls23andMe_2019"; prefix="PD_GWAS"; compute_ldscores=F; allow_missing_SNPs=T; chrom="all"
  
  # 0. Create paths
  results_path <- file.path(dirname(Directory_info(dataset_name = dataset, "fullSS.local")), "_genome_wide")
  PF.output.path <- file.path(results_path, "PolyFun")
  dir.create(PF.output.path, showWarnings = F, recursive = T) 
  out.path <- file.path(PF.output.path,"output")
  output_prefix <- file.path(out.path, prefix)
  dir.create(out.path, showWarnings = F, recursive = T)
  
  # 1. Munge summary stats
  printer("PolyFun:: [1]  Create a munged summary statistics file in a PolyFun-friendly parquet format.")
  POLYFUN.munge_summ_stats(polyfun=polyfun,
                         dataset="Nalls23andMe_2019",
                         sample.size=1474097, 
                         min_INFO = 0.6,
                         min_MAF = 0.05)  
  
  # 2. 
  ## If compute_ldscores == F:
  # This will create 2 output files for each chromosome: output/testrun.<CHR>.snpvar_ridge.gz and output/testrun.<CHR>.snpvar_ridge_constrained.gz. The first contains estimated per-SNP heritabilities for all SNPs (which can be used for downstream analysis with PolyFun; see below), and the second contains truncated per-SNP heritabilities, which can be used directly as prior causal probabilities in fine-mapping.  
  printer("PolyFun:: [2] Run PolyFun with L2-regularized S-LDSC")
  cmd2 <- paste("python3",file.path(polyfun,"polyfun.py"),
                "--compute-h2-L2",
               # Approach 2 = Parametric = no partitions = T
               # Approach 3 = Non-parametric = partitions = F 
                ifelse(compute_ldscores,"","--no-partitions"),
                "--output-prefix",output_prefix,
                "--sumstats",munged.path,
                "--ref-ld-chr",annotations.path,
                "--w-ld-chr",weights.path,
                ifelse(allow_missing_SNPs,"--allow-missing",""))
  print(cmd2)
  system(cmd2)
  
  # Computationally intensive: can parallelize by chromosomes
  if(compute_ldscores){
    # 3.
    printer("PolyFun:: [3] Compute LD-scores for each SNP bin")
    cmd3 <- paste("python3",file.path(polyfun,"polyfun.py"),
                  "--compute-ldscores",
                  "--output-prefix",output_prefix,
                  "--bfile-chr",file.path(dirname(annotations.path),"reference."),
                  ifelse(chrom=="all","",paste("--chr",chrom)),
                  ifelse(allow_missing_SNPs,"--allow-missing",""))
    print(cmd3)
    system(cmd3)
    
    # 4.
    printer("PolyFun:: [4] Re-estimate per-SNP heritabilities via S-LDSC")
    cmd4 <- paste("python3",file.path(polyfun,"polyfun.py"),
                  "--compute-h2-bins",
                  "--output-prefix",output_prefix,
                  "--sumstats",munged.path,
                  "--w-ld-chr",weights.path, 
                  ifelse(allow_missing_SNPs,"--allow-missing",""))
    print(cmd4)
    system(cmd4)
    
    printer("PolyFun:: Results directory =",dirname(output_prefix))
    printer("PolyFun:: Results files:")
    printer("          *.snpvar_ridge.gz")
    printer("          *.snpvar_ridge_constrained.gz") 
    # The output of the partitioned LDSC has the suffix: .snpvar_constrained.gz (one per chrom)
    LDSC.files <- list.files(out.path, pattern = "*.snpvar_constrained.gz", full.names = T)
    # pd_ldsc <- data.table::fread(PS_LDSC.files[1], nThread = 4) 
    # ldscore <- POLYFUN.read_parquet(file.path(out.path,"PD_GWAS.1.l2.ldscore.parquet")) 
    # bin.1 <- POLYFUN.read_parquet(file.path(out.path,"PD_GWAS.2.bins.parquet"))
    #rowSums(bin.1[,-c(1:5)]) # each SNP belongs to only 1 bin 
  } else { LDSC.files <- list.files(out.path, pattern = "_ridge_constrained.gz", full.names = T) } 
  return(LDSC.files)
}


POLYFUN.merge_ldsc_files <- function(ldsc.files){
  priors <- lapply(ldsc.files, function(x){
    printer(x)
    pri <- data.table::fread(x, nThread = 4) 
    return(pri)
  }) %>% data.table::rbindlist()
  return(priors)
}


# %%%%%%%%%%%%%%%% Run PolyFun+SUSIE %%%%%%%%%%%%%%%% 
POLYFUN.SUSIE <- function(polyfun="./echolocatoR/tools/polyfun", 
                          dataset="Nalls23andMe_2019",
                          finemap_DT=NULL, 
                          polyfun_priors=c("precomputed","parametric","non-parametric"),
                          locus="_genome_wide"){
  
  # polyfun="./echolocatoR/tools/polyfun";  results_path="./Data/GWAS/Nalls23andMe_2019/_genome_wide"; dataset="Nalls23andMe_2019"; locus="LRRK2"; finemap_DT=NULL; polyfun_priors="parametric"
  results_path <- file.path(dirname(Directory_info("Nalls23andMe_2019")),locus)
  finemap_DT <- data.table::fread(file.path(results_path, "Multi-finemap/Multi-finemap_results.txt"))
  chrom <- unique(finemap_DT$CHR)
   
  # Import priors
  # ~~~~~~~~ Approach 1 ~~~~~~~~ 
  if (polyfun_priors=="precomputed"){
    priors <- POLYFUN.get_precomputed_priors(results_path=results_path,
                                             dataset=dataset,
                                             finemap_DT=finemap_DT,
                                             locus=locus)
    
  # ~~~~~~~~ Approach 2 ~~~~~~~~ 
  } else if (polyfun_priors=="parametric"){
    ldsc.files <- list.files(out.path, pattern = "*.snpvar_ridge_constrained.gz", full.names = T) %>% 
      grep(pattern = paste0(".",chrom,"."), value = T)
    priors <- POLYFUN.merge_ldsc_files(ldsc.files) 
    # ~~~~~~~~ Approach 3 ~~~~~~~~ 
  } else if (polyfun_priors=="non-parametric"){
    ldsc.files <- list.files(out.path, pattern = "*.snpvar_constrained.gz", full.names = T) %>% 
      grep(pattern = paste0(".",chrom,"."), value = T)
    priors <- POLYFUN.merge_ldsc_files(ldsc.files) 
  } 
  # Prepare data
  merged_DT <- data.table::merge.data.table(finemap_DT, 
                                            dplyr::select(priors, SNP, PolyFun.priors=SNPVAR) %>% 
                                              data.table::data.table(), 
                                            by="SNP")
  LD_matrix <- readRDS(file.path(results_path,"plink/LD_matrix.RData"))
  LD_matrix <- LD_matrix[merged_DT$SNP, merged_DT$SNP]
  # Run SUSIE
  subset_DT <- SUSIE(merged_DT, 
                     LD_matrix, 
                     dataset_type="GWAS",
                     n_causal=5,
                     sample_size=NA, 
                     var_y="estimate",
                     prior_weights=merged_DT$PolyFun.priors) 
  subset_DT <- subset_DT %>% dplyr::rename(PolyFun_SUSIE.Probability=Probability, 
                                           PolyFun_SUSIE.Credible_Set=Credible_Set) %>% 
    data.table::data.table()
  return(subset_DT)
}




# %%%%%%%%%%%%%%%% Run PolyFun+SUSIE %%%%%%%%%%%%%%%% 
POLYFUN.plot <- function(subset_DT,
                         locus=NULL,
                         subtitle="PolyFun Comparison",
                         conditions=c("SUSIE","PolyFun_SUSIE","FINEMAP","PAINTOR","PAINTOR_Fairfax")){
  # # Get r2 
  lead.snp <- subset(subset_DT, leadSNP)$SNP
  r2 <- data.table::data.table(SNP=names(LD_matrix[lead.snp,]), r2=LD_matrix[lead.snp,]^2)
  dat <- data.table::merge.data.table(subset_DT, r2, by="SNP")
  dat <- dplyr::mutate(dat, Mb=round(POS/1000000,3)) 
  
  
  library(patchwork)
  # GWAS
  gg <- ggplot(dat, aes(x=Mb, y=-log10(P), color=r2)) + 
    scale_color_gradient(low="blue",high="red", breaks=c(0,.5,1)) + 
    geom_point() + 
    labs(y="GWAS -log10(P)") + 
    # PolyFun priors
    ggplot(dat, aes(x=Mb, y=PolyFun.priors, color=PolyFun.priors)) + 
    scale_color_viridis_c(limits=c(0,1), breaks=c(0,.5,1)) +
    geom_point() + 
    # ylim(0,1) +
    # PolyFun+SUSIE PP
    ggplot(dat, aes(x=Mb, y=PolyFun_SUSIE.Probability, color=PolyFun_SUSIE.Probability)) + 
    geom_point() + 
    # SUSIE PP
    ggplot(dat, aes(x=Mb, y=SUSIE.Probability, color=SUSIE.Probability)) + 
    geom_point() + 
    # Overall layers
    patchwork::plot_layout(ncol = 1) + 
    patchwork::plot_annotation(title = "LRRK2", 
                               subtitle = subtitle,
                               theme =  theme(plot.title = element_text(hjust = 0.5),
                                              plot.subtitle = element_text(hjust = 0.5)))  
  print(gg)  
  ggsave(file.path(results_path,'PolyFun',"PolyFun.plot.png"), plot = gg,dpi = 400, height = 8, width = 7)
}





