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

 
# GitHub Notes:
## How to pull changes from original repo into the forked repo
# https://digitaldrummerj.me/git-syncing-fork-with-original-repo/


read_parquet <- function(parquet.file){
          # library(reticulate) 
          # conda_path="/hpc/packages/minerva-centos7/anaconda3/2018.12"
          reticulate::use_condaenv("polyfun_venv")
          # reticulate::conda_list() 
          pd <- reticulate::import("pandas")
          a_df <- pd$read_parquet(parquet.file)
          return(a_df) 
}

POLYFUN.help <- function(){
  cmd <- paste("python3",file.path(polyfun,"polyfun.py"),
               "--help")
  system(cmd)
}

POLYFUN.install_dependencies <- function(libraries = c("numpy",
                                                       "scipy",
                                                       "scikit-learn",
                                                       "pandas",
                                                       "tqdm",
                                                       "pyarrow",
                                                       "bitarray",
                                                       "networkx",
                                                       "rpy2")){  
  # Python
  system(paste("ml python/3.7.3 &&","pip install --user", paste(libraries,collapse=" "))) 
  system("pip install --user pandas && pip freeze | grep pandas")
  # R
  list.of.packages <- c("ggplot2", "Ckmeans.1d.dp","crayon")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages) 
}

 
POLYFUN.install_conda <- function(server=F){
  # Install anaconda
  # https://medium.com/ayuth/install-anaconda-on-macos-with-homebrew-c94437d63a37
  if(server){
    system("ml anaconda3")
  } else {
    system("brew cask install anaconda")
    system("conda init bash")
    try({system("conda init fish")}) 
    } 
  system("'export PATH='/usr/local/anaconda3/bin:$PATH' >> ~/.zshrc")
  system("source ~/.zshrc") 
}

POLYFUN.load_conda <- function(server=F){
  printer("POLYFUN:: Activating polyfun_venv...")
  if(server){ 
    reticulate::use_condaenv("polyfun_venv")
  }
  reticulate::use_condaenv("polyfun_venv",
                           conda = "/usr/local/anaconda3/condabin/conda")
  # conda_list("/usr/local/anaconda3/condabin/conda")
}


POLYFUN.conda_from_yaml <- function(yaml_path="./echolocatoR/tools/python_env.yml"){
  cmd <- paste("conda env create -f",yaml_path)
  print(cmd)
  system(cmd)
}

POLYFUN.conda_from_list <- function(libraries = c("numpy",
                                                  "scipy",
                                                  "scikit-learn",
                                                  "pandas",
                                                  "tqdm",
                                                  "pyarrow",
                                                  "bitarray",
                                                  "networkx",
                                                  "rpy2")){
  # NOTE: version specification must use quotes
  cmd <- paste("conda create -n polyfun_venv python=3.7.3",
               "numpy scipy scikit-learn 'pandas>=0.25.0'", 
               paste(libraries, collapse=" "))
  print(cmd)
  
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
  dataset <- basename(dirname(results_path))
  locus <- basename(results_path)
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
                                     min_MAF=0,
                                     server=F){  
  results_path <- file.path(dirname(Directory_info(dataset_name = dataset, "fullSS.local")), "_genome_wide")

  if(server){
    PF.output.path <- "/sc/orga/projects/pd-omics/tools/polyfun" 
  } else {
    PF.output.path <- file.path(results_path, "PolyFun")
  } 
  munged.path <- file.path(PF.output.path,
                           paste(dataset,".sumstats_munged.parquet",sep="")) 
  # data.table::fread("/sc/orga/projects/pd-omics/data/nallsEtAl2019/combined_meta/nallsEtAl2019_allSamples_allVariants.mod.txt", nrows = 2)
  fullSS.loc <- ifelse(server,"fullSS","fullSS.local")
  # Requires space-delimited file with the following columns (munging can recognize several variations of these names):
  ## SNP CHR BP ....and....
  ## either a p-value, an effect size estimate and its standard error, a Z-score or a p-value 
  if(!file.exists(munged.path)){ 
    printer("+ PolyFun:: Initiating data munging pipeline...")
    cmd <- paste("python", file.path(polyfun,"munge_polyfun_sumstats.py"),
                 "--sumstats",Directory_info(dataset_name = dataset, ifelse(server,"fullSS",fullSS.loc)),
                 "--n",sample.size, # Study sample size
                 "--out",munged.path,
                 "--min-info",0,#min_INFO,
                 "--min-maf", 0.001# min_MAF
                 )
    print(cmd)
    system(cmd)
  } else {printer("+ PolyFun:: Existing munged summary stats files detected.")} 
  return(munged.path)
}



POLYFUN.gather_ldscores <- function(output_prefix){
  ldscore.files <-  list.files(dirname(output_prefix), pattern = ".l2.ldscore.parquet", full.names = T)
  parquor <- POLYFUN.read_parquet(ldscore.files[1])
}
 
POLYFUN.gather_annotations <- function(chromosomes=c(1:22), 
                                       subset_SNPs=NULL,
                                       polyfun_annots="/sc/orga/projects/pd-omics/tools/polyfun/annotations/baselineLF2.2.UKB"){
  annot_DT <- lapply(chromosomes, function(chrom){
    printer("+ POLYFUN:: Chromosome",chrom,"...")
    parquet.file <- list.files(polyfun_annots, 
                               pattern = paste0("*\\.",chrom,"\\.annot\\.parquet"), 
                               full.names = T)
    annot_df <- read_parquet(parquet.file)
    annot_df$BP <- as.integer(annot_df$BP)
    if(!is.null(subset_SNPs)){
      annot_df <- subset(annot_df, SNP %in% subset_SNPs)
    }
    return(data.table::data.table(annot_df))
  }) %>% data.table::rbindlist() 
  return(annot_DT)
}


POLYFUN.read_parquet <- function(parquet_path){
  # Converts parquet to data.table
  SparkR::sparkR.session()
  parquor <- SparkR::read.parquet(parquet_path)
  parquor <- SparkR::as.data.frame(parquor) %>% 
    data.table::data.table() 
  return(parquor)
}
# 
# 
# POLYFUN.ukbb_LD <- function(finemap_DT, 
#                             results_path,
#                             force_new_LD=F){
#   base_url  <- "./echolocatoR/tools/polyfun/LD_temp"
#   alkes_url <- "https://data.broadinstitute.org/alkesgroup/UKBB_LD"
#   chr <- unique(finemap_DT$CHR)
#   min.pos <- min(finemap_DT$POS)
#   max.pos <- max(finemap_DT$POS)
#   file.name <- paste0("chr",chr,"_","40000001_43000001")
#   gz.path <- file.path(base_url,paste0(file.name,".gz"))
#   npz.path <- file.path(base_url,paste0(file.name,".npz"))
#   UKBB.LD.file <- file.path(results_path,"plink/ukbb_LD.RDS")
#   
#   # 1. Download LD files if not present 
#   # system(paste("wget",file.path(alkes_url,file.name)))
#   
#   # Download all UKBB LD files
#   # "wget -nc -r -A '*.gz|*.npz' https://data.broadinstitute.org/alkesgroup/UKBB_LD"
#   # "wget https://data.broadinstitute.org/alkesgroup/UKBB_LD/chr12_40000001_43000001.npz" # For some reason, way faster withOUT the --no-parent flag
#   # RSIDs file
#   # rsids <- data.table::fread(gz.path, nThread = 4) 
#   if(file.exists(UKBB.LD.file) & force_new_LD!=T){
#     printer("POLYFUN:: Pre-existing UKBB LD file detected.")
#     ld_R <- readRDS(UKBB.LD.file)
#   } else {
#     POLYFUN.load_conda()
#     reticulate::source_python(file.path(polyfun,"load_ld.py"))
#     ld.out <- load_ld(ld_prefix = file.path(base_url,"chr12_40000001_43000001"))
#     # LD matrix
#     ld_R <- ld.out[[1]]
#     head(ld_R)[1:10]
#     # SNP info
#     ld_snps <- ld.out[[2]] 
#     rsids <- subset(ld_snps, rsid %in% finemap_DT$SNP)
#     ld.indices <- as.integer(row.names(rsids))
#     ld_R <- ld_R[ld.indices, ld.indices]
#     row.names(ld_R) <- rsids$rsid
#     colnames(ld_R) <- rsids$rsid
#     printer("LD matrix dimensions", paste(dim(ld_R),collapse=" x "))
#     printer("+ POLYFUN:: Saving LD =>",UKBB.LD.file)
#     saveRDS(ld_R, UKBB.LD.file)
#   } 
#   return(ld_R)
# }

POLYFUN.download_ref_files <- function(alkes_url="https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase1_plinkfiles.tgz", 
                                       output_path="/sc/orga/projects/pd-omics/data/1000_Genomes/Phase1"){
  file_name <- basename(alkes_url)
  cmd <- paste("wget",
               alkes_url,
               "--no-parent", # Makes everything waayyyyyy faster
               "& mv",file_name, output_path,
               "& tar zxvf",output_path,"--strip 1")
  paste(cmd)   
  ref.prefix <- list.files(output_path, pattern = "*.bim", full.names = T)[1]
  ref.prefix <- gsub(".?.bim","", ref.prefix)
  return(ref.prefix)
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
                                    allow_missing_SNPs=T,
                                    ref.prefix="/sc/orga/projects/pd-omics/data/1000_Genomes/Phase1/1000G.mac5eur."){
  # Quickstart:
  # polyfun="./echolocatoR/tools/polyfun"; parametric=T;  weights.path=file.path(polyfun,"example_data/weights."); annotations.path=file.path(polyfun,"example_data/annotations."); munged.path= "./Data/GWAS/Nalls23andMe_2019/_genome_wide/PolyFun/sumstats_munged.parquet"; parametric=T; dataset="Nalls23andMe_2019"; prefix="PD_GWAS"; compute_ldscores=F; allow_missing_SNPs=T; chrom="all"; finemap_DT=NULL; locus="LRRK2"; server=F; ref.prefix="/sc/orga/projects/pd-omics/data/1000_Genomes/Phase1/1000G.mac5eur.";
 
  POLYFUN.load_conda(server = server)
  if(server){
    annotations.path <-  "/sc/orga/projects/pd-omics/tools/polyfun/annotations/baselineLF2.2.UKB/baselineLF2.2.UKB."
    weights.path <-  "/sc/orga/projects/pd-omics/tools/polyfun/annotations/baselineLF2.2.UKB/weights.UKB."
  } 
   
  # 0. Create paths
  results_path <- file.path(dirname(Directory_info(dataset_name = dataset, "fullSS.local")), "_genome_wide") 
  if(server){
    PF.output.path <- file.path("/sc/orga/projects/pd-omics/tools/polyfun")
  } else {
    PF.output.path <- file.path(results_path, "PolyFun")
  }
  dir.create(PF.output.path, showWarnings = F, recursive = T) 
  out.path <- file.path(PF.output.path,"output")
  output_prefix <- file.path(out.path, prefix, prefix)
  dir.create(out.path, showWarnings = F, recursive = T)
   
  
  # 1. Munge summary stats
  printer("PolyFun:: [1]  Create a munged summary statistics file in a PolyFun-friendly parquet format.")
  munged.path <- POLYFUN.munge_summ_stats(polyfun=polyfun,
                                           python = python,
                                           dataset="Nalls23andMe_2019",
                                           sample.size=1474097, 
                                           min_INFO = 0,
                                           min_MAF = 0.001, 
                                           server = server)   
  
  # 2. 
  ## If compute_ldscores == F:
  # This will create 2 output files for each chromosome: output/testrun.<CHR>.snpvar_ridge.gz and output/testrun.<CHR>.snpvar_ridge_constrained.gz. The first contains estimated per-SNP heritabilities for all SNPs (which can be used for downstream analysis with PolyFun; see below), and the second contains truncated per-SNP heritabilities, which can be used directly as prior causal probabilities in fine-mapping.
  # library(reticulate)
  # reticulate::use_virtualenv("polyfun_conda")
  # pd <- reticulate::import("pandas")
  # pd$read_csv("./Data/directories_table.csv")
  # reticulate::
  # source_python(file.path(polyfun,"polyfun.py"))
  
  # NOTE! if you're running without the "--no-partitions" flag, 
  ## you need to load R first `ml R`.
  printer("PolyFun:: [2] Run PolyFun with L2-regularized S-LDSC")
  cmd2 <- paste("python",file.path(polyfun,"polyfun.py"),
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
    # 3. Computationally intensive step
    printer("PolyFun:: [3] Compute LD-scores for each SNP bin")
    cmd3 <- paste("python",file.path(polyfun,"polyfun.py"),
                  "--compute-ldscores",
                  "--output-prefix",output_prefix,
                  "--bfile-chr",ref.prefix,
                  ifelse(chrom=="all","",paste("--chr",chrom)),
                  ifelse(allow_missing_SNPs,"--allow-missing","") )
    print(cmd3)
    system(cmd3) 
    # 4.
    printer("PolyFun:: [4] Re-estimate per-SNP heritabilities via S-LDSC")
    cmd4 <- paste("python",file.path(polyfun,"polyfun.py"),
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
    # The output of the PARTITIONED LDSC has the suffix: .snpvar_constrained.gz (one per chrom)
    LDSC.files <- list.files(out.path, 
                             pattern = "*.snpvar_constrained.gz", full.names = T)
    # pd_ldsc <- data.table::fread(PS_LDSC.files[1], nThread = 4) 
    # ldscore <- POLYFUN.read_parquet(file.path(out.path,"PD_GWAS.1.l2.ldscore.parquet")) 
    # bin.1 <- POLYFUN.read_parquet(file.path(out.path,"PD_GWAS.2.bins.parquet"))
    #rowSums(bin.1[,-c(1:5)]) # each SNP belongs to only 1 bin 
  } else { LDSC.files <- list.files(out.path, pattern = "_ridge_constrained.gz", full.names = T) } 
  return(LDSC.files)
}







#### Original LDSC (extended by PolyFun) 
POLYFUN.run_ldsc <- function(polyfun="./echolocatoR/tools/polyfun",
                             PF.output.path,
                             munged.path, 
                             min_INFO = 0.6,
                             min_MAF = 0.05,
                             annotations.path=file.path(polyfun,"example_data/annotations."),
                             weights.path=file.path(polyfun,"example_data/weights."), 
                             prefix="PD_GWAS_LDSC",
                             chrom="all",
                             compute_ldscores=F, 
                             allow_missing_SNPs=T,
                             munged_path="/sc/orga/projects/pd-omics/tools/polyfun/Nalls23andMe_2019.sumstats_munged.parquet",
                             ref.prefix="/sc/orga/projects/pd-omics/data/1000_Genomes/Phase1/1000G.mac5eur.",
                             freq.prefix="/sc/orga/projects/pd-omics/tools/polyfun/1000G_frq/1000G.mac5eur."){
  
  POLYFUN.load_conda(server = server)
  if(server){
    annotations.path <-  "/sc/orga/projects/pd-omics/tools/polyfun/annotations/baselineLF2.2.UKB/baselineLF2.2.UKB."
    weights.path <-  "/sc/orga/projects/pd-omics/tools/polyfun/annotations/baselineLF2.2.UKB/weights.UKB."
  } 
  
  # 0. Create paths
  results_path <- file.path(dirname(Directory_info(dataset_name = dataset, "fullSS.local")), "_genome_wide") 
  if(server){
    PF.output.path <- file.path("/sc/orga/projects/pd-omics/tools/polyfun")
  } else {
    PF.output.path <- file.path(results_path, "PolyFun")
  }
  dir.create(PF.output.path, showWarnings = F, recursive = T) 
  out.path <- file.path(PF.output.path,"output")
  output_prefix <- file.path(out.path, prefix, prefix)
  dir.create(out.path, showWarnings = F, recursive = T)
  # https://github.com/bulik/ldsc/wiki/Partitioned-Heritability
  cmd <- paste("python",file.path(polyfun,"ldsc.py"),
                "--h2",munged_path,
                "--ref-ld-chr",annotations.path,
                "--w-ld-chr",weights.path,
                "--overlap-annot",
                "--frqfile-chr",freq.prefix,
                "--not-M-5-50", # Important! enrichment estimates will be provided with MAF>0.1% SNPs instead of MAF>5% SNPs.
                "--out",output_prefix)
  # help_cmd <- paste("python",file.path(polyfun,"ldsc.py -h"))
  print(cmd)
  system(cmd)
  
}






# %%%%%%%%%%%%%%%% Run PolyFun+SUSIE %%%%%%%%%%%%%%%% 
POLYFUN.SUSIE <- function(results_path,
                          polyfun="./echolocatoR/tools/polyfun",
                          finemap_DT=NULL,
                          LD_matrix=NULL,
                          polyfun_approach="non-parametric",
                          dataset_type="GWAS",
                          n_causal=5, 
                          sample_size=NA){
  
  # polyfun="./echolocatoR/tools/polyfun";  results_path="./Data/GWAS/Nalls23andMe_2019/_genome_wide"; dataset="Nalls23andMe_2019"; locus="LRRK2"; finemap_DT=NULL; polyfun_priors="parametric"; sample.size=1474097; min_INFO=0; min_MAF=0; server=T;
  out.path <- file.path(dirname(results_path),"_genome_wide/PolyFun/output") 
  chrom <- unique(finemap_DT$CHR)

  
   
  # Import priors
  # ~~~~~~~~ Approach 1 ~~~~~~~~ 
  if (polyfun_approach=="precomputed"){
    priors <- POLYFUN.get_precomputed_priors(results_path=results_path, 
                                             finemap_DT=finemap_DT)
    
  # ~~~~~~~~ Approach 2 ~~~~~~~~ 
  } else if (polyfun_approach=="parametric"){
    ldsc.files <- list.files(out.path, pattern = "*.snpvar_ridge_constrained.gz", full.names = T) %>% 
      grep(pattern = paste0(".",chrom,"."), value = T)
    h2 <- rbind.file.list(ldsc.files) 
    # ~~~~~~~~ Approach 3 ~~~~~~~~ 
  } else if (polyfun_approach=="non-parametric"){
    ldsc.files <- list.files(out.path, pattern = "*.snpvar_constrained.gz", full.names = T) %>%  grep(pattern = paste0(".",chrom,"."), value = T)
    h2 <- rbind.file.list(ldsc.files) 
  } 
  # Prepare data
  merged_DT <- data.table::merge.data.table(finemap_DT, 
                                            dplyr::select(h2, SNP, POLYFUN.h2=SNPVAR) %>% 
                                              data.table::data.table(), 
                                            by="SNP")
  if(is.null(LD_matrix)){
    # LD_matrix <- readRDS(file.path(results_path,"plink/LD_matrix.RData"))
    LD_matrix <- POLYFUN.ukbb_LD(finemap_DT, 
                                 results_path,
                                 force_new_LD=F)
  } 
  LD_matrix <- LD_matrix[merged_DT$SNP, merged_DT$SNP]
  # Run SUSIE
  subset_DT <- SUSIE(merged_DT, 
                     LD_matrix, 
                     dataset_type=dataset_type,
                     n_causal=n_causal,
                     sample_size=sample_size, 
                     var_y="estimate",
                     prior_weights=merged_DT$POLYFUN.h2) 
  # subset_DT <- subset_DT %>% dplyr::rename(POLYFUN_SUSIE.PP=Probability, 
  #                                          PolyFun_SUSIE.Credible_Set=Credible_Set) %>%  data.table::data.table()
  return(subset_DT)
}


POLYFUN.finemapper <- function(polyfun= "./echolocatoR/tools/polyfun",
                               munged.path=NULL,
                               locus=NULL,
                               results_path=NULL){ 
  # finemap_DT <- quick_finemap();  
  base_url  <- "./echolocatoR/tools/polyfun/LD_temp"
  chr <- unique(finemap_DT$CHR) 
  file.name <- paste0("chr",chr,"_","40000001_43000001")
  ld_path <- file.path(base_url,file.name)
  gz.path <- file.path(base_url,paste0(file.name,".gz"))
  npz.path <- file.path(base_url,paste0(file.name,".npz"))
   
  
  cmd <- paste("python",file.path(polyfun,"finemapper.py"), 
                "--ld",ld_path,
                "--sumstats", Directory_info(dataset_name = dataset, ifelse(server,"fullSS","fullSS.local")),
                # "--sumstats","./echolocatoR/tools/polyfun/"
                "--n",1474097,
                "--chr",chr,
                "--start",min(finemap_DT$POS),
                "--end",max(finemap_DT$POS),
                "--method susie",
                "--max-num-causal 5",
                "--threads 2",# use max detected cores if not specified
                "--out",file.path(results_path,paste0("finemap.UKBB.",locus,".gz")))
  print(cmd)
  system(cmd)
    
  
}



# %%%%%%%%%%%%%%%% Run PolyFun+SUSIE %%%%%%%%%%%%%%%% 
POLYFUN.plot <- function(subset_DT,
                         LD_matrix,
                         locus=NULL,
                         subtitle="PolyFun Comparison",
                         conditions=c("SUSIE","POLYFUN_SUSIE","FINEMAP","PAINTOR","PAINTOR_Fairfax")){
  # Quickstart
  # locus="LRRK2"; subtitle="PolyFun Comparison"; conditions=c("SUSIE","POLYFUN_SUSIE","FINEMAP","PAINTOR","PAINTOR_Fairfax")
  # # Get r2 
  
  if(plot_ld){
    lead.snp <- top_n(subset_DT,1,-P)$SNP #subset(subset_DT, leadSNP==T)$SNP
    r2 <- data.table::data.table(SNP=names(LD_matrix[lead.snp,]), 
                                 r2=LD_matrix[lead.snp,]^2)
    dat <- data.table::merge.data.table(subset_DT, r2, by="SNP") 
  } else{dat <- dplyr::mutate(subset_DT, r2=1)}
  dat <- dplyr::mutate(dat, Mb=round(POS/1000000,3)) 
 
  
  
  library(patchwork)
  # GWAS
  gg <- ggplot(dat, aes(x=Mb, y=-log10(P), color=r2)) + 
    scale_color_gradient(low="blue",high="red", breaks=c(0,.5,1), limits=c(0,1)) + 
    geom_point() + 
    labs(y="GWAS -log10(P)") + 
    ggrepel::geom_label_repel(data = top_n(dat,n=1,-P), 
                              aes(label=SNP),alpha=0.7) +
    geom_point(data=top_n(dat,n=1,-P), size=5, shape=1, color="red") + 
    scale_y_continuous(limits = c(0,max(-log10(dat$P))*1.1)) + 
    
    # PolyFun priors
    ggplot(dat, aes(x=Mb, y=POLYFUN.h2, color=POLYFUN.h2)) + 
    scale_color_viridis_c(limits=c(0,1), breaks=c(0,.5,1)) +
    geom_point() + 
    # ylim(0,1) +
    
    # PolyFun+SUSIE PP
    ggplot(dat, aes(x=Mb, y=POLYFUN_SUSIE.PP, color=POLYFUN_SUSIE.PP)) + 
    geom_point() + 
    ggrepel::geom_label_repel(data = subset(dat,POLYFUN_SUSIE.PP>=.5), 
                              aes(label=SNP),alpha=0.7, color='green') +
    geom_point(data=subset(dat, POLYFUN_SUSIE.Credible_Set>0),
               #subset(dat, PolyFun_SUSIE.PP>=.95), 
               size=5, shape=1, color="green") + 
    scale_y_continuous(breaks = c(0,.5,1), limits = c(0,1.1)) + 
    scale_color_continuous(breaks=c(0,.5,1), limits=c(0,1)) +
    
    # SUSIE PP
    ggplot(dat, aes(x=Mb, y=SUSIE.PP, color=SUSIE.PP)) + 
    geom_point() + 
    ggrepel::geom_label_repel(data = subset(dat,SUSIE.Credible_Set>0), 
                              aes(label=SNP),alpha=0.8, color='green') +
    geom_point(data=subset(dat,SUSIE.PP>=.95), 
               size=5, shape=1, color="green") + 
    scale_y_continuous(breaks = c(0,.5,1), limits = c(0,1.1)) + 
    scale_color_continuous(breaks=c(0,.5,1), limits=c(0,1)) +
    # Overall layers
    patchwork::plot_layout(ncol = 1) + 
    patchwork::plot_annotation(title = locus, 
                               subtitle = paste(nrow(dat),"SNPs"),#"PolyFun Comparison",
                               theme =  theme(plot.title = element_text(hjust = 0.5),
                                              plot.subtitle = element_text(hjust = 0.5)))  
  print(gg)  
  ggsave(file.path(results_path,'PolyFun',"PolyFun.plot.png"), plot = gg,
         dpi = 400, height = 10, width = 7)
}
 
















GGBIO.ucsc_tracks <- function(finemap_DT){ 
  # GLASS DATA: UCSC GB
  # https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr2:127770344-127983251&hgsid=778249165_ySowqECRKNxURRn6bafH0yewAiuf
  
  
  # finemap_DT <- quick_finemap()
  # UCSC Tracks   
  import.bw.filt <- function(bw.file, gr.dat){  
    bw.dat <- rtracklayer::BigWigSelection(ranges = gr.dat, 
                                           colnames = "score")  
    # import.bw:: https://rdrr.io/bioc/Gviz/src/R/Gviz.R#sym-.import.bw
    bw.filt <- rtracklayer::import.bw(con = bw.file, selection= gr.dat)
    # bw.filt <- Gviz:::.import.bw(bw.file, selection = gr.dat)
    # gr <- bw.filt@range
    # gr.filt <- gr[seqnames(gr) == chr & start(gr) > from & end(gr) < to] 
    plot(x = start(bw.filt), y=bw.filt$score)
    return(bw.filt)
  }
  # Import BigWig annotation files
  bigWigFiles <- readxl::read_excel("./echolocatoR/annotations/Glass_lab/Glass.snEpigenomics.xlsx") 
  bigWigFiles <- subset(bigWigFiles, marker!="-" &  cell_type!="peripheral microglia")
  
  # Convert finemap data to granges
  dat <- finemap_DT
  dat$seqnames <- paste0("chr",dat$CHR) 
  dat$start.end <- dat$POS
  gr.dat <- GenomicRanges::makeGRangesFromDataFrame(df = dat,
                                                    seqnames.field = "seqnames",
                                                    start.field = "start.end", 
                                                    end.field = "start.end", 
                                                    keep.extra.columns = T)   
  bw.grlist <- lapply(1:nrow(bigWigFiles), function(i){ 
    bw.file <- bigWigFiles$data_link[i]
    bw.name <- gsub("_pooled|pooled_","",bigWigFiles$name[i])
    printer("GVIZ:: Importing...",bw.name)
    bw.filt <- import.bw.filt(bw.file=bw.file, 
                              gr.dat = gr.dat)
    colnames(mcols(bw.filt))[1] <- bw.name
    # bw.filt$expt_name <- bw.name  
    # bw.filt$cell_type <-strsplit(bw.name, "_")[[1]][[1]]
    # bw.filt$assay <- strsplit(bw.name, "_")[[1]][[2]]
    return(bw.filt)
  })
  gr.snp <- Reduce(function(x, y) GenomicRanges::merge(x, y, all.x=T), 
                   append(bw.grlist, gr.dat))  
  
 
  
  # Fine-mapping tracks
  INIT_list <- list(   
    "GWAS"=ggbio::plotGrandLinear(obj = gr.snp, aes(y=-log10(P), fill=-log10(P))),
    "Fine_mapping"=ggbio::plotGrandLinear(obj = gr.snp, aes(y=mean.PP, color=mean.PP)) + 
      scale_color_viridis_c()
  )    
  # BigWig tracks from Glass
  bw.cols <- colnames(mcols(gr.snp))[1:12]  
  BW_tracks<- y(bw.cols, function(annot){
    gg.bw <- ggbio::plotGrandLinear(gr.snp, aes(y=eval(parse(text=annot)), 
                                                color=eval(parse(text=annot))) 
    ) + labs(y="Score", color="Score") + ylim(0,45)
    # if(bw.cols!=bw.cols[1]){
    #   gg.bw <- gg.bw + labs()
    # }
    return(gg.bw)                                                             
  }) 
  TRACKS_list <- append(INIT_list, BW_tracks)  
  names(TRACKS_list) <- c(names(INIT_list), gsub("_","\n",bw.cols))
  
  
  # gr.snp[start(gr.snp)!=Inf]
  # Fuse all tracks 
  gene <- "LRRK2"
  library(ggbio)
  params_list <- list(title = paste0(gene," [",length(seqnames(gr.snp))," SNPs]"), 
                      track.bg.color = "transparent",
                      track.plot.color = "transparent",
                      label.text.cex = .7, 
                      label.bg.fill = "grey12",
                      label.text.color = "white",
                      label.text.angle = 0,
                      label.width = unit(5.5, "lines"),
                      # xlim = c(min(start(gr.snp)), max(start(gr.snp))),
                      heights = c(rep(1,length(INIT_list)), rep(1,length(BW_tracks)) ))
  TRACKS_list <- append(TRACKS_list, params_list)
  trks <- suppressWarnings(do.call("tracks", TRACKS_list)) 
  
  # add lines
  lead.pos <- subset(finemap_DT,SNP=="rs76904798")$POS
  consensus.pos <- subset(finemap_DT, Consensus_SNP==T)$POS
  trks_plus_lines <- trks + 
    # theme_bw() + 
    # ggbio::theme_genome() +  
    # theme(strip.text.y = element_text(angle = 0), 
    #       strip.text = element_text(size=9 )) + 
    geom_vline(xintercept = lead.pos, color="red", alpha=1, size=.3, linetype='solid') +
    geom_vline(xintercept = consensus.pos, color="goldenrod2", alpha=1, size=.3, linetype='solid') 
   
  ggsave(file.path(results_path,"Annotation","Glass.snEpigenomics.png"), 
         plot = trks_plus_lines, dpi=400, height = 15, width = 8, 
         bg = "transparent")
  trks_plus_lines
  
  INIT_list[[1]] +
    geom_vline(xintercept = , color="red", alpha=1, size=.3, linetype='solid')   
   
  # Save merged Glass epigenomic data
  data.table::fwrite(GenomicRanges::as.data.frame(gr.snp[,bw.cols]),
                     file.path("./echolocatoR/annotations/Glass_lab/LRRK2.Glass.txt"),
                     sep="\t", nThread = 4)
}


NOTT_2019.superenhancers <- function(s6_path="./echolocatoR/annotations/Glass_lab/aay0793-Nott-Table-S6.xlsx"){ 
  s6 <- readxl::read_excel( , skip = 2)  
  annot_sub <- subset(s6, chr== paste0("chr",unique(finemap_DT$CHR)) & start>=min(finemap_DT$POS) & end<=max(finemap_DT$POS) )
  if(nrow(annot_sub)>0){
    merged_DT <- data.table::merge.data.table(finemap_DT %>% 
                                                dplyr::mutate(chr=paste0("chr",CHR), 
                                                              start=as.numeric(POS)) %>% 
                                                data.table::data.table(),
                                              data.table::data.table(s6), 
                                              by = c("chr","start"))
  }
}
  
NOTT_2019.enhancer_promoter_interactions <- function(finemap_DT, 
                                                     s5_path="./echolocatoR/annotations/Glass_lab/aay0793-Nott-Table-S5.xlsx"){  
  sheets_s5 <- readxl::excel_sheets(s5_path) 
  s5 <- readxl::read_excel(s5_path, sheet = sheets_s5[1], skip = 2) 
  s5 <- s5 %>% dplyr::rename(chr=Chr, start=Start, end=End)  
  # Subset to window
  annot_sub <- subset(s5, chr== paste0("chr",unique(finemap_DT$CHR)) & 
                        start>=min(finemap_DT$POS) & 
                        end<=max(finemap_DT$POS) )  
  return(annot_sub)
}

NOTT_2019.get_promoters <- function(annot_sub){ 
  promoter.cols <- grep("*_active_promoter",  colnames(annot_sub), value = T)
  promoter_celltypes <- gsub("\\_.*","", promoter.cols[as.logical(annot_sub[,promoter.cols])] )
  promoter_celltypes <- as.character(marker_key[promoter_celltypes])
  promoter_celltypes <- paste(promoter_celltypes,collapse="; ")
  return(promoter_celltypes)
}

NOTT_2019.get_interactions <- function(annot_sub){  
  interact.cols <- grep("*_interactions", colnames(annot_sub), value = T)   
  interact.DT <- lapply(interact.cols, function(column){
    coords <- strsplit(annot_sub[,column][[1]], ",")
    coord.dt <- lapply(coords, function(coord){
      data.table::data.table(Interaction=column, 
                             Cell_type=marker_key[gsub("\\_.*","",column)],
                             Coordinates=coord) %>% return() 
    }) %>% data.table::rbindlist()
    return(coord.dt)
  } )  %>% data.table::rbindlist()
  interact.DT <- subset(interact.DT, !is.na(Coordinates)) %>%  
    tidyr::separate(col = Coordinates, 
             into=c("chr","Start","End"), sep = ":|-") %>%
     separate(col = Interaction, into=c("Marker","Element",NA), sep="_", remove = F ) 
  interact.DT <- interact.DT %>% 
    dplyr::mutate(Cell_type_interaction=paste(Cell_type,"-",Element))
  interact.DT$Cell_type <- interact.DT$Cell_type %>% as.character()
  interact.DT$Start <- as.numeric(interact.DT$Start)
  interact.DT$End <- as.numeric(interact.DT$End)
  # Summarise distance from different celltype enhancer interactions
  summarise_top.consensus.dist <- interact.DT %>% 
    dplyr::mutate(top.consensus.dist=End - top.consensus.pos) %>% 
    dplyr::group_by(Cell_type) %>% 
    dplyr::summarise(top.consensus.dist = mean(top.consensus.dist))
  print(summarise_top.consensus.dist)
  return(interact.DT)
}

GGBIO.nott_etal_2019 <- function(){
  # finemap_DT <- quick_finemap()
  library(ggbio)
  marker_key <- list(PU1="Microglia",
                     Olig2="Oligodendrocytes",
                     NeuN="Neurons",
                     LHX2="Astrocytes")
  lead.pos <- subset(finemap_DT,SNP=="rs76904798")$POS  
  consensus.pos <- subset(finemap_DT, Consensus_SNP==T)$POS  
  top.consensus.pos <- top_n(subset(finemap_DT, Consensus_SNP==T), n=1, wt = mean.PP)$POS  
  # Subset to relevant region
  annot_sub <- NOTT_2019.enhancer_promoter_interactions(finemap_DT = finemap_DT)
  ## Extract active promoters
  promoter_celltypes <- NOTT_2019.get_promoters(annot_sub)
  ## Extract promoter interactions
  interact.DT <- NOTT_2019.get_interactions(annot_sub)
  
   
  # GWAS track
  GWAS_trk <- ggplot(data=finemap_DT, aes(x=POS, y=-log10(P), color=-log10(P))) +
    geom_point()  
  
  FM_trk <- ggplot(data=finemap_DT, aes(x=POS, y=mean.PP, color=mean.PP)) +
    geom_point() +  
    scale_color_viridis_c(breaks=c(0,.5,1), limits=c(0,1)) +
    ylim(0,1)
    
  # Nott tracks   
  # Nott_s5 <- ggbio::ggbio() + 
  #   ggbio::geom_rect(data = annot_sub, aes(xmin=start, xmax=end, ymin=0, ymax=1), 
  #                    fill="turquoise", alpha=.75) +   
  #   facet_grid(facets = Annotation~.) 
  # Nott_s5 <- invisible_legend(Nott_s5)
  
  # Nott:  interactions
  # Nott_interactions
  
  NOTT.interact_trk <- ggplot() +
     ggbio::geom_arch(data = interact.DT, aes(x=Start, xend=End, color=Cell_type_interaction)) + 
    scale_y_reverse() +
    scale_colour_brewer(palette = "Accent") + 
    labs(subtitle = paste0(annot_sub$Annotation[[1]]," - ",promoter_celltypes) ) 
    
  
  # Makes tracks list
  TRACKS_list <- list(
      "GWAS"=GWAS_trk,
      "Fine-mapping"=FM_trk,
      "Nott (2019)\nInteractome"=NOTT.interact_trk
      # "Nott et al. (2019)\nPromoter\ninteractome"=Nott_s5 
  )
  # Parameters
  gene <- "LRRK2"
  params_list <- list(title = paste0(gene), 
                      track.bg.color = "transparent",
                      track.plot.color = "transparent",
                      label.text.cex = .7, 
                      label.bg.fill = "grey12",
                      label.text.color = "white",
                      label.text.angle = 0,
                      label.width = unit(5.5, "lines"),
                      xlim = c(min(finemap_DT$POS), max(finemap_DT$POS))
                      # heights = c(rep(1,length(INIT_list)), rep(1,length(BW_tracks)) )
                      )
  TRACKS_list <- append(TRACKS_list, params_list)
  trks <- suppressWarnings(do.call("tracks", TRACKS_list)) 
   
  trks_plus_lines <- trks + 
    # Nott: promoter
    ggbio::geom_rect(data = annot_sub, aes(xmin=start, xmax=end, ymin=-0, ymax=Inf), 
                     fill="turquoise", alpha=.5, inherit.aes=F) +  
    # Lead GWAS line
    geom_vline(xintercept = lead.pos, color="red", alpha=1, size=.3, linetype='solid') +
    # Consensus line
    geom_vline(xintercept = consensus.pos, color="goldenrod2", alpha=1, size=.3, linetype='solid') + theme_classic() +
    theme(plot.subtitle = element_text(color = "turquoise", size = 8))
  
  ggsave(filename = file.path(results_path,"Annotation",paste0("Nott.sn-epigenomics_ggbio.png")), 
         plot = trks_plus_lines,
         height = 7, width = 7, dpi = 1000, bg = "transparent") 
  
  
  }
  

  


###############################################
########## POLYFUN H2 ENRICHMENT ##############
###############################################

POLYFUN.ldsc_annot_enrichment <- function(.results = "Data/GWAS/Nalls23andMe_2019/_genome_wide/PolyFun/output/PD_GWAS_LDSC/PD_GWAS_LDSC.results",
                                           show_plot=T,
                                           save_plot=F, 
                                           title = "LDSC Heritability Enrichment", 
                                           subtitle = "PD GWAS"){  
  res <- data.table::fread(.results)
  res$Category <- gsub("*_0$","", res$Category)
  res$Group <- gsub("^([^_]*_[^_]*)_.*$", "\\1", res$Category)
  # Get rid of absurd enrichment values
  res <- subset(res, Enrichment>-5000 & Enrichment<5000)
  res$Valence <- ifelse(res$Enrichment>=0,1,-1)
  res$p_adj <- p.adjust(res$Enrichment_p, method="fdr")
  
  
  POLYFUN.annot_enrichment_plot  <- function(res, 
                                             title = "LDSC Heritability Enrichment", 
                                             subtitle = "PD GWAS"){
    sig_res <- subset(res, p_adj<0.05)
    nudge_x <- ifelse(nrow(res)>150, -.5, -.03) 
    gp <- ggplot(res, aes(x=Category, y=Enrichment, fill=Group)) + 
      geom_col(show.legend = F) + 
      geom_errorbar(aes(ymin=Enrichment-Enrichment_std_error, ymax=Enrichment+Enrichment_std_error), 
                    width=.5, position=position_dodge(.9)) + 
      geom_text(data = sig_res, aes(x=Category, y=(Enrichment+Enrichment_std_error*Valence)+5*Valence), 
                label="*", nudge_x=nudge_x, color="magenta", size=7) +  
      coord_flip() +  
      labs(title = title, 
           subtitle = subtitle) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5))
    if(nrow(res)>150){ gp <- gp + theme(axis.text.y=element_text(size=7)) }
    return(gp)
  }

  if(show_plot){ 
    hist(res$Enrichment, breaks = 50)
    hist(res$Enrichment_p, breaks = 50)
    plot_dir <- dirname(.results)
    # Alphabetical
    gp.all <- POLYFUN.annot_enrichment_plot(res, title, subtitle)
    if(show_plot){print(gp.all)}
    if(save_plot){ggsave(gp.all, filename = file.path(plot_dir,"ldsc_annot_enrich_all.png"), 
                         dpi = 400, width=10, height=20)}
    # Top p-vals 
    gp.sig <- POLYFUN.annot_enrichment_plot(subset(res, p_adj<0.05), title, subtitle)
    if(show_plot){print(gp.sig)}
    if(save_plot){ggsave(gp.sig, filename = file.path(plot_dir,"ldsc_annot_enrich_sig.png"), 
                         dpi = 400, width=10,height=5)} 
  }
  
  
  POLYFUN.get_annot_refs <- function(res, 
                                     supp_file="./echolocatoR/tools/polyfun/SuppTables.xlsx", sheet="S1"){
    supp <- readxl::read_excel(supp_file,sheet = sheet)
    supp_merge <- data.table:::merge.data.table(data.table::data.table(supp), res,
                                  by.x = "Annotation", by.y = "Category") 
    return(supp_merge)
  }
  # supp_merge <- POLYFUN.get_annot_refs(sig_res)
  # createDT(supp_merge[,c("Annotation","Reference","Prop._SNPs","Prop._h2")])
 return(res)
}


POLYFUN.h2_enrichment <- function(h2_df, 
                                  target_SNPs=NULL){  
  # Only consider SNPs that overlap between LDCS and GWAS results to make things fair 
  # target_SNPs <- intersect(target_SNPs, h2_df$SNP)
  h2.target <- subset(h2_df, SNP %in% target_SNPs)
  if(nrow(h2.target)>0){
    # Calculate enrichment
    target_h2 <- sum(h2.target$SNPVAR, na.rm = T)
    total_h2 <- sum(h2_df$SNPVAR, na.rm = T)
    n_target_SNPs <- nrow(h2.target)
    n_total_SNPs <- nrow(h2_df)
    
    h2.enrichment <- (target_h2/total_h2) / (n_target_SNPs/n_total_SNPs)
  } else {
    h2.enrichment <- NA
  } 
  return(h2.enrichment)
}


POLYFUN.h2_enrichment_SNPgroups <- function(finemap_DT, 
                                            chrom="*", 
                                            ldsc_suffix="*.snpvar_constrained.gz",
                                            subtitle="", 
                                            show_plot=T,
                                            save_plot="./h2_enrichment.png",
                                            out.path="./Data/GWAS/Nalls23andMe_2019/_genome_wide/PolyFun/output"){
  # Quickstart
  # finemap_DT=quick_finemap(); chrom="*"; ldsc_suffix="*.snpvar_constrained.gz"; subtitle="";show_plot=T; save_plot="./h2_enrichment.png"; out.path="./Data/GWAS/Nalls23andMe_2019/_genome_wide/PolyFun/output"; conda_path="/hpc/packages/minerva-centos7/anaconda3/2018.12"; polyfun_annots="/sc/orga/projects/pd-omics/tools/polyfun/annotations/baselineLF2.2.UKB"
  # finemap_DT <- merge_finemapping_results(minimum_support = 0)
  # Gather your heritability
  ldsc.files <- list.files(out.path, pattern = ldsc_suffix, full.names = T) %>% 
    grep(pattern = paste0(".",chrom,"."), value = T) 
  h2_df <- rbind.file.list(ldsc.files)  
  # Subset to only common snps: h2 vs. finemap files
  # common.snps <- intersect(h2_df$SNP, finemap_DT$SNP)
  # h2_sub <- subset(h2_df, SNP %in% common.snps) 
  # finemap_DT <- subset(finemap_DT, SNP %in% common.snps) 
  
  
  
  # GWAS nominally sig hits
  GWAS.nom.sig <- POLYFUN.h2_enrichment(h2_df=h2_df, 
                                        target_SNPs=subset(finemap_DT, P<.05)$SNP )
  # GWAS sig hits
  GWAS.sig <- POLYFUN.h2_enrichment(h2_df=h2_df, 
                                    target_SNPs=subset(finemap_DT, P<5e-8)$SNP)
  # Credible Set
  Finemap.credset <- POLYFUN.h2_enrichment(h2_df=h2_df, 
                                           target_SNPs = subset(finemap_DT, Support>0)$SNP)
  # PAINTOR CS
  PAINTOR.credset <- POLYFUN.h2_enrichment(h2_df=h2_df, 
                                           target_SNPs = subset(finemap_DT, PAINTOR.Credible_Set>0)$SNP)
  FINEMAP.credset <- POLYFUN.h2_enrichment(h2_df=h2_df, 
                                           target_SNPs = subset(finemap_DT, FINEMAP.Credible_Set>0)$SNP)
  SUSIE.credset <- POLYFUN.h2_enrichment(h2_df=h2_df, 
                                           target_SNPs = subset(finemap_DT, FINEMAP.Credible_Set>0)$SNP)
  
  # Consenus SNP
  Finemap.consensus <- POLYFUN.h2_enrichment(h2_df=h2_df,
                                             target_SNPs = subset(finemap_DT, Support>1)$SNP)
  
  res <- data.frame(SNP.Group=c("GWAS_nom. sig.",
                                "GWAS_sig.",
                                "Fine-mapped_Credible Set",
                                "PAINTOR_Credible Set",
                                "FINEMAP_Credible Set",
                                "SUSIE_Credible Set",
                                "Fine-mapped_Consensus"),
                    h2.enrichment=c(GWAS.nom.sig,
                                    GWAS.sig, 
                                    Finemap.credset, 
                                    PAINTOR.credset,
                                    FINEMAP.credset,
                                    SUSIE.credset, 
                                    Finemap.consensus))
  
  if(show_plot){
    res$SNP.Group <- factor(gsub("_","\n",res$SNP.Group), 
                            levels = unique(gsub("_","\n",res$SNP.Group)))
    gp <- ggplot(data = res, aes(x=SNP.Group, y=h2.enrichment, fill= SNP.Group)) + 
      geom_col(show.legend = F) + 
      labs(title="PolyFun Heritability Enrichment",
           subtitle=subtitle,
           x="SNP Group") 
    print(gp)
    if(save_plot!=F){
      ggsave(plot = gp, filename = save_plot, width = 7, height=7)
    }
  } 
  return(res)
}
# merged_results <- merge_finemapping_results(minimum_support=0,
#                                             include_leadSNPs=T)
merged_results <- readxl::read_excel("./Data/annotated_finemapping_results.xlsx")
h2.enrich.df <- POLYFUN.h2_enrichment_SNPgroups(finemap_DT = merged_results,
                                         subtitle = "Genome-wide")
h2.enrich.annot.df <- POLYFUN.h2_enrichment_annot(finemap_DT = merged_results,
                                                  h2_df = h2_df)


# 
# POLYFUN.h2_enrichment_annot <- function(finemap_DT, 
#                                         h2_df,
#                                         annot_DT, 
#                                         locus=F){ 
#   if(locus!=F){ 
#     printer(locus)
#     finemap_DT <- subset(finemap_DT, Gene==locus) 
#   }
#   # Subset to only common snps: h2 vs. finemap files vs. annotation files
#   common.snps <- intersect(intersect(h2_df$SNP, finemap_DT$SNP), annot_DT$SNP)
#   annot_DT <- subset(annot_DT, SNP %in% common.snps)
#   finemap_DT <- subset(finemap_DT, SNP %in% common.snps)
#   h2_df.sub <- subset(h2_df, SNP %in% common.snps)
#   printer("+ Total number of SNPs:",nrow(h2_df.sub))
#   
#   # Calculate enrichment for each annotation
#   annot.names <- colnames(annot_DT[,-c(1:5)])
#   printer("+ POLYFUN:: Calculating enrichment for",length(annot.names),"annotations.")
#   annot.enrich <- lapply(annot.names, function(annot){
#     target_SNPs <- annot_DT[annot_DT[[annot]] == 1,]$SNP  
#     h2.enrich <- POLYFUN.h2_enrichment(h2_df = h2_df.sub, 
#                                        target_SNPs = target_SNPs)
#     res.dt <- data.table::data.table(Annotation=annot, h2.enrichment=h2.enrich)
#     return(res.dt)
#   }) %>% data.table::rbindlist() 
#   
#   return(annot.enrich)
# }


# 
# POLYFUN.h2_enrichment_annot_SNPgroups <- function(finemap_DT, 
#                                                   chrom="*", 
#                                                   ldsc_suffix="*.snpvar_ridge_constrained.gz",
#                                                   subtitle="", 
#                                                   show_plot=T,
#                                                   save_plot="./h2_enrichment_annotations.png",
#                                                   out.path="./Data/GWAS/Nalls23andMe_2019/_genome_wide/PolyFun/output"){
#   # Quickstart
#   # finemap_DT=quick_finemap(); chrom="*"; ldsc_suffix="*.snpvar_ridge_constrained.gz"; subtitle="";show_plot=T; save_plot="./h2_enrichment.png"; out.path="./Data/GWAS/Nalls23andMe_2019/_genome_wide/PolyFun/output"; conda_path="/hpc/packages/minerva-centos7/anaconda3/2018.12"; polyfun_annots="/sc/orga/projects/pd-omics/tools/polyfun/annotations/baselineLF2.2.UKB"
#   # Gather your heritability
#   ldsc.files <- list.files(out.path, pattern = ldsc_suffix, full.names = T) %>% 
#     grep(pattern = paste0(".",chrom,"."), value = T)
#   h2_df <- POLYFUN.merge_ldsc_files(ldsc.files) 
#   #
#   printer("POLYFUN:: Gathering binary annotatin files...")
#   chromosomes <- unique(finemap_DT$CHR)
#   annot_DT <- POLYFUN.gather_annotations(chromosomes = chromosomes, 
#                                          subset_SNPs = finemap_DT$SNP)
#   
#   # Calculate enrichment for each annotation
#   locus=F#"LRRK2"
#   GWAS.nom.sig <- POLYFUN.h2_enrichment_annot(subset(finemap_DT, P<0.05), 
#                                               h2_df,
#                                               annot_DT,
#                                               locus)
#   GWAS.nom.sig$SNP.Group <- "GWAS_nom. sig."
#   
#   GWAS.sig <- POLYFUN.h2_enrichment_annot(subset(finemap_DT, P<5e-8), 
#                                           h2_df,
#                                           annot_DT,
#                                           locus)
#   GWAS.sig$SNP.Group <- "GWAS_sig."
#   
#   Finemap.credset <- POLYFUN.h2_enrichment_annot(subset(finemap_DT, Support>0), 
#                                                  h2_df,
#                                                  annot_DT,
#                                                  locus)
#   Finemap.credset$SNP.Group <- "Fine-mapped_Credible Set"
#   
#   
#   Finemap.consensus <- POLYFUN.h2_enrichment_annot(subset(finemap_DT, Consensus_SNP==T), 
#                                                    h2_df,
#                                                    annot_DT,
#                                                    locus)
#   Finemap.consensus$SNP.Group <- "Fine-mapped_Consensus"
#   
#   res <- data.table::rbindlist(list(GWAS.nom.sig, GWAS.sig, Finemap.credset, Finemap.credset, Finemap.consensus)) 
#   
#   if(show_plot){
#     library(ggplot2)
#     gp <- ggplot(data = res, aes(x= Annotation, y=h2.enrichment, fill= SNP.Group)) + 
#       facet_grid(facets = SNP.Group~.) + 
#       geom_col(show.legend = F, position = "dodge") + 
#       labs(title="PolyFun-LDSC Enrichment",
#            subtitle=subtitle,
#            x="SNP Group") +
#       theme_classic() + 
#       theme(axis.text.x = element_text(angle = 45, hjust = 1, size=4))
#     print(gp)
#     if(save_plot!=F){
#       ggsave(plot = gp, filename  = save_plot, height=10, width=5)
#     }
#   }
# }
# 
# POLYFUN.finemapped_traits <- function(finemap_DT){  
#   fm_traits <- readxl::read_excel("./echolocatoR/tools/polyfun/SuppTables.xlsx", sheet = "S7")
#   fm_traits <- fm_traits %>% 
#     dplyr::mutate(POS=as.integer(BP)) %>%  
#     data.table::data.table() %>% 
#     subset(CHR %in% unique(finemap_DT$CHR) & POS >= min(finemap_DT$POS) & POS <= max(finemap_DT$POS))
#   fm_traits$`PIP (PolyFun + SuSiE)` <- as.numeric(gsub("<|>","",fm_traits$`PIP (PolyFun + SuSiE)`))
#   finemap_merged <- data.table:::merge.data.table(fm_traits, finemap_DT, 
#                                                   by=c("CHR","POS","SNP"), 
#                                                   all.y = T) 
#   library(patchwork)
#   ggplot(finemap_merged) + 
#     geom_point(aes(x=POS, y=mean.PP, color=mean.PP)) +  
#   ggplot(finemap_merged) +
#     geom_point(aes(x=POS, y=PolyFun_SUSIE.PP, 
#                    color=PolyFun_SUSIE.PP)) + 
#     patchwork::plot_layout(ncol = 1)
#   
# }




# 
# 
# GVIZ.ucsc_tracks <- function(){
#   # https://bioconductor.org/packages/release/bioc/html/Gviz.html
#   # Tutorial
#   # https://bioconductor.org/packages/release/bioc/vignettes/Gviz/inst/doc/Gviz.html
#   # Gviz::UcscTrack() 
#   library(Gviz)
#   library(GenomicRanges)
#   
#   subset_DT <- finemap_DT
#   from <- min(subset_DT$POS)
#   to <- max(subset_DT$POS)
#   chr <- paste0("chr",subset_DT$CHR[1])
#   ucsc_coords <- paste0(chr,":",from,"-",to)
#   gen <- "hg19" # "mm9"#
#    
#   # Finemapping Data track
#   d.track <- DataTrack(data=-log10(subset_DT$P),
#                        start=subset_DT$POS, end=subset_DT$POS-1, chr=chr,  
#                        genome=gen, name="GWAS", col="purple",
#                        type=c("p"))
#   # plotTracks(d.track, from=from, to=to, chr=chr) 
#   # availableDisplayPars(dtrack)
#   
#   # command above takes some time to fetch the data...
#   # plotTracks(nuc_track, from=from, to=to, genome=gen, chromosome=chr)
#   
#   # Data
#   # atrack <- AnnotationTrack(cpgIslands, name="CpG") 
#   # Axis
#   g.track <- GenomeAxisTrack()
#   # Gene Models
#   data(geneModels) 
#   gr.track <- GeneRegionTrack(geneModels, from=from, to=to,
#                               genome=gen, chromosome=chr, 
#                               name="Gene Model", transcriptAnnotation="symbol", 
#                               background.title="grey20") 
#   # Ideogram
#   i.track <- IdeogramTrack(genome=gen, chromosome=chr)
#   # Genomic Sequence
#   library(BSgenome.Hsapiens.UCSC.hg19)
#   s.track <- SequenceTrack(Hsapiens, chromosome=chr,name="DNA Sequence")  
#   
#   
#   # Plot
#   plotTracks(list(i.track, 
#                   g.track,  
#                   d.track, 
#                   nuc.track, 
#                   gr.track, 
#                   s.track),  
#              background.title="grey20", 
#              background.panel="transparent"
#              # rotation.title=0
#              )
# }




# rtracks <- function( ){
#   # https://www.bioconductor.org/packages/release/bioc/vignettes/rtracklayer/inst/doc/rtracklayer.pdf
#   library(rtracklayer)
#   library(GenomicRanges)
#   
#   mySession <- browserSession ()
#   genome(mySession) <- "hg19" 
#   track.names <- trackNames(ucscTableQuery(mySession)) 
#   df <- data.frame(track.names)
#   
#   
#   track.name <- "wgEncodeUwDgf"  
#   e2f3.grange <- GRanges("chr6", IRanges(20400587, 20403336)) 
#   mySession <- browserSession() 
#   tbl.k562.dgf.e2f3 <- getTable(ucscTableQuery(mySession, track=track.name,
#                                                range=e2f3.grange, table=table.name))
#   
#   tbl.k562.dgf.hg19 <- getTable(ucscTableQuery (mySession, track=track.name,
#                                                 table=table.name))
#   
# }




