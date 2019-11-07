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
# 1. Should be able to continue running with missing SNPs. Should also report the number of missing SNPs without having to go into the file.
# 2. Verbose options to let you know what's going on at each step (preference).
# 3. If a user supplies both SNP and BP/CHR, use all of these to merge (to avoid duplicate col names).
# 4. Additional python dependency for step 2 of Approach 2: bitarray
# 5. Line 69 of munge_polyfun_sumstats.py: "SNP" col wasn't being renamed.
# 6. polyfun.py script reduced the number of SNPs from 7764212 to 34851, printing "[WARNING]  number of SNPs is smaller than 200k; this is almost always bad."
# 7. Received the following error while running polyfun.py: "ValueError: not all SNPs in the sumstats file are also in the annotations file"

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
  snp.path <- file.path(PT.path,"snps_to_finemap.txt.gz")
  printer("+ PolyFun:: Writing SNP file ==>",snp.path)
  data.table::fwrite(PF.dat, file = snp.path, 
                     nThread = 4, sep = " ")
  return(snp.path)
}


POLYFUN.initialize <- function(finemap_DT, 
                               results_path, 
                               dataset="Nalls23andMe_2019", 
                               locus=NULL){
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
                                           locus=NULL){
  PF.output.path <- file.path(results_path, "PolyFun")
  finemap_DT <- POLYFUN.initialize(finemap_DT=finemap_DT, 
                                results_path=results_path, locus = "LRRK2")
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
  merged_DT <- data.table::merge.data.table(finemap_DT, 
                               dplyr::select(priors, SNP, PolyFun.priors=SNPVAR) %>% 
                                 data.table::data.table(), by="SNP")
 
  
  LD_matrix <- readRDS(file.path(results_path,"plink/LD_matrix.RData"))
  LD_matrix <- LD_matrix[merged_DT$SNP, merged_DT$SNP]
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
  return(priors)
}

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
  ggsave(file.path(results_path,'Multi-finemap',"PolyFun.plot.png"), plot = gg,dpi = 400, height = 8, width = 7)
}






# %%%%%%%%%%%%%%%% PolyFun approach 2 %%%%%%%%%%%%%%%% 
## Computing prior causal probabilities via an L2-regularized extension of S-LDSC

# POLYFUN.prepare_boltlmm <- function(finemap_DT,
#                                     PF.output.path){
#   printer("POLYFUN:: Creating input summary stats file...")
#   ss.file <- file.path(PF.output.path, "boltlmm_sumstats.gz")
#   input.dat <-  dplyr::select(finemap_DT, SNP, 
#                               CHR, BP=POS, 
#                               ALLELE1=A1, ALLELE0=A2, 
#                               P_LINREG=P, BETA=Effect, SE=StdErr)
#   data.table::fwrite(input.dat, ss.file, nThread = 4, sep=" ")
#   return(ss.file)
# }


POLYFUN.L2_regularized_LDSC <- function(polyfun="./echolocatoR/tools/polyfun",
                                        PF.output.path,
                                        munged.path, 
                                        annotations.path=file.path(polyfun,"example_data/annotations."),
                                        weights.path=file.path(polyfun,"example_data/weights."),
                                        parametric=T){
  # Approach 2 = Parametric = no partitions
  # Approach 3 = Non-parametric = partitions
  partitions <- ifelse(parametric, "--no-partitions", "")
  out.path <- file.path(PF.output.path,"output")
  dir.create(out.path, showWarnings = F, recursive = T)
  cmd <- paste("python3",file.path(polyfun,"polyfun.py"),
                "--compute-h2-L2",
                partitions,
                "--output-prefix",file.path(out.path,"testrun"),
                "--sumstats",munged.path,
                "--ref-ld-chr",annotations.path,
                "--w-ld-chr",weights.path)
  print(cmd)
  system(cmd)
}


POLYFUN.calculate_new_priors <- function(polyfun="./echolocatoR/tools/polyfun", 
                                         finemap_DT=NULL, 
                                         sample.size=NULL,
                                         dataset="Nalls23andMe_2019",
                                         sample.size=1474097){
  results_path <- file.path(dirname(Directory_info(dataset_name = dataset, "fullSS.local")), "_genome_wide")
  PF.output.path <- file.path(results_path, "PolyFun")
 
  finemap_DT <- POLYFUN.initialize(finemap_DT=finemap_DT,
                                   results_path = results_path,
                                   dataset = dataset)
  
  # Requires space-delimited file with the following columns (munging can recognize several variations of these names):
  ## SNP CHR BP ....and....
  ## either a p-value, an effect size estimate and its standard error, a Z-score or a p-value 
  
  # Munge summary stats
  munged.path <- file.path(PF.output.path,"sumstats_munged.parquet")
  cmd <- paste("python3",file.path(polyfun,"munge_polyfun_sumstats.py"),
                "--sumstats",Directory_info(dataset_name = dataset, "fullSS.local"),
                "--n",sample.size, # Study sample size
                "--out",munged.path,
                "--min-info",0,
                "--min-maf", 0)
  print(cmd)
  system(cmd)
}

