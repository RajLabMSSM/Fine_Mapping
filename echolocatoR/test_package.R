# script to test installation of all required packages
# Jack Humphrey

source("echolocatoR/R/MAIN.R")
reticulate::use_condaenv("echolocatoR")

finemap_results <- list()
coloc_results <- list()

Data_dirs <- read.csv("Data/directories_table.csv")

quick_finemap(locus = "LRRK2")


quit()


lead.snp <- subset(finemap_DT, leadSNP)
merged <-data.table:::merge.data.table(finemap_DT,
                              data.table::data.table(r2=LD_matrix[lead.snp$SNP,]^2,
                                                     SNP=names(LD_matrix[lead.snp$SNP,])),
                              by="SNP")
subset(merged, r2>0.4)$POS %>% summary()
min_POS=min(subset(merged, r2>0.4)$POS)
max_POS=max(subset(merged, r2>0.4)$POS)
finemap_DT %>% arrange(P)

dataset_name <- "Nalls23andMe_2019"
top_SNPs <- Nalls_top_SNPs <- import_topSNPs(
  topSS_path = Directory_info(dataset_name, "topSS"),
  chrom_col = "CHR", position_col = "BP", snp_col="SNP",
  pval_col="P, all studies", effect_col="Beta, all studies", gene_col="Nearest Gene",
  caption= "Nalls et al. (2018) w/ 23andMe PD GWAS Summary Stats",
  group_by_locus = T,
  locus_col = "Locus Number",
  remove_variants = "rs34637584")
# Manually add new SNP of interest
top_SNPs <- data.table::rbindlist(list(top_SNPs,
                                      data.frame(CHR=14,POS=55360203, SNP="rs3783642",
                                                   P=2e-308, Effect=1, Gene="ATG14", Locus=NA ),
                                      data.frame(CHR=12,POS=53797236, SNP="rs34656641",
                                                   P=2e-308, Effect=1, Gene="SP1", Locus=NA ),
                                      data.frame(CHR=5,POS=126181862, SNP="rs959573",
                                                   P=2e-308, Effect=1, Gene="LMNB1", Locus=NA ),
                                      data.frame(CHR=17,POS=40609405, SNP="rs9897702",
                                                   P=2e-308, Effect=1, Gene="ATP6V0A1", Locus=NA ),
                                      data.frame(CHR=12,POS=40388109, SNP="rs140722239",
                                                   P=2.691e-37, Effect=0.3869, Gene="SLC2A13", Locus=NA )
                                      ), use.names = T)
