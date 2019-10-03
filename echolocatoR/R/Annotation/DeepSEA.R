## ~~~~~~~~~~~~~~~~~~        ~~~~~~~~~~~~~~~~~~ 
# ~~~~~~~~~~~~~~~~~~ DeepSEA ~~~~~~~~~~~~~~~~~~ 
## ~~~~~~~~~~~~~~~~~~        ~~~~~~~~~~~~~~~~~~ 

# Results generated using DeepSEA online server(http://deepsea.princeton.edu/help/)
## or the selene python module (https://github.com/FunctionLab/selene/blob/master/tutorials/analyzing_mutations_with_trained_models/analyzing_mutations_with_trained_models.ipynb).


DeepSEA.subset_kunkle <- function(){
  # Summary stats
  # kunk <- data.table::fread("./Data/GWAS/Kunkle_2019/Kunkle_etal_Stage1_results.txt.gz", nThread = 4) 
  kunk <- data.table::fread("/sc/orga/projects/ad-omics/data/AD_GWAS/Kunkle_2019/Kunkle_etal_Stage1_results.txt", 
                            nThread = 4) 
  kunk$Pvalue <- as.numeric(kunk$Pvalue) 
  kunk.sig <- subset(kunk, Pvalue<=0.05)
  # Locus coordinates 
  kunkle_loci <- readxl::read_excel("./Data/GWAS/Kunkle_2019/Kunkle2019_supplementary_tables.xlsx", 
                                    sheet = "Supplementary Table 8") %>%
    dplyr::select(Locus=`Lead SNV Gene`,SNP=`Top Associated SNV`, LD_block=`LD Block (GRCh37)`) %>% 
    data.table::data.table() %>% 
    tidyr::separate(LD_block, c("CHR","Start","End"), sep=":|-")
  
  kunkle.merge <- lapply(kunkle_loci$Locus, function(locus){
    printer("Merging data for locus:",locus)
    locus.sub <- subset(kunkle_loci, Locus==locus)
    kunk.sub <- subset(kunk.sig, (Chromosome==locus.sub$CHR) & (Position >= locus.sub$Start) & (Position <= locus.sub$End))
    if(nrow(kunk.sub)>0){
      kunk.sub <- cbind(Locus=locus.sub$Locus, kunk.sub)
      return(kunk.sub)
    } else {printer("No SNPs identified"); return(NULL)}
  }) %>% data.table::rbindlist(fill=T)
  
  # Tally before sig subset
  kunkle.merge %>% dplyr::group_by(Locus) %>% tally()
  kunkle.merge <- kunkle.merge %>% dplyr::mutate(FDR=p.adjust(Pvalue, method = "fdr", n = nrow(kunk)))
  # Tally after sig subset 
  subset(kunkle.merge, FDR<=0.05)  %>% dplyr::group_by(Locus) %>% tally() %>% data.frame()
  kunkle.merge.sig <- subset(kunkle.merge, FDR<=0.05) 
  
  data.table::fwrite(kunkle.merge.sig,"./Data/GWAS/Kunkle_2019/Kunkle_Stage1_sig.txt", sep="\t")
}




DeepSEA.vcf_subset <- function(){
  locus_list <- c("BIN1",
                 "TREM2",
                 # "NYAP1",
                 "PTK2B",
                 "SPI1",
                 "MS4A2",
                 "ABCA7",
                 "LRRK2")
  kunkle.sig <- data.table::fread("./Data/GWAS/Kunkle_2019/Kunkle_Stage1_sig.txt")
  kunkle_loci <- readxl::read_excel("./Data/GWAS/Kunkle_2019/Kunkle2019_supplementary_tables.xlsx", 
                     sheet = "Supplementary Table 8")
  if(locally){
    vcf_folder <- "~/Desktop/DeepSEA/ROSMAP_vcf"
  } else{
    vcf_folder <- "/sc/orga/projects/ad-omics/brian/ROSMAP"
  }
  
  for(locus in locus_list){  
    printer("DeepSEA:: Preparing locus",locus,"...")
    regions.file <- file.path(vcf_folder, paste0(locus,".regions.txt"))
    if(locus=="LRRK2"){
      finemap_DT <- data.table::fread("./Data/GWAS/Nalls23andMe_2019/LRRK2/Multi-finemap/Multi-finemap_results.txt")
      coords <- data.frame(CHROM=unique(finemap_DT$CHR),POS=sort(subset(finemap_DT, P<=5e-8)$POS) )
    }else { 
      kunkle.sub <- subset(kunkle.merge.sig, Locus==locus) 
      coords <- data.frame(CHROM=unique(kunkle.sub$Chromosome)[1], POS=kunkle.sub$Position) 
    }
    data.table::fwrite(coords, regions.file, sep="\t", col.names = F)
    
    chr <- coords$CHROM[1]
    vcf.gz <- file.path(vcf_folder, paste0("DEJ_11898_B01_GRM_WGS_2017-05-15_",chr,".recalibrated_variants.vcf.gz"))
    locus.vcf <- file.path(vcf_folder, paste0(locus,".vcf"))
    # Subset vcf with tabix ==> 
    ## ==> Force multiallelic sites to become biallelic ==> 
    ## filter snps with alt alleles longer than 100bp (DeepSEA's upper limit) ==>
    ## ==> Remove all but the first 5 cols  ==>
    ## ==> replace * with "'==>

    ## Important: need to include '-h' flag in Tabix to include the header so that bcftools can understand the file format.  
    cmd <- capture.output( 
      cat("tabix",vcf.gz,"-h -R",regions.file,
                  "| bcftools norm -m-", 
                  # "| bcftools filter -e '(STRLEN(ALT[0])>=90)'",
                    # "| bcftools query -H -f '%CHROM\\t%POS\\t%TYPE\\t%REF\\t%ALT\\n'", # Needs to be last bcftool command bc erases header
                  "| bcftools view -e '(STRLEN(REF)>100) | (STRLEN(ALT)>100)'",
                  "| cut -f1-5",
                  "| sed 's/\\*//g'",
                  "> ",locus.vcf)  
    )
    system(cmd) 
    cat("\n\n")
  } 
   
  results_links <- list(LRRK2="http://deepsea.princeton.edu/job/analysis/results/ea111803-8f88-4299-8a0e-872d8594b854",
                        ABCA7="http://deepsea.princeton.edu/job/analysis/results/5fa99587-1051-4a29-b961-5182a12ec81c",
                        BIN1="http://deepsea.princeton.edu/job/analysis/results/8753189d-75f7-47af-afd3-6ee7a21b6b66",
                        MS4A2="http://deepsea.princeton.edu/job/analysis/results/1732d2a7-0e0b-4c0b-aa71-cdc8e677aa79",
                        NYAP1="http://deepsea.princeton.edu/job/analysis/results/8682bd9e-6431-4b12-9c57-de9ab812d51d",
                        PTK2B="http://deepsea.princeton.edu/job/analysis/results/e2c72379-f163-4571-9e70-a19b828943b5",
                        SPI1="http://deepsea.princeton.edu/job/analysis/results/67c215c2-7474-46b9-b732-766f99391f35",
                        TREM1="http://deepsea.princeton.edu/job/analysis/results/19f228f9-ca64-4a62-bda2-f3e23d7b374f") 
}





DeepSEA.prepare_data <- function(results_path){
  
  DeepSEA.merge_data <- function(dat){
    dat.merged <- dat %>% dplyr::mutate(CHR=gsub("chr","",chr), POS=pos) %>% 
      merge.data.frame(finemap_DT,#[,c("CHR","POS","SNP")], 
                       by=c("CHR","POS")) %>% 
      data.table::data.table() %>% arrange(SNP)
    return(dat.merged)
  }
  # Fine-mapping results
  finemap_DT <- data.table::fread(file.path(results_path,"Multi-finemap/Multi-finemap_results.txt"))
  # SNP-wise 'Functional Significant Scores'
  funsig <- data.table::fread(file.path(results_path, "DeepSEA/infile.vcf.out.funsig")) %>% DeepSEA.merge_data() 
  # SNP-wise probabilites of being eQTL, GWAS, or HGMD (The Human Gene Mutation Database) hit
  snpclass <- data.table::fread(file.path(results_path, "DeepSEA/infile.vcf.out.snpclass")) %>% DeepSEA.merge_data()
  deepsea.dat <- cbind(snpclass, `Functional significance score`=funsig$`Functional significance score`)
  return(deepsea.dat)
}
# deepsea.dat <- DeepSEA.prepare_data(results_path = "./Data/GWAS/Nalls23andMe_2019/LRRK2")


DeepSEA.corrplot <- function(deepsea.dat, save_path=F){ 
  PP.cols <- grep(".Probability",colnames(deepsea.dat), value = T)
  cor.dat <- deepsea.dat[,c(PP.cols,"mean.PP","eQTL-probability","GWAS-probability","HGMD-probability","Functional significance score")]
  cor.df <- cor(cor.dat)
  cor.df[is.na(cor(cor.dat))] <- 0 
  if(save_path!=F){png(file = file.path(results_path,"DeepSEA","DeepSEA.corrplot.png") )}
  corrplot::corrplot(cor.df,
                     # title = "Correlations Between Fine-mapping PPs and DeepSEA Predictions",
                     # bg = "black",  
                     order = "hclust",
                     # addrect = 2,
                     type = "upper",
                     tl.col = "black", 
                     method = "ellipse",
                     addCoef.col = T,
                     number.cex = 1,
                     tl.srt = 45,
                     tl.cex = .5
                     # col = colorRampPalette(c("red","white","blue"))(200)
                     # col = RColorBrewer::brewer.pal(5,"RdBu")
                     ) 
  dev.off()
}

DeepSEA.plot_predictions <- function(deepsea.dat){ 
  deepsea.melt <- data.table::melt.data.table(data.table::data.table(deepsea.dat.all), 
                                             id.vars = "SNP",
                                             measure.vars = c("eQTL-probability",
                                                              "GWAS-probability",
                                                              "HGMD-probability",
                                                              "Functional significance score"),
                                             variable.name = "Prediction",
                                             value.name = "Probability")
  deepsea.melt$Prediction <- gsub("-| ","\n",deepsea.melt$Prediction)
  deepsea.melt$Prediction <- factor(deepsea.melt$Prediction, levels = unique(deepsea.melt$Prediction) )
  gp <- ggplot(data=deepsea.melt, aes(x=SNP, y=Probability, fill=Probability)) +
    geom_col() + 
    facet_grid(Prediction~.) + 
    theme_dark() + 
    ylim(c(0,1)) +   
    scale_fill_viridis_c() +
    geom_hline(yintercept = 0) + 
    geom_text(stat = 'identity',aes(label=SNP, hjust=-.2, vjust= -.5, srt=45),size=3 ) + 
    theme(axis.title.x=element_blank(), 
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          rect = element_rect(fill = "transparent"),
          panel.background = element_rect(fill = "transparent"),
          strip.text.y = element_text(angle = 0),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +  
    labs(title="DeepSEA Predictions")
  print(gp)
  if(save_path!=F){
    ggsave(filename = file.path(results_path,"DeepSEA","DeepSEA.predictions.png"),
           plot = gp, dpi = 600, height=10, width=10)
  } 
  return(gp)
}



