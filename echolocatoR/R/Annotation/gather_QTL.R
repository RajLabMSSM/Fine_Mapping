 


################## QTL Data ################## 
# V- psychENCODE 
# V- Fairfax
#  - MESA
#  - Cardiogenics
#  - ImmVar
#  - STARNET
#  - GTEx


<<<<<<< HEAD

=======
>>>>>>> d832a3f42935b381dafc8800aa6af42762b78e31
psychENCODE.QTL_overlap <- function(FM_all=merge_finemapping_results(minimum_support = 0),
                                    consensus_only=T,
                                    local_files=T,
                                    force_new_subset=F){
<<<<<<< HEAD
  root <- ifelse(local_files, 
                 "./Data/QTL/psychENCODE/",
=======
  root <- ifelse(local_files, "./echolocatoR/tools/Annotations/psychENCODE",
>>>>>>> d832a3f42935b381dafc8800aa6af42762b78e31
                 "/sc/orga/projects/ad-omics/data/psychENCODE")
  ASSAY_files <- file.path(root,
                           c("DER-08a_hg19_eQTL.significant.txt.gz",
                             "DER-09_hg19_cQTL.significant.txt.gz",
                             "DER-10a_hg19_isoQTL.significant.txt.gz",
                             "INT-16_HiC_EP_linkages_cross_assembly.csv.gz"))
  ASSAY_files = setNames(ASSAY_files, nm = c("eQTL","cQTL","isoQTL","HiC"))
  
  if(consensus_only){
    FM_all <- subset(FM_all, Consensus_SNP==T)
  }
  # Add SNP_id column
  FM_all <- FM_all %>% 
    dplyr::mutate(SNP_id=paste0(CHR,":",POS)) %>% 
    data.table::data.table()
  
  for(assay in  c("eQTL","cQTL","isoQTL")){
    printer("psychENCODE:: Checking for overlap with",assay)
    output_path <- file.path(dirname(ASSAY_files[[assay]]),paste0("psychENCODE.",assay,".finemap.txt"))
    
    if(file.exists(output_path) & force_new_subset==F){
      printer("psychENCODE:: Importing pre-existing file.")
      QTL.sub <- data.table::fread(output_path)
    } else {
      QTL <- data.table::fread(ASSAY_files[[assay]], nThread = 4)
      # Subset QTL data
      QTL.sub <- subset(QTL, SNP_id %in% FM_all$SNP_id)
      # QTLs have can multiple probes per SNP location. 
      ## Pick only the best one (lowest FDR and highest Effect) keep each row as a unique genomic position:
      QTL.sub <- subset(QTL.sub, select=c("SNP_id","regression_slope","FDR")) %>% 
        dplyr::group_by(SNP_id) %>%
        arrange(FDR, desc(regression_slope)) %>% 
        dplyr::slice(1) %>% 
        data.table::data.table()
      # Select and rename columns
      QTL.sub <-  QTL.sub %>%  
        `colnames<-`(c("SNP_id",
                       paste0("psychENCODE.",assay,".Effect"),
                       paste0("psychENCODE.",assay,".FDR"))) %>% unique() %>% 
        data.table::data.table()
<<<<<<< HEAD
      github_output <- file.path("./Data/QTL/psychENCODE",assay,paste0("psychENCODE.",assay,".finemap.txt"))
      dir.create(dirname(github_output), showWarnings = F, recursive = T)
      data.table::fwrite(QTL.sub, github_output, sep="\t", nThread = 4)
=======
      github_output <- file.path("./Data/eQTL/psychENCODE",assay,paste0("psychENCODE.",assay,".finemap.txt"))
      dir.create(dirname(github_output), showWarnings = F, recursive = T)
      data.table::fwrite(QTL.sub, github_output)
>>>>>>> d832a3f42935b381dafc8800aa6af42762b78e31
    }
    
    FM_all <- data.table:::merge.data.table(FM_all,
                                            QTL.sub,
                                            by = "SNP_id",
                                            all.x = T, allow.cartesian = T)
    # dat.merge[paste0("psychENCODE.",assay)] <- ifelse(subset(dat.merge, select=paste0("psychENCODE.",assay,".FDR")) <= 0.05, "Y","N")
  }
  return(FM_all)
}
<<<<<<< HEAD

=======
>>>>>>> d832a3f42935b381dafc8800aa6af42762b78e31


# Fairfax: eQTL
Fairfax.QTL_overlap <- function(FM_all, conditions=c("CD14","IFN","LPS2","LPS24") ){
  FM_all <- dplyr::distinct(FM_all)
  for(condition in conditions){
    dataset <- paste0("Fairfax_2014_",condition)
    printer("+ Processing:",dataset)
    data_dirs <-list_Data_dirs(writeCSV = F)
    file_path <- subset(data_dirs, Dataset==dataset)$fullSumStats
    output_path <- file.path("./Data/QTL/Fairfax_2014",condition,paste0(dataset,".finemap.txt"))
    dir.create(dirname(output_path),showWarnings = F, recursive = T)
    if(!file.exists(output_path)){
      system(paste0("grep -E '", paste(unique(FM_all$SNP_id), collapse="|"),"' " ,file_path," > ",output_path))
    }
    # Import data 
    dat <- data.table::fread(output_path, 
                             col.names = colnames(data.table::fread(file_path, nrow=0)), 
                             nThread = 4)#file_path
<<<<<<< HEAD
    ## Fix the header
    data.table::fwrite(dat, output_path,sep = "\t", nThread = 4)
=======
    ## Write after fixing the header
    data.table::fwrite(dat, output_path)
>>>>>>> d832a3f42935b381dafc8800aa6af42762b78e31
    dat <- subset(dat, SNP %in% FM_all$SNP_id)
    # Take only the top QTL per SNP location
    ## Rename columns
    dat.sub <- dplyr::select(dat, SNP, beta, FDR) %>% 
                      dplyr::group_by(SNP) %>%
                      arrange(FDR, desc(beta)) %>% 
                      dplyr::slice(1) %>% 
                      `colnames<-`(c("SNP_id", paste0(dataset,".Effect"), paste0(dataset,".FDR") )) %>%
                      data.table::data.table()
    FM_all <- data.table:::merge.data.table(FM_all, dat.sub, 
                                            by="SNP_id", 
                                            all.x = T, allow.cartesian = T) 
  } 
 return(FM_all) 
}



<<<<<<< HEAD
####----------- Gather QTL Overlap -----------####
=======
#~~~~~~~~~~~~~~~~~~~~~~~~~ Gather! ~~~~~~~~~~~~~~~~~~~~~~~~~# 
>>>>>>> d832a3f42935b381dafc8800aa6af42762b78e31

# Gather all Fine-mapping results
FM_all <- merge_finemapping_results(minimum_support = 0, 
                                    include_leadSNPs = T, 
<<<<<<< HEAD
                                    dataset = "./Data/GWAS/Nalls23andMe_2019")
# psychENCODE: "eQTL","cQTL","isoQTL"
FM_all <- psychENCODE.QTL_overlap(FM_all=FM_all, 
                                  consensus_only=F, 
                                  local_files = T, 
                                  force_new_subset = T)
# Fairfax eQTL
FM_all <- Fairfax.QTL_overlap(FM_all, )


=======
                                    dataset = "./Data/GWAS/Nalls23andMe_2019",
                                    xlsx_path = F)

################################
# psychENCODE eQTL, cQTL, isoQTL
FM_all <- psychENCODE.QTL_overlap(FM_all=FM_all, 
                                  consensus_only=F, 
                                  local_files = F, 
                                  force_new_subset = T)
#############################################
# Fairfax eQTL: Naive?, CD14, IFN, LPS2, LPS24
FM_all <- Fairfax.QTL_overlap(FM_all)
>>>>>>> d832a3f42935b381dafc8800aa6af42762b78e31



