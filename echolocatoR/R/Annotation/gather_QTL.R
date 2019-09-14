 


################## QTL Data ################## 
# - psychENCODE
# - Fairfax
# - MESA
# - Cardiogenics
# - 


# Gather all Fine-mapping results
FM_all <- merge_finemapping_results(minimum_support = 0, 
                                    include_leadSNPs = T, 
                                    dataset = "./Data/GWAS/Nalls23andMe_2019")


# psychENCODE: "eQTL","cQTL","isoQTL"
FM_all <- psychENCODE.QTL_overlap(FM_all=FM_all, 
                                  consensus_only=F, 
                                  local_files = F, 
                                  force_new_subset = T)


# Fairfax: eQTL
Fairfax.QTL_overlap <- function(FM_all, conditions=c("CD14","IFN","LPS2","LPS24") ){
  FM_all <- dplyr::distinct(FM_all)
  for(condition in conditions){
    dataset <- paste0("Fairfax_2014_",condition)
    printer("+ Processing:",dataset)
    data_dirs <-list_Data_dirs(writeCSV = F)
    file_path <- subset(data_dirs, Dataset==dataset)$fullSumStats
    output_path <- file.path("./Data/eQTL/Fairfax_2014",condition,paste0(dataset,".finemap.txt"))
    dir.create(dirname(output_path),showWarnings = F, recursive = T)
    if(!file.exists(output_path)){
      system(paste0("grep -E '", paste(unique(FM_all$SNP_id), collapse="|"),"' " ,file_path," > ",output_path))
    }
    # Import data 
    dat <- data.table::fread(output_path, 
                             col.names = colnames(data.table::fread(file_path, nrow=0)), 
                             nThread = 4)#file_path
    ## Fix the header
    data.table::fwrite(dat, output_path)
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

FM_all <- Fairfax.QTL_overlap(FM_all)

