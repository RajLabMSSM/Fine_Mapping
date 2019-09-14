 


################## QTL Data ################## 
# - psychENCODE
# - Fairfax
# - MESA
# - Cardiogenics
# - 


# Gather all Fine-mapping results
FM_all <- merge_finemapping_results(minimum_support = 1)


# psychENCODE: "eQTL","cQTL","isoQTL"
FM_all <- psychENCODE.QTL_overlap(FM_all=FM_all, consensus_only=F, local_files = F)


# Fairfax: eQTL
Fairfax.QTL_overlap <- function(FM_all, conditions=c("CD14","IFN","LPS2","LPS24") ){
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
                             col.names = colnames(data.table::fread(file_path, nrow=0)))#file_path
    dat <- subset(dat, SNP %in% FM_all$SNP_id)
    dat.sub <- dplyr::select(dat, SNP, beta, FDR) %>% `colnames<-`(c("SNP_id",
                                                                         paste0(dataset,".Effect"),
                                                                         paste0(dataset,".FDR") ))
    FM_all <- data.table:::merge.data.table(FM_all, dat.sub, 
                                            by="SNP_id", 
                                            all.x = T, allow.cartesian = T) 
  }
 
 return(FM_all) 
}

FM_all <- Fairfax.QTL_overlap(FM_all)

