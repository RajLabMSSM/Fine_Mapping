 


################## QTL Data ################## 
# - psychENCODE
# - Fairfax
# - MESA
# - 


# Gather all Fine-mapping results
FM_all <- merge_finemapping_results(minimum_support = 0)


# psychENCODE: "eQTL","cQTL","isoQTL"
FM_all <- psychENCODE.QTL_overlap(FM_all=FM_all, consensus_only=F)

paste(FM_all$SNP_id, collapse = "|")

# Fairfax: eQTL
Fairfax.QTL_overlap <- function(){
  data_dirs <-list_Data_dirs(writeCSV = F)
  subset(data_dirs, Dataset==)
  base_url <- "schilb03@data4.hpc.mssm.edu:/sc/orga/projects/pd-omics"
 
  system(paste("scp", ""))
  
}


