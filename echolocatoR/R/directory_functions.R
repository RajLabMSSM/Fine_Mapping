# %%%%%%%%%%%%%%%%% #
####### DIRECTORY ####### 
# %%%%%%%%%%%%%%%%% #

Data_dirs_to_table <- function(Data_dirs, writeCSV=F){
  df <- data.frame()
  for(n in names(Data_dirs)){
    newRow <- data.frame(Data_dirs[n], stringsAsFactors = F)
    colnames(newRow) <- c("Type","topSumStats", "fullSumStats","Reference")
    df <- rbind(df, cbind(Dataset=n,newRow))
  }  
  createDT(df)
  if(writeCSV!=F){ 
    data.table::fwrite(df, writeCSV, quote = F, sep = ",", row.names = F) 
  }
  return(df)
}


list_Data_dirs <- function(){
  root <- "~/../../../sc/orga/projects"
  Data_dirs = list(
    # ++++++++ GWAS SUMMARY STATS ++++++++ # 
    # Nall et al. (2019) w/ 23andMe
    "Nalls_23andMe" = list(type="Parkinson's GWAS", 
                           topSS="Data/Parkinsons/Nalls_23andMe/Nalls2019_TableS2.xlsx",
                           # fullSS=file.path(root,"pd-omics/data/nallsEtAl2019/23andme/PD_all_post30APRIL2015_5.2_extended.txt")),
                           fullSS=file.path(root,"pd-omics/data/nallsEtAl2019/combined_meta/nallsEtAl2019_allSamples_allVariants.mod.txt"),
                           reference="https://www.biorxiv.org/content/10.1101/388165v3"),  
    
    ## IGAP 
    "Lambert_2013" = list(type="Alzheimer's GWAS", 
                          topSS="Data/Alzheimers/Lambert_2013/Lambert_2019_AD_GWAS.xlsx",
                          fullSS=file.path(root,"ad-omics/data/AD_GWAS/Lambert_2013/IGAP_stage_1.1000G.phase3.20130502.tsv"),
                          reference="https://www.nature.com/articles/ng.2802"),
    
    ## Marioni et al. (2018) 
    "Marioni_2018" = list(type="Alzheimer's GWAS", 
                          topSS="Data/Alzheimers/Marioni_2018/Marioni2018_supplementary_tables.xlsm",
                          fullSS=file.path(root,"ad-omics/data/AD_GWAS/Marioni_2018/Marioni2018.4_UK_Biobank_IGAP_17May2018.1000G.phase3.20130502.tsv"),
                          reference="https://www.nature.com/articles/s41398-018-0150-6"),
    
    ## Jansen et al. (2018) 
    "Posthuma_2018" = list(type="Alzheimer's GWAS", 
                           topSS="Data/Alzheimers/Posthuma_2018/Posthuma2018_suppTables.xlsx", 
                           fullSS=file.path(root,"ad-omics/data/AD_GWAS/Posthuma_2018/phase3.beta.se.hrc.txt"),
                           reference="https://www.nature.com/articles/s41588-018-0311-9"),
    
    ## Kunkle et al. (2018) Alzheimer's GWAS 
    "Kunkle_2019" = list(type="Alzheimer's GWAS",
                         topSS="Data/Alzheimers/Kunkle_2019/Kunkle2019_supplementary_tables.xlsx", 
                         fullSS=file.path(root,"ad-omics/data/AD_GWAS/Kunkle_2019/Kunkle_etal_Stage1_results.txt"),
                         # fullSS_stage2=file.path(root,"ad-omics/data/AD_GWAS/Kunkle_etal_Stage2_results.txt"),
                         reference="https://www.nature.com/articles/s41588-019-0358-2"),
    
    # ++++++++ eQTL SUMMARY STATS ++++++++ #
    ## MESA eQTLs: African Americans
    "MESA_AFA" = list(type="eQTL",
                      topSS="Data/eQTL/MESA/AFA_eQTL_PTK2B.txt",
                      fullSS=file.path(root,"ad-omics/data/mesa/AFA_cis_eqtl_summary_statistics.txt"), 
                      reference="https://www.nhlbi.nih.gov/science/multi-ethnic-study-atherosclerosis-mesa"),
    ## MESA eQTLs: Caucasians
    "MESA_CAU" = list(type="eQTL",
                      topSS="Data/eQTL/MESA/CAU_eQTL_PTK2B.txt",
                      fullSS=file.path(root,"ad-omics/data/mesa/CAU_cis_eqtl_summary_statistics.txt"),
                      reference="https://www.nhlbi.nih.gov/science/multi-ethnic-study-atherosclerosis-mesa"),
    ## MESA eQTLs: Hispanics
    "MESA_HIS" = list(type="eQTL",
                      topSS="Data/eQTL/MESA/HIS_eQTL_PTK2B.txt",
                      fullSS=file.path(root,"ad-omics/data/mesa/HIS_cis_eqtl_summary_statistics.txt"),
                      reference="https://www.nhlbi.nih.gov/science/multi-ethnic-study-atherosclerosis-mesa") 
  ) 
  Data_dirs_table <- Data_dirs_to_table(Data_dirs, writeCSV = "Results/directories_table.csv")
  return(Data_dirs)
}


Directory_info <- function(dataset_name, data="topSumStats"){
  directory = subset(Data_dirs, Dataset==dataset_name, select=c(data))[1,] %>% as.character()
  return(directory)
}
