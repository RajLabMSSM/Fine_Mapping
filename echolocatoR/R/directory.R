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


list_Data_dirs <- function(writeCSV = "Results/directories_table.csv"){
  root <- "/sc/orga/projects"
  Data_dirs = list(
    # ++++++++ GWAS SUMMARY STATS ++++++++ # 
    # Nall et al. (2019) w/ 23andMe
    "Nalls23andMe_2019" = list(type="Parkinson's GWAS", 
                           topSS="Data/GWAS/Nalls23andMe_2019/Nalls2019_TableS2.xlsx",
                           # fullSS=file.path(root,"pd-omics/data/nallsEtAl2019/23andme/PD_all_post30APRIL2015_5.2_extended.txt")),
                           fullSS=file.path(root,"pd-omics/data/nallsEtAl2019/combined_meta/nallsEtAl2019_allSamples_allVariants.mod.txt"),
                           reference="https://www.biorxiv.org/content/10.1101/388165v3"),  
    
    ## IGAP 
    "Lambert_2013" = list(type="Alzheimer's GWAS", 
                          topSS="Data/GWAS/Lambert_2013/Lambert_2019_AD_GWAS.xlsx",
                          fullSS=file.path(root,"ad-omics/data/AD_GWAS/Lambert_2013/IGAP_stage_1.1000G.phase3.20130502.tsv"),
                          reference="https://www.nature.com/articles/ng.2802"),
    
    ## Marioni et al. (2018) 
    "Marioni_2018" = list(type="Alzheimer's GWAS", 
                          topSS="Data/GWAS/Marioni_2018/Marioni2018_supplementary_tables.xlsm",
                          fullSS=file.path(root,"ad-omics/data/AD_GWAS/Marioni_2018/Marioni2018.4_UK_Biobank_IGAP_17May2018.1000G.phase3.20130502.tsv"),
                          reference="https://www.nature.com/articles/s41398-018-0150-6"),
    
    ## Jansen et al. (2018) 
    "Posthuma_2018" = list(type="Alzheimer's GWAS", 
                           topSS="Data/GWAS/Posthuma_2018/Posthuma2018_suppTables.xlsx", 
                           fullSS=file.path(root,"ad-omics/data/AD_GWAS/Posthuma_2018/phase3.beta.se.hrc.txt"),
                           reference="https://www.nature.com/articles/s41588-018-0311-9"),
    
    ## Kunkle et al. (2018) Alzheimer's GWAS 
    "Kunkle_2019" = list(type="Alzheimer's GWAS",
                         topSS="Data/GWAS/Kunkle_2019/Kunkle2019_supplementary_tables.xlsx", 
                         fullSS=file.path(root,"ad-omics/data/AD_GWAS/Kunkle_2019/Kunkle_etal_Stage1_results.txt"),
                         # fullSS_stage2=file.path(root,"ad-omics/data/AD_GWAS/Kunkle_etal_Stage2_results.txt"),
                         reference="https://www.nature.com/articles/s41588-019-0358-2"),
    
    # ++++++++ eQTL SUMMARY STATS ++++++++ #
    ## MESA eQTLs: African Americans
    "MESA_AFA" = list(type="QTL",
                      topSS="Data/eQTL/MESA/AFA/AFA_eQTL_PTK2B.txt",
                      fullSS=file.path(root,"ad-omics/data/MESA/AFA_cis_eqtl_summary_statistics.txt"), 
                      reference="https://www.nhlbi.nih.gov/science/multi-ethnic-study-atherosclerosis-mesa"),
    ## MESA eQTLs: Caucasians
    "MESA_CAU" = list(type="QTL",
                      topSS="Data/eQTL/MESA/CAU/CAU_eQTL_PTK2B.txt",
                      fullSS=file.path(root,"ad-omics/data/MESA/CAU_cis_eqtl_summary_statistics.txt"),
                      reference="https://www.nhlbi.nih.gov/science/multi-ethnic-study-atherosclerosis-mesa"),
    ## MESA eQTLs: Hispanics
    "MESA_HIS" = list(type="QTL",
                      topSS="Data/eQTL/MESA/HIS/HIS_eQTL_PTK2B.txt",
                      fullSS=file.path(root,"ad-omics/data/MESA/HIS_cis_eqtl_summary_statistics.txt"),
                      reference="https://www.nhlbi.nih.gov/science/multi-ethnic-study-atherosclerosis-mesa"),
    
    ## Fairfax eQTLs: CD14
    "Fairfax_2014_CD14" = list(type="QTL",
                      topSS=NA,
                      fullSS=file.path(root,"ad-omics/data/fairfax/sumstats/cis.eqtls.fairfax.all.chr.CD14.47231.414.b.qced.f.txt"),
                      reference="https://science.sciencemag.org/content/343/6175/1246949"), 
    ## Fairfax eQTLs: IFN
    "Fairfax_2014_IFN" = list(type="QTL",
                               topSS=NA,
                               fullSS=file.path(root,"ad-omics/data/fairfax/sumstats/cis.eqtls.fairfax.all.chr.IFN.47231.367.b.qced.f.txt"),
                               reference="https://science.sciencemag.org/content/343/6175/1246949"),
     
  ## Fairfax eQTLs: IFN
  "Fairfax_2014_LPS2" = list(type="QTL",
                            topSS=NA,
                            fullSS=file.path(root,"ad-omics/data/fairfax/sumstats/cis.eqtls.fairfax.all.chr.LPS2.47231.261.b.qced.f.txt"),
                            reference="https://science.sciencemag.org/content/343/6175/1246949"), 
  ## Fairfax eQTLs: IFN
  "Fairfax_2014_LPS24" = list(type="QTL",
                             topSS=NA,
                             fullSS=file.path(root,"ad-omics/data/fairfax/sumstats/cis.eqtls.fairfax.all.chr.LPS24.47231.322.b.qced.f.txt"),
                             reference="https://science.sciencemag.org/content/343/6175/1246949"),   
  
  ## Cardiogenics: Macrophages
  "Cardiogenics_macrophages" = list(type="QTL",
                              topSS=NA,
                              fullSS=file.path(root,"ad-omics/data/cardiogenics/Cardiogenics/Macrophages.REPORT.fdr-0.5.tab"),
                              reference="https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003240"),
  ## Cardiogenics: Monocytes
  "Cardiogenics_monocytes" = list(type="QTL",
                                    topSS=NA,
                                    fullSS=file.path(root,"ad-omics/data/cardiogenics/Cardiogenics/Monocites.REPORT.fdr-0.5.tab"),
                                    reference="https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003240")  
          
  ) 
  Data_dirs_table <- Data_dirs_to_table(Data_dirs, writeCSV)
  return(Data_dirs_table)
}


Directory_info <- function(dataset_name, variable="topSumStats"){
  Data_dirs <- list_Data_dirs(writeCSV = F)
  directory = subset(Data_dirs, Dataset==dataset_name)[1,variable] %>% as.character() 
  return(directory)
}


get_dataset_name <- function(file_path){
  dataset_name <- tail(strsplit(dirname(file_path), "/")[[1]],n = 1)
  return(dataset_name)
}

delete_subset <- function (force_new_subset, subset_path){
  # Force new file to be made
  if(force_new_subset==T){
    printer("\n + Removing existing summary stats subset...\n")
    # dataset_name <- get_dataset_name(subset_path)
    # subset_path <- paste(dirname(subset_path),"/",gene,"_",superpopulation,"_",dataset_name,"_subset.txt",sep="") 
    suppressWarnings(file.remove(subset_path))
  }
}


get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}
# get_os() 


#### Make paths for results and subsets
make_results_path <- function(dataset_name, dataset_type, gene){ 
  results_path <- file.path("Data",dataset_type, dataset_name, gene)
  dir.create(results_path, showWarnings = F, recursive = T)
  return(results_path)
}

get_subset_path <- function(results_path, gene, subset_path="auto"){
  # Specify subset file name 
  if(subset_path=="auto"){  
    dataset_name <- basename(dirname(results_path))
    created_sub_path <- file.path(results_path, paste(gene,"_",dataset_name,"_subset.txt",sep="") )
    return(created_sub_path)
  } else{return(subset_path)}
}


