 


################## QTL Data ################## 
# V- psychENCODE
# V- Fairfax
# V- MESA
# V- Cardiogenics
# - ROSMAP
# - ImmVar
# - STARNET
# V- GTEx (49 tissues)



psychENCODE.QTL_overlap <- function(FM_all=merge_finemapping_results(minimum_support = 0),
                                    consensus_only=T,
                                    local_files=T,
                                    force_new_subset=F){
  root <- ifelse(local_files, 
                 "./Data/QTL/psychENCODE",
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
    output_path <- file.path(dirname(ASSAY_files[[assay]]), assay, paste0("psychENCODE.",assay,".finemap.txt"))
    
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
      github_output <- file.path("./Data/QTL/psychENCODE",assay,paste0("psychENCODE.",assay,".finemap.txt"))
      dir.create(dirname(github_output), showWarnings = F, recursive = T)
      data.table::fwrite(QTL.sub, github_output, sep="\t", nThread = 4)
    }
    
    FM_all <- data.table:::merge.data.table(FM_all,
                                            QTL.sub,
                                            by = "SNP_id",
                                            all.x = T, allow.cartesian = T)
    # dat.merge[paste0("psychENCODE.",assay)] <- ifelse(subset(dat.merge, select=paste0("psychENCODE.",assay,".FDR")) <= 0.05, "Y","N")
  }
  return(FM_all)
}



# Fairfax: eQTL
Fairfax.QTL_overlap <- function(FM_all, conditions=c("CD14","IFN","LPS2","LPS24") ){
  FM_all <- dplyr::distinct(FM_all)
  for(condition in conditions){
    dataset <- paste0("Fairfax_2014_",condition)
    printer("Fairfax:: Processing:",dataset)
    server_path <- Directory_info("Fairfax_2014_CD14",variable = "fullSumStats") 
    output_path <- file.path("./Data/QTL/Fairfax_2014",condition,paste0(dataset,".finemap.txt"))
    dir.create(dirname(output_path),showWarnings = F, recursive = T)
    
    if(file.exists(output_path)){
      dat <- data.table::fread(output_path, nThread = 4)
    } else { 
      system(paste0("grep -E '", paste(unique(FM_all$SNP_id), collapse="|"),"' " ,server_path," > ",output_path))
      # Import data 
      dat <- data.table::fread(output_path, 
                               col.names = colnames(data.table::fread(server_path, nrow=0)), 
                               nThread = 4)#file_path
      ## Overwrite subset w/ fixed header
      data.table::fwrite(dat, output_path, sep = "\t", nThread = 4) 
    }  
    
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


MESA.QTL_overlap <- function(FM_all, force_new_subset=F){ 
  for(pop in c("AFA","CAU","HIS")){
    printer("MESA:: Extracting",pop,"data.")
    dataset <- paste0("MESA_",pop)
    output_path <- file.path("./Data/QTL/MESA",pop,paste0(dataset,".finemap.txt.gz"))
    
    if(file.exists(output_path) & force_new_subset==F){
      printer("MESA:: Pre-existing file detect. Importing...")
      # for some reason, saved files lose their header. Have to explicitly provide them instead...
      dat <- data.table::fread(output_path, nThread = 4, 
                               col.names = c("snps","gene","statistic","pvalue","FDR","beta",
                                             "chr","gene_name","start","end","gene_type","pos_snps","ref","alt"))
    } else{
      server_path <- Directory_info(paste0("MESA_",pop), "fullSumStats")
      output_path_txt <- gsub(".gz","",output_path)
      system(paste0("grep -E '", paste(unique(FM_all$Gene), collapse="|"),"' " ,server_path," > ",output_path_txt))
      # Import data 
      dat <- data.table::fread(output_path_txt, 
                               col.names = colnames(data.table::fread(server_path, nrow=0)), 
                               nThread = 4)#file_path
      dat <- subset(dat, snps %in% FM_all$SNP)
      ## Overwrite subset w/ fixed header (and smaller size)
      dir.create(dirname(output_path), recursive = T, showWarnings = F)
      data.table::fwrite(dat, output_path, sep = "\t", col.names = T)
      R.utils::gzip(output_path_txt, overwrite=T, remove=T)
    }  
      # Take only the top QTL per SNP location
      ## Rename columns
      dat.sub <- dplyr::select(dat, snps, beta, FDR, gene_name) %>% 
        dplyr::group_by(snps) %>%
        arrange(FDR, desc(beta)) %>% 
        dplyr::slice(1) %>% 
        `colnames<-`(c("SNP", paste0(dataset,".Effect"), paste0(dataset,".FDR"), paste0("MESA.",pop,".gene") )) %>%
        data.table::data.table()
   
    FM_all <- data.table:::merge.data.table(FM_all, dat.sub, 
                                            by="SNP", 
                                            all.x = T) 
  }
  return(FM_all)
}


Cardiogenics.QTL_overlap <- function(FM_all, force_new_subset=F, cis_only=T){
  for (celltype in c("macrophages","monocytes")){
    printer("Cardiogenics:: Processing",celltype,"data...")
    dataset <- paste0("Cardiogenics_",celltype)
    output_path <- file.path("./Data/QTL/Cardiogenics",celltype,paste0(dataset,".finemap.txt"))
    if(file.exists(output_path) & force_new_subset==F){
      print("Cardiogenics:: Pre-existing file detected. Importing...")
      dat.sub <- data.table::fread(output_path, nThread = 4)
    } else{ 
        server_path <- Directory_info(paste0("Cardiogenics_",celltype), "fullSumStats")
        dir.create(dirname(output_path), recursive = T, showWarnings = F)
        DAT <- data.table::fread(server_path,  nThread = 4)
        dat.sub <- subset(DAT, SNPID %in% unique(FM_all$SNP))
        data.table::fwrite(dat.sub, output_path, sep="\t", nThread = 4)
    }
    
    if(cis_only){dat.sub <- subset(dat.sub, relativePosition=="cis")}
    dat.sub <- dplyr::select(dat.sub, SNPID, beta, FDR, reporterID) %>% 
                dplyr::group_by(SNPID) %>%
                arrange(FDR, desc(beta)) %>% 
                dplyr::slice(1) %>% 
                `colnames<-`(c("SNP", 
                               paste0(dataset,".Effect"), 
                               paste0(dataset,".FDR"), 
                               paste0(dataset,".probe") )) %>%
                data.table::data.table()
    
    FM_all <- data.table:::merge.data.table(FM_all, dat.sub, 
                                            by="SNP", 
                                            all.x = T) 
  }
  return(FM_all)
}

GTEx.QTL_overlap <- function(FM_all, fuzzy_search=F){
  output_dir <- "./Data/QTL/GTEx"
  tissues_ref <- file.path(output_dir,"GTEx_tissues.txt") 
  server_dir <- Directory_info("GTEx","fullSumStats")
  # Identify which tissues we have eQTL data for
  if(file.exists(tissues_ref)){
    tissues_df <- data.table::fread(tissues_ref)
  } else{ 
     printer("GTEx:: Constructing reference file of available single-tissue eQTL files...")
     all_tissues <- lapply(list.files(server_dir), function(e){strsplit(e,"[.]")[[1]][1] }) %>% unlist() %>% unique()
     dir.create(output_dir, showWarnings = F, recursive = T)
     tissues_df <- data.frame(tissue=all_tissues, 
                              egene_file=list.files(server_dir, pattern = "*.egenes.*", full.names = T),
                              signif_pairs_file=list.files(server_dir, pattern = "*.signif_variant_gene_pairs.*", full.names = T) )
     data.table::fwrite(tissues_df, tissues_ref, sep="\t")
  }
  if(fuzzy_search!=F){
    tissues_df = tissues_df[grep(fuzzy_search, tissues_df$tissue)]
  }
  
  # Gather QTL data
  for(tiss in tissues_df$tissue){
    printer("GTEx:: Processsing eQTL for",tiss)
    dataset <- paste0("GTEx_",tiss)
    dat <- fread(subset(tissues_df,tissue==tiss)$egene_file %>% as.character(), nThread = 4)
    dat.sub <- subset(dat, (rs_id_dbSNP151_GRCh38p7 %in% unique(FM_all$SNP)))
    # Subset
    # **NOTE!!**: There's a A LOOOTT of different variables in the gtex files. 
    ## Need to go back and check if I'm using the right ones.
    dat.sub <- dplyr::select(dat.sub, rs_id_dbSNP151_GRCh38p7, beta_shape1, qval, gene_name) %>% 
      dplyr::group_by(rs_id_dbSNP151_GRCh38p7) %>%
      arrange(qval, desc(beta_shape1)) %>% 
      # dplyr::slice(1) %>% 
      `colnames<-`(c("SNP", 
                     paste0(dataset,".Effect"), 
                     paste0(dataset,".FDR"), 
                     paste0(dataset,".gene") )) %>%
      data.table::data.table() 
    FM_all <- data.table:::merge.data.table(FM_all, dat.sub, 
                                            by="SNP", 
                                            all.x = T, allow.cartesian = T)  
  }
 return(FM_all)
}


####----------- Gather QTL Overlap -----------####

merge_all_QTLs <- function(){
  # Gather all Fine-mapping results
  FM_all <- merge_finemapping_results(minimum_support = 0, 
                                      include_leadSNPs = T, 
                                      dataset = "./Data/GWAS/Nalls23andMe_2019")
  
  # psychENCODE eQTL, cQTL, isoQTL, HiC: DLPFC
  FM_merge <- psychENCODE.QTL_overlap(FM_all=FM_all, consensus_only = F,  local_files = T,  force_new_subset = F)
  # Fairfax eQTL: monocytes
  FM_merge <- Fairfax.QTL_overlap(FM_merge)
  # MESA eQTL: monocytes
  FM_merge <- MESA.QTL_overlap(FM_merge, force_new_subset = F) 
  # Cardiogenics eQTL: macrophages, monocytes
  FM_merge <- Cardiogenics.QTL_overlap(FM_merge, force_new_subset = F, cis_only = T)
  # GTEx eQTL: 49 different tissues
  FM_merge.final <- GTEx.QTL_overlap(FM_merge, fuzzy_search=F)
  
  ## Write file and compress
  QTL_merged_path <- file.path("./Data/GWAS/Nalls23andMe_2019/_genome_wide","Nalls23andMe_2019.QTL_overlaps.txt")
  data.table::fwrite(FM_merge.final, QTL_merged_path, sep="\t", nThread = 4)
  R.utils::gzip(QTL_merged_path, overwrite=T, remove=T)
  return(FM_merge)
}










