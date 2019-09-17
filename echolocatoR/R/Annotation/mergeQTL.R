 


################## QTL Data ################## 
# V- psychENCODE
# V- Fairfax
# V- MESA
# V- Cardiogenics
# - ROSMAP
# - ImmVar
# - STARNET
# V- GTEx (49 tissues)



mergeQTL.psychENCODE <- function(FM_all=merge_finemapping_results(minimum_support = 0),
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
                       paste0("psychENCODE_",assay,".Effect"),
                       paste0("psychENCODE_",assay,".FDR"))) %>% unique() %>% 
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
mergeQTL.Fairfax <- function(FM_all, conditions=c("CD14","IFN","LPS2","LPS24") ){
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


mergeQTL.MESA <- function(FM_all, force_new_subset=F){ 
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


mergeQTL.Cardiogenics <- function(FM_all, force_new_subset=F, cis_only=T){
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

mergeQTL.GTEx <- function(FM_all, fuzzy_search=F){
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

mergeQTL.merge_all <- function(){
  # Gather all Fine-mapping results
  FM_all <- merge_finemapping_results(minimum_support = 0, 
                                      include_leadSNPs = T, 
                                      dataset = "./Data/GWAS/Nalls23andMe_2019")
  
  # psychENCODE eQTL, cQTL, isoQTL, HiC: DLPFC
  FM_merge <- mergeQTL.psychENCODE(FM_all=FM_all, consensus_only = F,  local_files = T,  force_new_subset = T)
  # Fairfax eQTL: monocytes
  FM_merge <- mergeQTL.Fairfax(FM_merge)
  # MESA eQTL: monocytes
  FM_merge <- mergeQTL.MESA(FM_merge, force_new_subset = F) 
  # Cardiogenics eQTL: macrophages, monocytes
  FM_merge <- mergeQTL.Cardiogenics(FM_merge, force_new_subset = F, cis_only = T)
  # GTEx eQTL: 49 different tissues
  FM_merge.final <- mergeQTL.GTEx(FM_merge, fuzzy_search=F)
  
  ## Write file and compress
  QTL_merged_path <- file.path("./Data/GWAS/Nalls23andMe_2019/_genome_wide",
                               "Nalls23andMe_2019.QTL_overlaps.txt")
  data.table::fwrite(FM_merge.final, QTL_merged_path, sep="\t", nThread = 4)
  R.utils::gzip(QTL_merged_path, overwrite=T, remove=T)
  return(FM_merge)
}


mergeQTL.melt_FDR <- function(FM_merge){ 
  dim(FM_merge) 
  Effect.cols <- grep(".Effect", colnames(FM_merge), value = T)
  FDR.cols <- grep(".FDR", colnames(FM_merge), value = T)
  id.vars <- c("Gene","SNP","CHR","POS","P","Consensus_SNP","Support","leadSNP")
  FM_sub <- subset(FM_merge, select=c(id.vars, FDR.cols))
  colnames(FM_sub)[colnames(FM_sub) %in% FDR.cols] <- gsub(".FDR","",FDR.cols)
  FM_melt <- data.table::melt.data.table(FM_sub, id.vars = id.vars, 
                                         variable.name = "QTL.Source",
                                         value.name = "FDR") 
  FM_melt[is.infinite(FM_melt$FDR),"FDR"] <- NA 
  FM_melt[FM_melt$FDR %in% c(Inf,-Inf),"FDR"] <- NA 
  FM_melt[FM_melt$FDR==0,"FDR"] <- .Machine$double.xmin 
  return(FM_melt)
}

mergeQTL.count_overlap <- function(){
  library(dplyr)
  FM_merge <- data.table::fread(file.path("./Data/GWAS/Nalls23andMe_2019/_genome_wide",
                                          "Nalls23andMe_2019.QTL_overlaps.txt.gz"), nThread = 4)
  FM_melt <- mergeQTL.melt_FDR(FM_merge) 
  
  mergedQTL.get_count <- function(SNP.subset){
    count.df <- subset(SNP.subset, (FDR<0.05), drop=F) %>% 
                      dplyr::group_by(QTL.Source, .drop=F) %>% 
                      tally(sort = F) %>% 
                      arrange(QTL.Source)
    count.df$all.SNPs <- nrow(SNP.subset)
    count.df$Proportion <- count.df$n / count.df$all.SNPs
    return(count.df)
  }
  
  consensus <- mergedQTL.get_count(subset(FM_melt, (Consensus_SNP==T)))
  cred.set <- mergedQTL.get_count(subset(FM_melt, (Support>0)))
  leadSNP <- mergedQTL.get_count(subset(FM_melt, (leadSNP==T)))
   
  
  merged.data <- data.table::data.table(QTL.Source = consensus$QTL.Source,
                                        Consensus = consensus$Proportion,
                                        Credible.Set = cred.set$Proportion,
                                        GWAS.lead.SNP = leadSNP$Proportion) %>% 
                data.table::melt.data.table(id.vars = "QTL.Source", 
                                            variable.name = "SNP.Group", 
                                            value.name = "Proportion.Overlap")
  library(ggplot2)
  ggplot(merged.data, aes(x=QTL.Source, y=Proportion.Overlap, fill=SNP.Group)) + 
    geom_col(position = position_dodge(preserve = "single")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) +
    labs(title="Proportion of Overlapping QTL")
   
    
}





######### PLOTS #############


mergeQTL.stacked_manhattan_plot <- function(GENE_df, SNP.Group){  
  op <- ggplot(GENE_df, aes(x=POS, y=-log10(FDR), color=-log10(FDR))) + 
    geom_point() + 
    facet_grid(QTL.Source~., drop = F) + 
    theme(strip.text.y = element_text(angle=0), 
          strip.background = element_rect(fill="slategray1"),
          strip.text = element_text(colour = 'black'),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) +  
    labs(subtitle=paste0(SNP.Group,": Significant QTL Overlap")) + 
    geom_hline(yintercept = -log10(0.05), linetype="dashed", alpha=.9, color="purple", size=.2) + 
    ylim(c(0, -log10(min(GENE_df$FDR))*1.2  )) + 
    geom_point(data = subset(GENE_df, FDR>=0.05), aes(x=POS, y=-log10(FDR)), color="gray")
  return(op)
}

mergeQTL.QTL_distributions_plot <- function(){
  FM_merge <- data.table::fread(file.path("./Data/GWAS/Nalls23andMe_2019/_genome_wide",
                                          "Nalls23andMe_2019.QTL_overlaps.txt.gz"), 
                                nThread = 4) 
  FM_melt <- melt_FDR.QTL_overlap(FM_merge)
  
  
  plot.QTL_distributions.subset <- function(SNP.Group="Consensus_SNP", cardiogenics=F){
    FM_melt <-  FM_melt %>% arrange(QTL.Source, FDR)
    if(cardiogenics==F){
      FM_melt <- subset(FM_melt, !(QTL.Source %in% c("Cardiogenics_macrophages","Cardiogenics_monocytes")))
    }
    
    if(SNP.Group=="Consensus_SNP"){
      GENE_df <- subset(FM_melt, (Consensus_SNP==T) & (FDR<0.05), drop=F) 
    }
    if(SNP.Group=="Credible_Set"){
      GENE_df <- subset(FM_melt, (Support>0) & (FDR<0.05), drop=F)  
    }
    if(SNP.Group=="leadSNP"){
      GENE_df <- subset(FM_melt, (leadSNP==T) & (FDR<0.05), drop=F)  
    } 
    GENE_df$SNP <- factor(GENE_df$SNP, levels = unique(GENE_df$SNP), ordered = T)
    
    # op <- plot.QTL_overlap(GENE_df, SNP.Group = SNP.Group)
    # print(op)
    
    gp <- ggplot(GENE_df , aes(x=QTL.Source, y=-log10(FDR), fill=QTL.Source, group=SNP)) + 
      geom_col(position = position_dodge(preserve = "single"), show.legend = F, aes(group=SNP), alpha=.7) + 
      geom_hline(yintercept = -log10(0.05), linetype="dashed", alpha=.9, color="purple", size=.2) + 
      theme_classic() + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5)) + 
      labs(title=NULL, subtitle=paste0(SNP.Group, " : Significant QTL Overlap"), x=NULL) + 
      scale_fill_discrete(drop=F) + 
      scale_x_discrete(drop=F)  
    # ylim(c(0,-log10( min(FM_melt$FDR, na.rm = T)) ))
    # ylim(c(0,150))
    # print(gp)
    return(gp) 
  } 
  consensus <- plot.QTL_distributions.subset(SNP.Group = "Consensus_SNP")
  cred.set <- plot.QTL_distributions.subset(SNP.Group = "Credible_Set")
  lead.snp <- plot.QTL_distributions.subset(SNP.Group = "lead_SNP")
  
  cowplot::plot_grid(consensus, cred.set, lead.snp, ncol = 1, align = "h")
  
  
}



cell <- readxl::read_excel("/sc/orga/projects/ad-omics/data/Brain_xQTL_Serve/CellSpecifictyeQTLs.xlsx", skip = 2)
data.table::fwrite(cell, "/sc/orga/projects/ad-omics/data/Brain_xQTL_Serve/cell-specificity-eQTLs.tsv", sep="\t")




