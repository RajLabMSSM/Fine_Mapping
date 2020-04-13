
# IMPACT #
# https://github.com/immunogenomics/IMPACT



IMPACT.get_annotation_key <- function(URL="https://github.com/immunogenomics/IMPACT/raw/master/IMPACT707/IMPACT_annotation_key.txt",
                                      save_path="./echolocatoR/annotations/IMPACT/IMPACT_annotation_key.txt.gz",
                                      force_new_download=F){
  if(file.exists(save_path) & force_new_download==F){
    print("+ IMPACT:: Importing local anotation key...")
    annot.key <- data.table::fread(save_path)
  } else {
    print("+ IMPACT:: Downloading annotation key from GitHub...")
    annot.key <- data.table::fread(URL)
    data.table::fwrite(annot.key, save_path, sep="\t")
    R.utils::gzip(save_path)
  }
  annot.key$Annot <- as.factor(paste0("Annot",annot.key$IMPACT))
  return(annot.key)
}


IMPACT.get_annotations <- function(baseURL="https://github.com/immunogenomics/IMPACT/raw/master/IMPACT707/Annotations",
                                   chrom=NULL, 
                                   subset_DT=NULL,
                                   nThread=4){ 
  # These are large files stored via GitHub's large file storage (lfs) 
  # Install LFS: https://git-lfs.github.com
  # Getting started with LFS: https://www.atlassian.com/git/tutorials/git-lfs
  # Download metadata / raw data for specific files with curl:  https://docs.gitlab.com/ee/api/repository_files.html#get-file-from-repository
  
  if(!is.null(subset_DT)){
    chrom <- subset_DT$CHR[1]
  }
  # https://github.com/immunogenomics/IMPACT.git
  # "curl --header 'https://git-lfs.github.com/spec/v1/projects/13083/repository/files/app%2Fmodels%2Fkey%2Erb/raw?ref=master'"
  
  
  URL <- file.path(baseURL, paste0("IMPACT707_EAS_chr",chrom,".annot.gz"))
  annot <- data.table::fread(URL, nThread = nThread)
  
  if(!is.null(subset_DT)){ 
    annot_merge <- data.table:::merge.data.table(data.table::data.table(subset_DT),
                                                 annot,
                                                 by.x = c("SNP","CHR","POS"), 
                                                 by.y = c("SNP","CHR","BP"),
                                                 all.x = T)
  } else {annot_merge <- annot}
  # Merge with metadata
  annot.key <- IMPACT.get_annotation_key()
  annot_cols <- grep("^Annot*",colnames(annot_merge), value = T)
  annot_melt <- data.table:::melt.data.table(annot_merge, measure.vars = annot_cols,
                                             variable.name = "Annot",
                                             value.name = "IMPACT_score") %>%
    data.table:::merge.data.table(annot.key,
                                  by="Annot", 
                                  allow.cartesian = T)
  return(annot_melt)
} 



# quick_finemap()
# annot_melt <- IMPACT.get_annotations(baseURL = "/Volumes/Steelix/IMPACT/IMPACT707/Annotations", subset_DT = subset_DT)

IMPACT.iterate_get_annotations <- function(merged_DT,
                                           IMPACT_score_thresh=.1){
  IMPACT_score_thresh = .1
  no_no_loci <- c("HLA-DRB5","MAPT","ATG14","SP1","LMNB1","ATP6V0A1",
                  # Tau region
                  "RETREG3","UBTF","FAM171A2","MAP3K14","CRHR1","MAPT-AS1","KANSL1","NSF","WNT3")
  merged_DT <- merge_finemapping_results(minimum_support=1,
                                         include_leadSNPs=T,
                                         dataset = "./Data/GWAS/Nalls23andMe_2019",
                                         xlsx_path=F,
                                         from_storage=T,
                                         consensus_thresh = 2,
                                         haploreg_annotation=F,
                                         biomart_annotation=F,
                                         verbose = F) %>%
    dplyr::rename(Locus=Gene) %>%
    subset(!(Locus %in% no_no_loci))
 
  #
  
  ANNOT_MELT <- lapply(unique(merged_DT$Locus), function(locus){
    message("+ IMPACT:: Gathering annotations for Locus = ",locus)
    subset_DT <- subset(merged_DT, Locus==locus)
    annot_melt <- IMPACT.get_annotations(baseURL = "/Volumes/Steelix/IMPACT/IMPACT707/Annotations", 
                                         subset_DT = subset_DT,
                                         nThread = 4)
    annot_melt <- subset(annot_melt, IMPACT_score>=IMPACT_score_thresh)
    printer("+ IMPACT::",nrow(annot_melt),"annotations found at IMPACT_score ≥", IMPACT_score_thresh)
    return(annot_melt)
  }) %>% data.table::rbindlist()
 return(ANNOT_MELT)
}


IMPACT.get_top_annotations <- function(){
  topIMPACT <- annot_melt %>%
    dplyr::select(-Accession, -File) %>%
    dplyr::group_by(SNP) %>% 
    dplyr::top_n(n=1, wt=IMPACT_score) %>%
    data.table::data.table() %>% 
    unique()
  topIMPACT
  
  subset(annot_melt, SNP=="rs7294619" )
}

IMPACT.compute_enrichment <- function(annot_melt, locus=NULL){  
  sum.IMPACT <- sum(annot_melt$IMPACT_score, na.rm=T)
  len.SNPs <- n_distinct(annot_melt$SNP, na.rm = T)
  # SNP.groups <- list("leadGWAS"=subset(annot_melt, leadSNP),
  #                    "UCS"=subset(annot_melt, Consensus_SNP),
  #                    "ABF_CS"=subset(annot_melt, ABF.Credible_Set>0),
  #                    "FINEMAP_CS"=subset(annot_melt, FINEMAP.Credible_Set>0),
  #                    "SUSIE_CS"=subset(annot_melt, SUSIE.Credible_Set>0),
  #                    "POLYFUN_CS"=subset(annot_melt, POLYFUN_SUSIE.Credible_Set>0),
  #                    "Consensus"=subset(annot_melt, Consensus_SNP))
  # enrich <- lapply(names(SNP.groups), function(snp.group){
  #   print(snp.group)
  #   e <- SNP.groups[[snp.group]] %>% 
  #     dplyr::group_by(TF, Tissue, Cell, CellDeriv) %>%
  #     dplyr::summarise(enrichment = (sum(IMPACT_score, na.rm = T) / sum.IMPACT) /
  #                        (n_distinct(SNP, na.rm = T) / len.SNPs) ) %>% 
  #     dplyr::arrange(-enrichment) %>% 
  #     data.table::data.table()
  #   e <- cbind(SNP.group=snp.group, e)
  #   return(e)
  # }) %>% data.table::rbindlist()
  annot_melt[is.na(annot_melt$IMPACT_score),"IMPACT_score"] <- 0

  SNP.groups <- list(
    "leadGWAS" = annot_melt %>% 
      dplyr::group_by(TF, Tissue, Cell, CellDeriv) %>%
      dplyr::summarise(enrichment = (sum(IMPACT_score[leadSNP], na.rm = T) / sum(IMPACT_score, na.rm = T)) /
                         (n_distinct(SNP[leadSNP], na.rm = T) / n_distinct(SNP, na.rm = T)) ),
    "UCS" = annot_melt %>% 
      dplyr::group_by(TF, Tissue, Cell, CellDeriv) %>%
      dplyr::summarise(enrichment = (sum(IMPACT_score[Support>0], na.rm = T) / sum(IMPACT_score, na.rm = T)) /
                         (n_distinct(SNP[Support>0], na.rm = T) / n_distinct(SNP, na.rm = T)) ),
    "ABF_CS" = annot_melt %>% 
      dplyr::group_by(TF, Tissue, Cell, CellDeriv) %>%
      dplyr::summarise(enrichment = (sum(IMPACT_score[ABF.Credible_Set>0], na.rm = T) / sum(IMPACT_score, na.rm = T)) /
                         (n_distinct(SNP[ABF.Credible_Set>0], na.rm = T) / n_distinct(SNP, na.rm = T)) ), 
    "FINEMAP_CS" = annot_melt %>% 
      dplyr::group_by(TF, Tissue, Cell, CellDeriv) %>%
      dplyr::summarise(enrichment = (sum(IMPACT_score[FINEMAP.Credible_Set>0], na.rm = T) / sum(IMPACT_score, na.rm = T)) /
                         (n_distinct(SNP[FINEMAP.Credible_Set>0], na.rm = T) / n_distinct(SNP, na.rm = T)) ), 
    "SUSIE_CS" = annot_melt %>% 
      dplyr::group_by(TF, Tissue, Cell, CellDeriv) %>%
      dplyr::summarise(enrichment = (sum(IMPACT_score[SUSIE.Credible_Set>0], na.rm = T) / sum(IMPACT_score, na.rm = T)) /
                         (n_distinct(SNP[SUSIE.Credible_Set>0], na.rm = T) / n_distinct(SNP, na.rm = T)) ),
    "POLYFUN_CS" = annot_melt %>% 
      dplyr::group_by(TF, Tissue, Cell, CellDeriv) %>%
      dplyr::summarise(enrichment = (sum(IMPACT_score[POLYFUN_SUSIE.Credible_Set>0], na.rm = T) / sum(IMPACT_score, na.rm = T)) /
                         (n_distinct(SNP[POLYFUN_SUSIE.Credible_Set>0], na.rm = T) / n_distinct(SNP, na.rm = T)) ), 
    "Consensus" = annot_melt %>% 
      dplyr::group_by(TF, Tissue, Cell, CellDeriv) %>%
      dplyr::summarise(enrichment = (sum(IMPACT_score[Consensus_SNP], na.rm = T) / sum(IMPACT_score, na.rm = T)) /
                         (n_distinct(SNP[Consensus_SNP], na.rm = T) / n_distinct(SNP, na.rm = T)) ) 
  )
  enrich <- data.table::rbindlist(SNP.groups, idcol = "SNP.group") %>% dplyr::arrange(-enrichment)
  enrich <- cbind(Locus=locus, enrich)
  enrich$TF <- factor(enrich$TF, ordered = T)
  enrich$SNP.group <- factor(enrich$SNP.group, levels=names(SNP.groups), ordered = T)
  return(enrich)
}



IMPACT.iterate_enrichment <- function(gwas_paths,
                                      annot_baseURL="../../data/IMPACT/IMPACT707/Annotations"){
  # gwas_paths <- list.files(path = "./Data/GWAS/Nalls23andMe_2019", pattern = "Multi-finemap_results.txt", recursive = T, full.names = T)
  # no_no_loci <- c("HLA-DRB5","MAPT","ATG14","SP1","LMNB1","ATP6V0A1")
  # gwas_paths <- gwas_paths[!basename(dirname(dirname(gwas_paths))) %in% no_no_loci]
 
  ENRICH <- lapply(gwas_paths, function(x){
    locus <- basename(dirname(dirname(x)))
    message(locus)
    enrich <- NULL
    try({
      subset_DT <- data.table::fread(x, nThread = 4)
      if(!"Locus" %in% colnames(subset_DT)){
        subset_DT <- cbind(Locus=locus, subset_DT) 
      }
      subset_DT <- find_consensus_SNPs(finemap_DT = subset_DT)
      annot_melt <- IMPACT.get_annotations(baseURL = annot_baseURL, 
                                           subset_DT = subset_DT, 
                                           nThread = 4) 
      enrich <- IMPACT.compute_enrichment(annot_melt = annot_melt,
                                          locus = locus)
    }) 
    return(enrich)
  }) %>% data.table::rbindlist(fill=T)
  # ENRICH
  return(ENRICH)
}

IMPACT.genome_wide_enrichment <- function(ENRICH){
  
}


IMPACT.plot_enrichment <- function(ENRICH){
  enrich_dat <- subset(ENRICH, !is.na(enrichment) & enrichment>=1) %>% 
    dplyr::group_by(Tissue, CellDeriv) %>%
    dplyr::top_n(n=1, wt=enrichment)
  # ep <- ggplot(enrich_dat, aes(x=TF, y=enrichment, fill=SNP.group)) +
  #   geom_col(position = "dodge") +
  #   geom_hline(yintercept = 1, linetype="dashed", alpha=.8) +
  #   facet_grid(facets = SNP.group ~ Tissue + CellDeriv, 
  #              switch = "y", space = "free_x",
  #              scales = "free_x") +
  #   theme_bw() + 
  #   theme(strip.text = element_text(angle=0), 
  #         axis.text.x = element_text(angle=45, hjust=1))
  ep <- ggplot(enrich_dat, aes(x=TF, y=SNP.group, fill=enrichment)) +
    geom_col(position = "dodge") +
    geom_hline(yintercept = 1, linetype="dashed", alpha=.8) +
    facet_grid(facets = . ~ Tissue + CellDeriv, 
               switch = "y", space = "free_x") +
    theme_bw() + 
    theme(
      # strip.text = element_text(angle=0), 
          axis.text.x = element_text(angle=45, hjust=1))
  print(ep)
}


IMPACT.plot_impact_score <- function(annot_melt, 
                                     save_path=F,
                                     show_plot=T){
  library(patchwork)
  library(ggridges)
  
  annot_melt$Mb <- annot_melt$POS/1000000
  # Get the SNP w/ the top impact score for each annotation
  annot_top <- annot_melt %>% 
    dplyr::group_by(Tissue, Cell, CellDeriv, TF) %>% 
    top_n(n=1, wt=IMPACT_score) %>% 
    data.table::data.table()
  # subset(annot_top,SNP %in% unique(subset(annot_melt, Consensus_SNP)$SNP))
  # annot_top
  
  
  # Reduce to smaller df to make plotting faster
  finemap_cols <- grep("*.PP$|*.Credible_Set$",colnames(annot_melt),value=T)
  annot_snp <- subset(annot_melt, select=c("SNP","CHR","POS","Mb","P","Consensus_SNP","leadSNP","Support",finemap_cols)) %>% unique()
  annot_snp <- dplyr::mutate(annot_snp, SNP.Group = ifelse(Consensus_SNP,"Consensus SNP",ifelse(leadSNP,"Lead GWAS SNP",ifelse(Support>0,"Credible Set SNP",NA))))
  labelSNPs <- construct_SNPs_labels(DT = annot_snp, lead=T, method=T, consensus=T)
  leader_SNP <- subset(labelSNPs, type=="Lead SNP")
  CS_set <- subset(labelSNPs, type=="Credible Set")
  # ggb <- GGBIO.plot(finemap_DT = annot_snp, LD_matrix = LD_matrix, 
  #            XGR_libnames = NULL, 
  #            save_plot=F,
  #            Nott_sn_epigenome=F) 
  # GWAS row
  gwas <- ggplot(annot_snp, aes(x=Mb, y=-log10(P), color=-log10(P))) +
    geom_point(size=1) +
    geom_point(data=leader_SNP, pch=18, fill=NA, size=2.5, color=leader_SNP$color) +
    # Green rings aronud Credible Set SNPs
    geom_point(data=CS_set, pch=21, fill=NA, size=2.5, color=CS_set$color, stroke=1, alpha=0.8) +
    ### Background color label
    ggrepel::geom_label_repel(data=labelSNPs,
                              aes(label=SNP),
                              color=NA,
                              # nudge_x = .5,
                              fill="black",
                              box.padding = .25,
                              label.padding = .25,
                              label.size=NA,
                              alpha=.6,
                              seed = 1,
                              size = 3,
                              min.segment.length = 1) +
    ### Foreground color label
    ggrepel::geom_label_repel(data=labelSNPs,
                              aes(label=SNP),
                              color=labelSNPs$color,
                              segment.alpha = .5,
                              # nudge_x = .5,
                              box.padding = .25,
                              label.padding = .25,
                              segment.size = 1,
                              fill = NA,
                              alpha=1,
                              seed = 1,
                              size = 3) + 
    ylim(c(0,max(-log10(annot_snp$P)))*1.1) +
    geom_vline(xintercept = unique(subset(labelSNPs, Consensus_SNP)$Mb), color="goldenrod2",
               alpha=1, size=.3, linetype='solid') +
    geom_vline(xintercept = unique(subset(labelSNPs, leadSNP)$Mb), color="red",
               alpha=1, size=.3, linetype='solid') +
    theme_bw()
  
  # finemaping rows
  finemap <- ggplot(annot_snp, aes(Mb, y=POLYFUN_SUSIE.PP, color=POLYFUN_SUSIE.PP)) +
    geom_point(size=1) +  
    scale_color_viridis_c( breaks=c(0,.5,1)) +
    geom_point(data=subset(annot_snp,POLYFUN_SUSIE.Credible_Set>0), pch=21, fill=NA, size=2.5, 
               color="green3", stroke=1, alpha=0.8) +
    ylim(c(0,1.1)) +
    geom_vline(xintercept = unique(subset(labelSNPs, Consensus_SNP)$Mb), color="goldenrod2",
               alpha=1, size=.3, linetype='solid') +
    geom_vline(xintercept = unique(subset(labelSNPs, leadSNP)$Mb), color="red",
               alpha=1, size=.3, linetype='solid') + 
    theme_bw()
  # print(finemap)
  
  # IMPACT rows
  impact <- ggplot(subset(annot_melt, IMPACT_score>0.5), 
                   aes(x=Mb, y=IMPACT_score, color=TF)) + 
    geom_point(show.legend = T) +
    # geom_col(position = "identity", show.legend = T) +
    facet_grid(facets = Tissue ~ ., switch = "y") +
    # ggridges::geom_ridgeline(aes(height = IMPACT_score), na.rm = T, size=.1, show.legend = F) +
    # ggridges::theme_ridges() +
    theme_bw() + 
    labs(y="IMPACT score per tissue") + 
    theme(strip.text.y = element_text(angle = 0), 
          panel.grid = element_blank(), 
          axis.title.y = element_text(vjust = .5))
  # impact 
  
  impact_plot <- gwas + finemap + impact + patchwork::plot_layout(ncol = 1, heights = c(.2,.2,1)) +
    geom_vline(xintercept = unique(subset(labelSNPs, Consensus_SNP)$Mb), color="goldenrod2",
               alpha=1, size=.3, linetype='solid') +
    geom_vline(xintercept = unique(subset(labelSNPs, leadSNP)$Mb), color="red",
               alpha=1, size=.3, linetype='solid')
  # impact_plot
  
  if(show_plot){print(impact_plot)}
  
  if(save_path!=F){
    # save_path="./Data/GWAS/Nalls23andMe_2019/LRRK2/IMPACT/LRRK2_IMPACT_plot.png"
    dir.create(dirname(save_path), showWarnings = F, recursive = T)
    printer("IMPACT:: Saving plot ==>",save_path)
    ggsave(save_path, plot=impact_plot, height=10, width=10)
  }
  return(impact_plot)
}




IMPACT.get_ldscores <- function(chrom=NULL, 
                                subset_DT=NULL,
                                nThread=4){ 
  warning("LDSCores do not include any SNPs with MAF<0.5%, as they are restricted to HapMap3 SNPs. \
This may affect subsequent analyss (e.g. fine-mapping).")
  if(!is.null(subset_DT)){
    chrom <- subset_DT$CHR[1]
  }
  baseURL <- "https://github.com/immunogenomics/IMPACT/raw/master/IMPACT707/LDscores"
  URL <- file.path(baseURL, paste0("IMPACT707_EAS_chr",chrom,".l2.ldscore.gz"))
  ldscore <- data.table::fread(URL, nThread = nThread)
  
  if(!is.null(subset_DT)){ 
    ldscore_merge <- data.table:::merge.data.table(data.table::data.table(subset_DT),
                                                   ldscore,
                                                   by.x = c("SNP","CHR","POS"), 
                                                   by.y = c("SNP","CHR","BP"),
                                                   all = F)  
    return(ldscore_merge)
  } else {
    return(ldscore)
  } 
}

