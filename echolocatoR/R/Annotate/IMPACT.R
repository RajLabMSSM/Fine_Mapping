
# IMPACT #
# https://github.com/immunogenomics/IMPACT


IMPACT.get_annotations <- function(chrom=NULL, 
                                subset_DT=NULL,
                                nThread=4){ 
  if(!is.null(subset_DT)){
    chrom <- subset_DT$CHR[1]
  }
  baseURL <- "https://github.com/immunogenomics/IMPACT/raw/master/IMPACT707/Annotations"
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
   

IMPACT.get_annotation_key <- function(URL="https://github.com/immunogenomics/IMPACT/raw/master/IMPACT707/IMPACT_annotation_key.txt"){
  annot.key <- data.table::fread(URL)
  annot.key$Annot <- as.factor(paste0("Annot",annot.key$IMPACT))
  return(annot.key)
}

# quick_finemap()
# annot_melt <- IMPACT.get_annotations(subset_DT = subset_DT)
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
  impact <- ggplot(annot_melt, aes(x=Mb, y=Tissue, fill=Tissue, height = IMPACT_score)) + 
    # facet_grid(facets = Tissue ~.) +
    # geom_point(aes(color=IMPACT_score)) +
    ggridges::geom_ridgeline(na.rm = T, size=.1, show.legend = F) +
    # ggridges::theme_ridges() +
    theme_bw() + 
    labs(y="IMPACT score per tissue") + 
    theme(strip.text.y = element_text(angle = 0), 
          panel.grid = element_blank(), 
          axis.title.y = element_text(vjust = .5))
  
  impact_plot <- gwas + finemap + impact + patchwork::plot_layout(ncol = 1, heights = c(.2,.2,1)) +
    geom_vline(xintercept = unique(subset(labelSNPs, Consensus_SNP)$Mb), color="goldenrod2",
               alpha=1, size=.3, linetype='solid') +
    geom_vline(xintercept = unique(subset(labelSNPs, leadSNP)$Mb), color="red",
               alpha=1, size=.3, linetype='solid')
  
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

 