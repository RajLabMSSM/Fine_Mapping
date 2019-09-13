# ----------------------- #
# ----- ggbio plots ------#
# ----------------------- #

### ALL LOCI
 
ggbio.all_loci <- function(FM = merge_finemapping_results()){
  # Convert to GRange object
  DT <- FM %>% dplyr::mutate(SEQnames = paste0("chr",CHR)) 
  gr.snp <- biovizBase::transformDfToGr(DT, seqnames = "SEQnames", start = "POS", end = "POS")
  gr.snp_CHR <- biovizBase::transformDfToGr(FM, seqnames = "CHR", start = "POS", end = "POS")
  
  
  snp.track <- SNP_track(gr.snp, r2 = NULL, labels_subset = c("Lead SNP", "Consensus SNP"))
  snp.track 
}


### SINGLE LOCUS

GR.name_filter_convert <- function(GR.final, GR.names, min_hits=1){
  names(GR.final) <- GR.names
  grl <- GR.final[!as.logical(lapply(GR.final, is.null))]
  # Filter to those that had at least N hits 
  grl <- grl[as.logical(lapply(grl, function(g, min_hits.=min_hits){length(seqnames(g)) >= min_hits.}))]
  # Convert to GRangesList (important)
  grl <- GRangesList(grl)
  return(grl)
}


SNP_track <- function(gr.snp, 
                      method = "original",
                      labels_subset = c("Lead SNP", "Credible Set", "Consensus SNP"), 
                      r2=NULL){
  # Format data
  if(method!="original"){
    mcols(gr.snp)[,"PP"] <- mcols(gr.snp)[,paste0(method,".Probability")] 
  }
  dat <- as.data.frame(gr.snp)
  ### Label set
  labelSNPs <- construct_SNPs_labels(dat, lead=T, method=T, consensus=T) 
  leader_SNP <- subset(labelSNPs, type=="Lead SNP")
  CS_set <- subset(labelSNPs, type=="Credible Set")  
  label_tags <- subset(labelSNPs, (type %in% labels_subset))
  
  ## Make track
  if(method=="original"){
    cutoff <- -log10(5e-8)
    r2_multiply <- 100
    a1 <- plotGrandLinear(gr.snp, 
                   geom = "point", 
                   coord = "genome",
                   aes(y = -log10(P), x=POS, color=r2))
  } else {
    cutoff <- 0.5
    r2_multiply <- 5
    # Further filter label tags if plotting fine-mapping results
    label_tags <- label_tags[ (label_tags[paste0(method,".Credible_Set")]>0),] # IMPORTANT!: >0 (not TRUE)
    a1 <- plotGrandLinear(gr.snp, 
                   geom = "point", 
                   coord = "genome", 
                   aes(y = PP, x=POS, color=r2),
                   legend = FALSE)
  }
  
  a1 <- a1 + 
    scale_color_gradient(low="blue", high="red", limits = c(0,1)) + 
    geom_hline(yintercept=0,alpha=.5, linetype=1, size=.5) + 
    geom_hline(yintercept=cutoff, alpha=.5, linetype=2, size=.5, color="black")
  if(is.null(r2)){
    a1 <- a1 + stat_smooth(data=dat, aes(x=POS, y=r2*r2_multiply, fill=r2),
                color="turquoise",  se = F, formula = y ~ x, 
                method = 'loess', span=.1, size=.5, alpha=.5)
  }
  a1 <- a1 + 
    # Add diamond overtop leadSNP
    geom_point(data=leader_SNP, pch=18, fill=NA, size=4, color=leader_SNP$color) +
    # Green rings aronud Credible Set SNPs
    geom_point(data=CS_set, pch=21, fill=NA, size=4, color=CS_set$color, stroke=2, alpha=0.8) +
    ### Background color label
    geom_label_repel(data=label_tags, 
                     aes(label=SNP),
                     color=NA,
                     nudge_x = .5,
                     fill="black",
                     box.padding = .25,
                     label.padding = .25,
                     label.size=NA,
                     alpha=.8,
                     seed = 1, 
                     size = 3,
                     min.segment.length = 1) +
    ### Foreground color label
    geom_label_repel(data=label_tags, 
                     aes(label=SNP),
                     color=label_tags$color,
                     segment.alpha = .5,
                     nudge_x = .5,
                     box.padding = .25,
                     label.padding = .25,
                     segment.size = 1,
                     fill = NA,
                     alpha=1,
                     seed = 1,
                     size = 3,
                     min.segment.length = 1) + 
    theme_classic() 
  return(a1)
}


Roadmap_tabix <- function(results_path, chrom, min_pos, max_pos, eid, convert_to_GRanges=T){ 
  tbx_start = Sys.time()
  printer("++ Downloading Roadmap Chromatin Marks:",eid)
  fname <- paste0(eid,"_15_coreMarks_dense.bed.bgz")
  URL <- file.path("https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final",
                   fname) # _15_coreMarks_stateno.bed.gz
  cmd <- paste0("cd ",results_path,"&& tabix -p bed --begin 2 --end 3 ", URL," ",
                chrom,":",min_pos,"-",max_pos)
  out <- system(cmd, intern = T)
  dat <- data.table::fread(text = out, sep = "\t", select = 1:4, col.names = c("Chrom","Start","End","State"))
  dat$EID <- eid
  dat$File <- fname
  if(convert_to_GRanges){
    dat <- biovizBase::transformDfToGr(dat, seqnames = "Chrom", start = "Start", end="End")
  }
  tbx_end =  Sys.time()
  printer("BED subset downloaded in",round(tbx_end-tbx_start,3),"seconds")
  return(dat)
}

ROADMAP_track <- function(results_path, gr.snp, limit_files=NA){
  rm_start = Sys.time()
  RoadMap_ref <- GS_construct_reference() 
  if(!is.na(limit_files)){
    RoadMap_ref <- RoadMap_ref[1:limit_files,]
  }
  # Download via tabix (fast)
  counter <- 1
  gr.roadmap <- lapply(unique(RoadMap_ref$EID), function(eid, 
                                                         gr.snp.=gr.snp, 
                                                         results_path.=results_path){
    printer("+ Querying subset from Roadmap API:", eid," - ",counter,"/",length(unique(RoadMap_ref$EID)))
    counter <<- counter+1
    dat <- GRanges()
    try({
      dat <- Roadmap_tabix(results_path=results_path, 
                           chrom = unique(seqnames(gr.snp.)), 
                           min_pos = min(gr.snp.$POS), 
                           max_pos = max(gr.snp.$POS), 
                           eid=eid,
                           convert_to_GRanges=T)
    }) 
    if(length(seqnames(dat))>0){
      return(dat)
    } else{return(NULL)}
  }) 
  remove(counter)
  grl.roadmap <- GR.name_filter_convert(gr.roadmap, RoadMap_ref$Epigenome.name, min_hits=1)
  rm_end = Sys.time()
  printer("All downloads complete in",round(rm_end-rm_start,1),"minutes")
  return(grl.roadmap)
}



####### XGR track
XGR_track <- function(gr.snp, 
                      anno_data_path=file.path("echolocatoR/tools/Annotations", paste0("XGR_",lib.name,".rds")) , 
                      lib.name, 
                      save_xgr=T){ 
  library(GenomicRanges)
  if(file.exists(anno_data_path)){
    printer("")
    printer("+ Saved annotation file detected. Loading...")
    GR.annotations <- readRDS(anno_data_path) 
  } else {
    printer("")
    printer("+ XGR: Downloading...",lib.name)
    GR.annotations <- XGR::xRDataLoader(lib.name)
    if(save_xgr){
      saveRDS(GR.annotations, file = anno_data_path)
    } 
  } 
  GR.orig <- unlist(GR.annotations) 
  
  gr.xgr <- lapply(names(GR.orig), function(g, gr.snp. = gr.snp){
    printer("Finding overlap for:", g)  
    GR.overlap <- subsetByOverlaps(GR.orig[[g]], gr.snp.) 
    len <- length(seqnames(GR.overlap) ) 
    printer("   - Overlapping annotations = ",len)
    if(len>0){
      return(GR.overlap)
    } else{return(NULL)}
  }) 
  grl.xgr <- GR.name_filter_convert(gr.xgr, names(GR.orig), min_hits=1)
  return(grl.xgr)
}


save_annotations <- function(gr, anno_path, libName){
  dir.create(dirname(anno_path), showWarnings = F, recursive = T)
  saveRDS(gr, file.path(anno_path))
}



############ PLOT ALL TRACKS TOGETHER ############

ggbio_plot <- function(finemap_DT, 
                       LD_matrix,
                       gene,
                       results_path,
                       method_list=c("SUSIE","FINEMAP"),
                       XGR_libnames=c("ENCODE_TFBS_ClusteredV3_CellTypes",
                                      "ENCODE_DNaseI_ClusteredV3_CellTypes",
                                      "Broad_Histone")
                       ){
  # http://bioconductor.org/packages/release/bioc/vignettes/ggbio/inst/doc/ggbio.pdf
  library(ggbio)
  require(GenomicRanges)
  require(biovizBase)
 
  # finemap_DT <- data.table::fread("Data/GWAS/Nalls23andMe_2019/LRRK2/Multi-finemap/Multi-finemap_results.txt", sep="\t")
  # load("Data/GWAS/Nalls23andMe_2019/LRRK2/plink/LD_matrix.RData")
  ## To merge a GRangesList into a single GRanges object:
  # GR.merged <- unlist(grl)
  
  # Set up data 
  TRACKS_list <- NULL
  # Add LD into the DT
  LD_SNP <- subset(finemap_DT, leadSNP==T)$SNP
  LD_sub <- LD_with_leadSNP(LD_matrix, LD_SNP)
  DT <- data.table:::merge.data.table(finemap_DT, LD_sub, by = "SNP")
  # Convert to GRange object
  DT <- DT %>% dplyr::mutate(SEQnames = paste0("chr",CHR))
  gr.snp <- biovizBase::transformDfToGr(DT, seqnames = "SEQnames", start = "POS", end = "POS")  
  gr.snp_CHR <- biovizBase::transformDfToGr(DT, seqnames = "CHR", start = "POS", end = "POS")
  
  # Track 1: GWAS
  track.gwas <- SNP_track(gr.snp, 
                          method = "original", 
                          labels_subset = c("Lead SNP", "Consensus SNP"))
  TRACKS_list <- append(TRACKS_list, track.gwas)
  names(TRACKS_list)[1] <- "GWAS"
  # Tracks 1n: Fine-mapping
  for(m in method_list){
    track.finemapping <- SNP_track(gr.snp, method = m, 
                             labels_subset = c("Lead SNP", "Credible Set")) 
   
    TRACKS_list <- append(TRACKS_list, track.finemapping)
    names(TRACKS_list)[length(TRACKS_list)] <- m
  }
  
  # Track 2: Genes
  # library(EnsDb.Hsapiens.v75) 
  track.genes <- autoplot(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75, which = gr.snp_CHR, names.expr = "gene_name")
  TRACKS_list <- append(TRACKS_list, track.genes)
  names(TRACKS_list)[length(TRACKS_list)] <- "Gene Track"
  
  # Track 3: Annotation - XGR Annotations
  ## Download
  for(lib in XGR_libnames){
    anno_data_path <- file.path("echolocatoR/tools/Annotations", paste0("XGR_",lib,".rds")) 
    grl.xgr <- XGR_track(gr.snp, 
                          anno_data_path = anno_data_path, 
                          lib.name = lib, 
                          save_xgr=T)
    # grl.xgr <- check_saved_XGR(results_path, lib) 
    ## Make track
    track.xgr <- autoplot(grl.xgr, which = gr.snp, 
                           fill = "magenta",
                           color = NA, 
                           geom = "density",
                           alpha = .1) + theme_bw()
    TRACKS_list <- append(TRACKS_list, track.xgr)
    new_name <- paste(strsplit(lib,"_")[[1]], collapse="\n")
    names(TRACKS_list)[length(TRACKS_list)] <- new_name
  }
 
  # Track 4: Roadmap Chromatin Marks API
  ## Download 
  lib <- "Roadmap_ChromatinMarks_CellTypes"
  anno_path <- file.path(results_path, "Annotation",paste0("GRanges_",lib,".rds"))
  if(file.exists(anno_path)){
    printer("+ Saved annotation file detected. Loading...")
    grl.roadmap <- readRDS(anno_path)
  } else {
    grl.roadmap <- ROADMAP_track(results_path = results_path,
                                 gr.snp = gr.snp, 
                                 limit_files=NA)
    save_annotations(gr = grl.roadmap, anno_path = anno_path, libName = lib)
  } 
  ## Make track
  track.roadmap <- autoplot(grl.roadmap, which = gr.snp, 
                            fill="blue",
                            color=NA, 
                            geom = "density", 
                            alpha=.1) + theme_bw()
  TRACKS_list <- append(TRACKS_list, track.roadmap)
  names(TRACKS_list)[length(TRACKS_list)] <- "Roadmap\nChromatinMarks\nCellTypes"
  
  
  
  # Fuse all tracks  
  params_list <- list(title = paste0(gene," [",length(seqnames(gr.snp))," SNPs]"), 
                      label.text.cex = .7, 
                      label.bg.fill = "grey5",
                      label.text.color = "white",
                      label.text.angle = 0,
                      label.width = unit(5.5, "lines"),
                      xlim = c(min(gr.snp$POS), max(gr.snp$POS)))
  TRACKS_list <- append(TRACKS_list, params_list)
  trks <- suppressWarnings(do.call("tracks", TRACKS_list))
  ggsave(file.path(results_path,"Multi-finemap",paste0(gene,"_ggbio.png")))
  print(trks) 
  return(trks)
}




