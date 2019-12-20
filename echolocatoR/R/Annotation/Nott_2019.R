# 
# ^^^^^^^^^^^^^^ Nott et al. (2019) # ^^^^^^^^^^^^^^
# ^^^^^^^^^^^^^^ single-nucleus human brain epigenomic dataset # ^^^^^^^^^^^^^^ 
# Nott, Alexi, Inge R. Holtman, Nicole G. Coufal, Johannes C.M. Schlachetzki, Miao Yu, Rong Hu, Claudia Z. Han, et al. “Cell Type-Specific Enhancer-Promoter Connectivity Maps in the Human Brain and Disease Risk Association.” Science 0793, no. November (2019): 778183. https://doi.org/10.1101/778183.
#  ^^^^^^^^^^^^^^ ^^^^^^^^^^^^^^ ^^^^^^^^^^^^^^ ^^^^^^^^^^^^^^ ^^^^^^^^^^^^^^


GGBIO.ucsc_tracks <- function(finemap_DT){ 
  # GLASS DATA: UCSC GB
  # https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr2:127770344-127983251&hgsid=778249165_ySowqECRKNxURRn6bafH0yewAiuf
   
  # finemap_DT <- quick_finemap("TRIM40")
  
  # UCSC Tracks   
  import.bw.filt <- function(bw.file, gr.dat){  
    bw.dat <- rtracklayer::BigWigSelection(ranges = gr.dat, 
                                           colnames = "score")  
    # import.bw:: https://rdrr.io/bioc/Gviz/src/R/Gviz.R#sym-.import.bw
    bw.filt <- rtracklayer::import.bw(con = bw.file, selection= gr.dat)
    # bw.filt <- Gviz:::.import.bw(bw.file, selection = gr.dat)
    # gr <- bw.filt@range
    # gr.filt <- gr[seqnames(gr) == chr & start(gr) > from & end(gr) < to] 
    plot(x = start(bw.filt), y=bw.filt$score)
    return(bw.filt)
  }
  # Import BigWig annotation files
  bigWigFiles <- readxl::read_excel("./echolocatoR/annotations/Glass_lab/Glass.snEpigenomics.xlsx") 
  bigWigFiles <- subset(bigWigFiles, marker!="-" &  cell_type!="peripheral microglia")
  
  # Convert finemap data to granges
  dat <- finemap_DT
  dat$seqnames <- paste0("chr",dat$CHR) 
  dat$start.end <- dat$POS
  gr.dat <- GenomicRanges::makeGRangesFromDataFrame(df = dat,
                                                    seqnames.field = "seqnames",
                                                    start.field = "start.end", 
                                                    end.field = "start.end", 
                                                    keep.extra.columns = T)   
  bw.grlist <- lapply(1:nrow(bigWigFiles), function(i){ 
    bw.file <- bigWigFiles$data_link[i]
    bw.name <- gsub("_pooled|pooled_","",bigWigFiles$name[i])
    printer("GVIZ:: Importing...",bw.name)
    bw.filt <- import.bw.filt(bw.file=bw.file, 
                              gr.dat = gr.dat)
    colnames(mcols(bw.filt))[1] <- bw.name
    # bw.filt$expt_name <- bw.name  
    # bw.filt$cell_type <-strsplit(bw.name, "_")[[1]][[1]]
    # bw.filt$assay <- strsplit(bw.name, "_")[[1]][[2]]
    return(bw.filt)
  })
  gr.snp <- Reduce(function(x, y) GenomicRanges::merge(x, y, all.x=T), 
                   append(bw.grlist, gr.dat))  
  
  
  
  # Fine-mapping tracks
  INIT_list <- list(   
    "GWAS"=ggbio::plotGrandLinear(obj = gr.snp, aes(y=-log10(P), fill=-log10(P))),
    "Fine_mapping"=ggbio::plotGrandLinear(obj = gr.snp, aes(y=mean.PP, color=mean.PP)) + 
      scale_color_viridis_c()
  )    
  # BigWig tracks from Glass
  bw.cols <- colnames(mcols(gr.snp))[1:12]  
  BW_tracks<- y(bw.cols, function(annot){
    gg.bw <- ggbio::plotGrandLinear(gr.snp, aes(y=eval(parse(text=annot)), 
                                                color=eval(parse(text=annot))) 
    ) + labs(y="Score", color="Score") + ylim(0,45)
    # if(bw.cols!=bw.cols[1]){
    #   gg.bw <- gg.bw + labs()
    # }
    return(gg.bw)                                                             
  }) 
  TRACKS_list <- append(INIT_list, BW_tracks)  
  names(TRACKS_list) <- c(names(INIT_list), gsub("_","\n",bw.cols))
  
  
  # gr.snp[start(gr.snp)!=Inf]
  # Fuse all tracks
  # gene <- "LRRK2"
  library(ggbio)
  params_list <- list(title = paste0(gene," [",length(seqnames(gr.snp))," SNPs]"), 
                      track.bg.color = "transparent",
                      track.plot.color = "transparent",
                      label.text.cex = .7, 
                      label.bg.fill = "grey12",
                      label.text.color = "white",
                      label.text.angle = 0,
                      label.width = unit(5.5, "lines"),
                      # xlim = c(min(start(gr.snp)), max(start(gr.snp))),
                      heights = c(rep(1,length(INIT_list)), rep(1,length(BW_tracks)) ))
  TRACKS_list <- append(TRACKS_list, params_list)
  trks <- suppressWarnings(do.call("tracks", TRACKS_list)) 
  
  # add lines
  lead.pos <- subset(finemap_DT,SNP=="rs76904798")$POS
  consensus.pos <- subset(finemap_DT, Consensus_SNP==T)$POS
  trks_plus_lines <- trks + 
    # theme_bw() + 
    # ggbio::theme_genome() +  
    # theme(strip.text.y = element_text(angle = 0), 
    #       strip.text = element_text(size=9 )) + 
    geom_vline(xintercept = lead.pos, color="red", alpha=1, size=.3, linetype='solid') +
    geom_vline(xintercept = consensus.pos, color="goldenrod2", alpha=1, size=.3, linetype='solid') 
  
  ggsave(file.path(results_path,"Annotation","Glass.snEpigenomics.png"), 
         plot = trks_plus_lines, dpi=400, height = 15, width = 8, 
         bg = "transparent")
  trks_plus_lines
  
  INIT_list[[1]] +
    geom_vline(xintercept = , color="red", alpha=1, size=.3, linetype='solid')   
  
  # Save merged Glass epigenomic data
  data.table::fwrite(GenomicRanges::as.data.frame(gr.snp[,bw.cols]),
                     file.path("./echolocatoR/annotations/Glass_lab/LRRK2.Glass.txt"),
                     sep="\t", nThread = 4)
}


NOTT_2019.superenhancers <- function(s6_path="./echolocatoR/annotations/Glass_lab/aay0793-Nott-Table-S6.xlsx"){ 
  s6 <- readxl::read_excel( , skip = 2)  
  annot_sub <- subset(s6, chr== paste0("chr",unique(finemap_DT$CHR)) & start>=min(finemap_DT$POS) & end<=max(finemap_DT$POS) )
  if(nrow(annot_sub)>0){
    merged_DT <- data.table:::merge.data.table(finemap_DT %>% 
                                                 dplyr::mutate(chr=paste0("chr",CHR), 
                                                               start=as.numeric(POS)) %>% 
                                                 data.table::data.table(),
                                               data.table::data.table(s6), 
                                               by = c("chr","start"))
  }
}

NOTT_2019.enhancer_promoter_interactions <- function(finemap_DT, 
                                                     s5_path="./echolocatoR/annotations/Glass_lab/aay0793-Nott-Table-S5.xlsx"){  
  sheets_s5 <- readxl::excel_sheets(s5_path) 
  s5 <- readxl::read_excel(s5_path, sheet = sheets_s5[1], skip = 2) 
  s5 <- s5 %>% dplyr::rename(chr=Chr, start=Start, end=End)  
  # Subset to window
  annot_sub <- subset(s5, chr== paste0("chr",unique(finemap_DT$CHR)) & 
                        start>=min(finemap_DT$POS) & 
                        end<=max(finemap_DT$POS) )  
  return(annot_sub)
}

NOTT_2019.get_promoters <- function(annot_sub, marker_key){  
  promoter.cols <- grep("*_active_promoter",  colnames(annot_sub), value = T)
  logical.list <- colSums(annot_sub[,promoter.cols])>0
  promoter_celltypes <- gsub("\\_.*","", promoter.cols[as.logical(logical.list)] )
  promoter_celltypes <- as.character(marker_key[promoter_celltypes])
  promoter_celltypes <- paste(promoter_celltypes,collapse="; ")
  return(promoter_celltypes)
}

NOTT_2019.get_interactions <- function(annot_sub, top.consensus.pos, marker_key){  
  interact.cols <- grep("*_interactions", colnames(annot_sub), value = T)   
  interact.DT <- lapply(interact.cols, function(column){
    coords <- strsplit(annot_sub[,column][[1]], ",")
    coord.dt <- lapply(coords, function(coord){
      data.table::data.table(Interaction=column, 
                             Cell_type=marker_key[gsub("\\_.*","",column)],
                             Coordinates=coord) %>% return() 
    }) %>% data.table::rbindlist()
    return(coord.dt)
  } )  %>% data.table::rbindlist()
  interact.DT <- subset(interact.DT, !is.na(Coordinates) & Coordinates!="") %>%  
    tidyr::separate(col = Coordinates, 
                    into=c("chr","Start","End"), sep = ":|-") %>%
    separate(col = Interaction, into=c("Marker","Element",NA), sep="_", remove = F ) 
  interact.DT <- interact.DT %>% 
    dplyr::mutate(Cell_type_interaction=paste(Cell_type,"-",Element))
  interact.DT$Cell_type <- interact.DT$Cell_type %>% as.character()
  interact.DT$Start <- as.numeric(interact.DT$Start)
  interact.DT$End <- as.numeric(interact.DT$End)
  # Summarise distance from different celltype enhancer interactions
  summarise_top.consensus.dist <- interact.DT %>% 
    dplyr::mutate(top.consensus.dist=End - top.consensus.pos) %>% 
    dplyr::group_by(Cell_type) %>% 
    dplyr::summarise(top.consensus.dist = mean(top.consensus.dist))
  print(summarise_top.consensus.dist)
  return(interact.DT)
}



# ***************************** #
GGBIO.nott_etal_2019 <- function(finemap_DT=NULL,
                                 locus="LRRK2",
                                 print_plot=T, 
                                 save_plot=T, 
                                 return_interaction_track=F){
  # finemap_DT <- quick_finemap()
  library(ggbio)
  marker_key <- list(PU1="Microglia",
                     Olig2="Oligodendrocytes",
                     NeuN="Neurons",
                     LHX2="Astrocytes")
  lead.pos <- top_n(finemap_DT, n = 1, wt = -P)$POS
  #subset(finemap_DT,SNP=="rs76904798")$POS  
  consensus.pos <- subset(finemap_DT, Consensus_SNP==T)$POS  
  if(length(consensus.pos)>0){
    top.consensus.pos <- (top_n(subset(finemap_DT, Consensus_SNP==T), 
                               n=1, wt = mean.PP) %>% top_n(1,wt=Effect))$POS[1]
  } else {
    top.consensus.pos <- (top_n(subset(finemap_DT, Support>0), 
                               n=1, wt = mean.PP )%>% top_n(1,wt=Effect))$POS[1]
  } 
  # Subset to relevant region
  annot_sub <- NOTT_2019.enhancer_promoter_interactions(finemap_DT = finemap_DT)
  ## Extract active promoters
  promoter_celltypes <- NOTT_2019.get_promoters(annot_sub = annot_sub, 
                                                marker_key = marker_key)
  ## Extract promoter interactions
  interact.DT <- NOTT_2019.get_interactions(annot_sub = annot_sub,
                                            top.consensus.pos =  top.consensus.pos, 
                                            marker_key = marker_key)
  
  
  # GWAS track 
  GWAS_trk <- ggplot(data=finemap_DT, aes(x=POS, y=-log10(P), color=-log10(P))) +
    geom_point()  
  
  FM_trk <- ggplot(data=finemap_DT, aes(x=POS, y=mean.PP, color=mean.PP)) +
    geom_point() +  
    scale_color_viridis_c(breaks=c(0,.5,1), limits=c(0,1)) +
    ylim(0,1)
  
  # Nott tracks   
  # Nott_s5 <- ggbio::ggbio() + 
  #   ggbio::geom_rect(data = annot_sub, aes(xmin=start, xmax=end, ymin=0, ymax=1), 
  #                    fill="turquoise", alpha=.75) +   
  #   facet_grid(facets = Annotation~.) 
  # Nott_s5 <- invisible_legend(Nott_s5)
  
  # Nott:  interactions
  # Nott_interactions 
  NOTT.interact_trk <- ggbio::ggbio() +
    ggbio::geom_arch(data = interact.DT, aes(x=Start, xend=End, color=Cell_type_interaction)) + 
    scale_y_reverse() +
    scale_colour_brewer(palette = "Accent") + 
    labs(subtitle = paste0(annot_sub$Annotation[[1]]," - ",promoter_celltypes) ) +
    theme(legend.key.width=unit(1.5,"line"),
          legend.key.height=unit(1.5,"line"))
  if(return_interaction_track){
    printer("++ Nott sn-epigenomics:: Returning PLAC-seq track only.")
    return(NOTT.interact_trk +  
             ggbio::geom_rect(data = annot_sub, 
                              aes(xmin=start, xmax=end, ymin=-0, ymax=Inf),
                              fill="turquoise", alpha=.5, inherit.aes=F) + 
             theme_bw())
  } else{
    # Makes tracks list
    TRACKS_list <- list(
      "GWAS"=GWAS_trk,
      "Fine-mapping"=FM_trk,
      "Nott (2019)\nInteractome"=NOTT.interact_trk
      # "Nott et al. (2019)\nPromoter\ninteractome"=Nott_s5 
    )
    # Parameters 
    params_list <- list(title = paste0(locus), 
                        track.bg.color = "transparent",
                        track.plot.color = "transparent",
                        label.text.cex = .7, 
                        label.bg.fill = "grey12",
                        label.text.color = "white",
                        label.text.angle = 0,
                        label.width = unit(5.5, "lines"),
                        xlim = c(min(finemap_DT$POS), max(finemap_DT$POS))
                        # heights = c(rep(1,length(INIT_list)), rep(1,length(BW_tracks)) )
    )
    TRACKS_list <- append(TRACKS_list, params_list)
    trks <- suppressWarnings(do.call("tracks", TRACKS_list)) 
    
    trks_plus_lines <- trks + 
      # Nott: promoter
      ggbio::geom_rect(data = annot_sub, aes(xmin=start, xmax=end, ymin=-0, ymax=Inf), 
                       fill="turquoise", alpha=.5, inherit.aes=F) +  
      # Lead GWAS line
      geom_vline(xintercept = lead.pos, color="red", alpha=1, size=.3, linetype='solid') +
      # Consensus line
      geom_vline(xintercept = consensus.pos, color="goldenrod2", alpha=1, size=.3, linetype='solid') + 
      theme_classic() +
      theme(plot.subtitle = element_text(color = "turquoise", size = 8))
    
    if(print_plot){print(trks_plus_lines)}
    # SAVE PLOT
    if(save_plot){
      plot.path <- file.path(results_path,"Multi-finemap",
                             paste0("Nott.sn-epigenomics_ggbio.png"))
      dir.create(dirname(plot.path), showWarnings = F, recursive = T)
      ggsave(filename = plot.path, 
             plot = trks_plus_lines,
             height = 7, width = 7, dpi = 1000, bg = "transparent")
    }
    return(trks_plus_lines)
  }
}



