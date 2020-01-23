
# DeepBlueR: Alternate source of annotation files
## https://bioconductor.org/packages/release/bioc/vignettes/DeepBlueR/inst/doc/DeepBlueR.html#listing-experiments


# TUTORIALS
# http://xgr.r-forge.r-project.org/#xenrichergenes
# http://galahad.well.ox.ac.uk:3020/R/ds2

# SNP-based enrichment analysis using GWAS Catalog traits (mapped to EF) 
# data <- merged_results$SNP %>% unique()

xgr_snpEnrich <- function(data_path="./Data/GWAS/Nalls23andMe_2019/LRRK2/Multi-finemap/Multi-finemap_results.txt"){
  library(XGR) 
  fullDat <- data.table::fread(data_path, sep="\t", stringsAsFactors = F)
  data <- fullDat[P<=5e-8,]$SNP %>% unique()
  
  eTerm <- xEnricherSNPs(data=data, 
                         ontology="EF",
                         path.mode=c("all_paths"), 
                         min.overlap = 1)
  xEnrichViewer(eTerm)
  # e) barplot of significant enrichment results
  bp <- xEnrichBarplot(eTerm, top_num="auto", displayBy="adjp")
  print(bp)
  # f) visualise the top 10 significant terms in the ontology hierarchy
  # color-code terms according to the adjust p-values (taking the form of 10-based negative logarithm)
  xEnrichDAGplot(eTerm, top_num=10, 
                 displayBy="adjp",
                 node.info=c("full_term_name"))
  # color-code terms according to the z-scores
  xEnrichDAGplot(eTerm, top_num=10, 
                 displayBy="zscore",
                 node.info=c("full_term_name"))
  # Circos plot
  library(RCircos)
  SNP.g <- xSocialiserSNPs(data, include.LD=NA, LD.r2 = 1) 
  xCircos(g=SNP.g, 
          entity="SNP")
  # b') optionally, enrichment analysis for input SNPs plus their LD SNPs
  ## LD based on European population (EUR) with r2>=0.8
  # eTerm_LD <- xEnricherSNPs(data=data, 
  #                        include.LD="EUR", 
  #                        LD.r2=0.8,  
  #                        min.overlap = 1)
  # xEnrichViewer(eTerm_LD)
}

xgr_geneEnrich <- function(top_SNPs){
  library(XGR) 
  # GENE-LEVEL ENRICHMENT 
  data <- top_SNPs$Gene %>% unique()
  eTerm_gene <- xEnricherGenes(data = data, 
                               ontology = "MsigdbC2REACTOME", 
                               min.overlap = 1) 
  bp <- xEnrichBarplot(eTerm_gene, top_num="auto", displayBy="adjp")
  print(bp)
  
  xEnrichDAGplot(eTerm, top_num=10, displayBy="zscore",
                 node.info=c("full_term_name"), graph.node.attrs=list(fontsize=200)) 
  } 


xgr_annoEnrich <- function(data_path="./Data/GWAS/Nalls23andMe_2019/LRRK2/Multi-finemap/Multi-finemap_results.txt"){
  library(XGR) 
  RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
  fullDat <- data.table::fread(data_path, sep="\t", stringsAsFactors = F)
  CS <- subset(fullDat, )
  
  ## a) provide input data
  # data.file <- "http://galahad.well.ox.ac.uk/bigdata/ImmunoBase_GWAS.bed" 
  data.file_full <- fullDat %>% dplyr::mutate(chrom = paste0("chr",CHR), 
                                         chromStart = POS, 
                                         chromEnd = POS, 
                                         name = SNP) %>% 
    dplyr::select(chrom, chromStart, chromEnd, name)
  data.file_CS <- CS %>% dplyr::mutate(chrom = paste0("chr",CHR), 
                                            chromStart = POS, 
                                            chromEnd = POS, 
                                            name = SNP) %>% 
    dplyr::select(chrom, chromStart, chromEnd, name)
  
  
  
  # Built-in database
  ## b) perform enrichment analysis using FANTOM expressed enhancers
  ### one-tail p-value calculation (by default) 
  eTerm_CS <- xGRviaGenomicAnno(data.file = data.file_CS, 
                                background.file = data.file_full,
                                format.file = "data.frame",
                                GR.annotation = "FANTOM5_Enhancer_Cell")
  xEnrichViewer(eTerm_CS, 10)
  # Downloaded database
  ## Download
  gr_broad <- xRDataLoader(RData.customised="Broad_Histone")
  
  gr_reactome <- xRDataLoader(RData.customised ="org.Hs.egMsigdbC2REACTOME")
  
  
  # Run
  eTerm_broad <- xGRviaGenomicAnno(data.file = data.file_CS, 
                                     background.file = data.file_full,
                                     format.file="data.frame", 
                                     GR.annotation = gr_broad)
  subset(xEnrichViewer(eTerm_broad), adjp<=0.05) 
  
  
  # Gtex
  ImmunoBase <- xRDataLoader(RData.customised='ImmunoBase')
  grl <- xDefineGenomicAnno("Uniform_TFBS")
  
  gr_gtex <- xRDataLoader(RData.customised="org.Hs.egGTExV6")
  eTerm_ImmunoBase <- xGRviaGenomicAnno(data.file = data.file_CS, 
                                   background.file = data.file_full,
                                   format.file="data.frame",  
                                   GR.annotation = grl)
  subset(xEnrichViewer(eTerm_broad), adjp<=0.05) 
}  


## Download BED files via XGR
GRs.to.BED <- function(GR.annotations, output_path, sep="\t"){
  BED_paths <- lapply(names(GR.annotations), function(name){ 
    GR <- GR.annotations[[name]] 
    BED <- GR %>% as.data.table() %>% 
      dplyr::select(chrom=seqnames,
                    chromStart=start,
                    chromEnd=end,
                    strand) 
    BED_path <- file.path(output_path,paste0(gsub(":","-",name),".bed.txt"))
    dir.create(dirname(BED_path), recursive = T, showWarnings = F)
    data.table::fwrite(BED, BED_path, sep=sep, col.names = F, quote = F) 
    return(BED_path)
  }) %>% unlist()
  return(BED_paths)
}


xgr_annoEnrich_mass <- function(foreground_snps,
                                background_snps=NULL,
                                lib.selections = c("ENCODE_TFBS_ClusteredV3_CellTypes",
                                                   "ENCODE_DNaseI_ClusteredV3_CellTypes", 
                                                    # "Broad_Histone",
                                                    "FANTOM5_Enhancer",
                                                    "Segment_Combined_Gm12878",
                                                    "TFBS_Conserved",
                                                    "ReMap_PublicAndEncode_TFBS",
                                                    "Blueprint_VenousBlood_Histone",
                                                    "Blueprint_DNaseI",
                                                    # "Blueprint_Methylation_hyper",
                                                    # "Blueprint_Methylation_hypo",
                                                    # "Genic_anno",
                                                    "FANTOM5_CAT_Cell",
                                                    "FANTOM5_CAT_MESH",
                                                    "GWAScatalog_alltraits"),  
                                save_path=F){ 
  # Description of all datasets
  # https://www.rdocumentation.org/packages/XGR/versions/1.1.5/topics/xDefineGenomicAnno 
  
  data.file <- foreground_snps
  data.file_full <- background_snps
  
  database_results <- lapply(lib.selections, function(lib.name){ 
    tryCatch({
      # lib.name <- "ENCODE_TFBS_ClusteredV3_CellTypes"
      message("+ Testing enrichment: ",lib.name)
      # Import annotation files from server
      GR.annotations <- xRDataLoader(lib.name,
                                     RData.location=RData.location)
      ## Un-nest files list (e.g. some are organized per tissue per TF)
      GR.annotations <- unlist(GR.annotations)
      # Iterate over each file
      ls_df <- lapply(1:length(GR.annotations), function(i,
                                                         library = lib.name, 
                                                         data.file. = data.file){
        GR.annotation <- GR.annotations[i]
        message(sprintf("Analysing '%s' (%s) ...", names(GR.annotation),
                        as.character(Sys.time())), appendLF=T)
        df <- xGRviaGenomicAnno(data.file=data.file., 
                                background.file = data.file_full,
                                format.file="data.frame",
                                GR.annotation=GR.annotation, 
                                RData.location=RData.location, 
                                verbose=F)
        df$library <- library
        df$total_snps <- nrows(data.file.)
        return(df)
      })
      df <- do.call(rbind, ls_df) %>% data.table::as.data.table()
      return(df) 
    }, 
    error=function(e){data.table::data.table(name=NA, nAnno=NA, 
                                             nOverlap=NA,fc=NA,
                                             zscore=NA,pvalue=NA,
                                             adjp=NA,or=NA,
                                             CIl=NA,CIu=NA,
                                             total_snps=NA,
                                             library=lib.name)}) 
  })
  # Process results
  DT <- database_results %>% 
    data.table::rbindlist()
  # %>% subset(!is.na(name))
  # bonf <- 0.05 / nrow(DT)
  # subset(DT, adjp<=bonf & nOverlap > 1) %>% arrange(desc(nOverlap))
  
  if(save_path!=F){
    data.table::fwrite(DT, save_path, quote = F)
  }
  # Heatmap
  # gp <- xEnrichHeatmap(database_results, fdr.cutoff=0.05, displayBy="fdr",
  #                      reorder="both") 
  # gp
  # # Barplot
  # bp <- xEnrichBarplot(eTerm, top_num='auto', displayBy="fc")
  # bp
  # # Forest plot
  # gp <- xEnrichForest(eTerm)
  # gp
  return(DT)
}
#   
# gather_databases <- function(lib.selections){
#   databases <- lapply(lib.selections, function(lib.name){ 
#     tryCatch({
#       xRDataLoader(lib.name)
#     }, 
#     error=function(e) NULL) 
#   })
#   names(databases) <- paste0(lib.selections,"..")
#   dbs <- unlist(databases)
#   names(dbs)
#   return(dbs)
# }
# GR.annotations <- gather_databases(lib.selections)




