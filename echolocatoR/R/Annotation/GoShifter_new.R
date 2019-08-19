
GS_create_snpmap <- function(snp_df, GS_results, verbose=T){
  printer("++ GoShifter: Creating snpmap file...", v=verbose)  
  snp_df %>%
    dplyr::rename(Chrom = CHR, BP = POS) %>% 
    dplyr::mutate(Chrom = paste0("chr",Chrom)) %>%
    dplyr::select(SNP, Chrom, BP) %>%
    data.table::fwrite(file.path(GS_results,"snpmap.txt"),sep="\t") 
  printer("+++ snpmap written to :",file.path(GS_results,"snpmap.txt"))
}



GS_create_LD <- function(results_path, verbose=T){
  printer("++ GoShifter: Creating LD file(s)...", v=verbose)
  # (chrA\tposA\trsIdA\tposB\trsIdB\tRsquared\tDPrime)
  GS_results_path <- file.path(results_path,"GoShifter")
  plink.ld <- data.table::fread(file.path(results_path,"plink","saved","plink.ld"))
  ld_file <- plink.ld %>% dplyr::rename(chrA = CHR_A, 
                                        posA = BP_A,
                                        rsIdA = SNP_A,
                                        posB = BP_B,
                                        rsIdB = SNP_B, 
                                        Rsquared = R,
                                        DPrime = DP) %>% 
    dplyr::mutate(Rsquared = Rsquared^2, chrA = paste0("chr",chrA)) %>% 
    dplyr::select(chrA,posA,rsIdA,posB,rsIdB,Rsquared,DPrime)
  # Create tabix file(s) 
  LD_folder <- file.path(GS_results_path,"LD")
  dir.create(LD_folder, showWarnings = F, recursive = T)
  for(chr in unique(ld_file$chrA)){
    ld_path <- file.path(LD_folder, paste0(chr,".EUR.tsv"))
    gz_path <- paste0(ld_path,".gz") 
    file.remove(gz_path) 
    data.table::fwrite(subset(ld_file, chrA==chr), ld_path,
                       sep="\t", quote = F, col.names = F) 
    # gzip(ld_path, destname = gz_path)
    system(paste("bgzip",ld_path))
    cmd <- paste("tabix",
                 "--begin 2",
                 "--end 2",
                 "--force",
                 gz_path)
    system(cmd) 
  }
  printer("+++ LD file(s) written to :",LD_folder)
}


GS_construct_reference <- function(Roadmap_reference = "./echolocatoR/tools/goshifter/annotation_files/ROADMAP_Epigenomic.js",
                                   EID_filter = NA,
                                   GROUP_filter = NA,
                                   ANATOMY_filter = NA,
                                   GR_filter = NA,
                                   fuzzy_search = NA){
  printer("++ GoShifter: Searching for Roadmap annotation BED files...", v=verbose)
  RM_ref <- suppressWarnings(data.table::fread(Roadmap_reference, skip = 1 , header = F, 
                                               col.names = c("EID",
                                                             "GROUP",
                                                             "ANATOMY",
                                                             "GR",
                                                             "Epigenome.Mnemonic",
                                                             "Standardized.Epigenome.name",
                                                             "Epigenome.name",
                                                             "TYPE")) ) 
  if(!is.na(EID_filter)){RM_ref <- subset(RM_ref, EID %in% EID_filter)}
  if(!is.na(GROUP_filter)){RM_ref <- subset(RM_ref, GROUP %in% GROUP_filter)}
  if(!is.na(ANATOMY_filter)){RM_ref <- subset(RM_ref, ANATOMY %in% ANATOMY_filter)}
  if(!is.na(GR_filter)){RM_ref <- subset(RM_ref, GR %in% GR_filter)}
  if(!is.na(fuzzy_search)){
    RM_ref <- RM_ref[(tolower(GROUP) %like% fuzzy_search) | (tolower(ANATOMY)  %like% fuzzy_search) | (tolower(GR) %like% fuzzy_search) | (tolower(Standardized.Epigenome.name) %like% fuzzy_search) | (tolower(Epigenome.name) %like% fuzzy_search) | (tolower(TYPE) %like% fuzzy_search )]
  } 
  total_found <- dim(RM_ref)[1]
  if(total_found > 0){
    printer("++ GoShifter: Found",dim(RM_ref)[1],"annotation BED files matching search query.", v=verbose)
  } else {stop("No annotation BED files found :(")} 
  return(RM_ref)
}


GS_bed_names <- function(RM_ref){ 
  # Create URLs
  bed_names <- lapply(unique(RM_ref$EID), function(eid, v=verbose){ 
    return(paste0(eid,"_15_coreMarks_mnemonics.bed.gz")) 
  }) %>% unlist() 
  return(bed_names)
}

# GS_download_roadmap_annotations_mod <- function(finemap_DT,
#                                             output_dir = "./echolocatoR/tools/Annotations", 
#                                             bed_names,
#                                             chromatin_state = "TssA",
#                                             verbose = T){
#  
#   # Download 
#   bed_subset_paths <- lapply(bed_names, function(bed, 
#                                              v=verbose, 
#                                              output_dir.=output_dir,
#                                              check_overlap=T){  
#     bed_url <- file.path("https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final",
#                          bed)
#     output_path <- file.path(output_dir., bed)
#     chrom_name <- "allStates"
#     eid <- strsplit(bed, "_")[[1]][1]
#     bed_subset <- file.path(dirname(output_path), paste0(eid,"_",chrom_name,"_subset.bed"))
#     bed_subset.gz <- paste0(bed_subset,".gz")
#     # Check if subset exists
#     if(file.exists(bed_subset.gz)){ 
#      printer("++ .gz file already exists.")
#     } else {
#       # Import and save FULL bed file
#       if(file.exists(output_path)){
#         printer("+++ Importing Previously downloaded file:",bed,v=v);
#         dat <- data.table::fread(output_path, col.names = c("Chrom", "Start","End","State") )
#       } else{ 
#         printer("+++ Downloading annotation BED file from Roadmap server...", v=v)
#         printer("   ", bed) 
#         dat <- data.table::fread(bed_url, col.names = c("Chrom", "Start","End","State") )
#         dir.create(dirname(output_path), showWarnings = F, recursive = T) 
#         data.table::fwrite(dat, output_path, col.names = F, sep = " ") # SPACE delimited
#       } 
#       try({ 
#         # Create SUBSET bed file
#         ########################################
#         ## Check if there's any overlap whatsoever 
#         FM.dat <- finemap_DT %>% dplyr::mutate(SEQnames = paste0("chr",finemap_DT$CHR))
#         gr.FM <- biovizBase::transformDfToGr(FM.dat, seqnames = "SEQnames", start = "POS", end = "POS")
#         gr.ANNO <-  biovizBase::transformDfToGr(dat, seqnames = "Chrom", start = "Start", end = "End")
#         library(GenomicRanges)
#         gr.overlap <- subsetByOverlaps(gr.ANNO, gr.FM)
#         dat <- as.data.frame(gr.overlap)
#         ########################################
#         ## Subset by chromatin state(s) 
#         if(!is.na(chromatin_state)){
#           printer("+++ Subsetting bed file: ",chromatin_state, v=v)
#           dat[, c("Num", "State") := tstrsplit(State, "_", fixed=TRUE)] # inplace transform
#           dat <- dat[State == chromatin_state,] 
#           chrom_name <- chromatin_state
#         }  
#         # If there's anything left, write it to a file
#         if(nrow(dat)>0){ 
#           printer("+++ Final bed dimensions:", paste( dim(dat), collapse = " x ")) 
#           data.table::fwrite(dat, bed_subset, col.names = F,sep = " ") # SPACE delimited
#           # zip it afterwards 
#           gzip(bed_subset, destname=bed_subset.gz, overwrite = T)
#         }  else {bed_subset.gz <- NULL}
#       }) # end try
#     } 
#     return(bed_subset.gz)
#   }) %>% unlist()
#   
#   # bed_subset_paths <- bed_subset_paths[!lapply(bed_subset_paths, is.null)%>%unlist()] %>% unique()
#   bed_subset_paths <- unique()
#   return(bed_subset_paths)
# }
GS_download_roadmap_annotations <- function(goshifter_path = "./echolocatoR/tools/goshifter", 
                                            bed_names,
                                            chromatin_state = "TssA",
                                            verbose = T){  
  # Term key for chromatin states
  # chromState_key <- data.table::fread(file.path(goshifter_path,
  #                                                 "annotation_files",
  #                                                 "ROADMAP_chromatinState_HMM.tsv"))  
  
  output_paths <- lapply(bed_names, function(bed, v=verbose){ 
    output_path <- file.path(goshifter_path,"annotation_files",bed)
    bed_url <- file.path("https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final",
                         bed)
    if(!file.exists(output_path)){
      printer("+++ Downloading annotation BED file from Roadmap server...", v=v)
      # download.file(bed_url, destfile = output_path)
      dat <- data.table::fread(bed_url, col.names = c("Chrom", "Start","End","State") )
      data.table::fwrite(dat, bed, col.names = F, sep = " ")
    } else{printer("+++ Importing Previously downloaded file:",bed,v=v);
      dat <- data.table::fread(output_path, col.names = c("Chrom", "Start","End","State") )
    }
    
    # Subset data to just the chromatin state(s) you want to check enrichment for
    printer("+++ Subsetting bed file: ",chromatin_state)
    dat[, c("Num", "State") := tstrsplit(State, "_", fixed=TRUE)] # inplace transform
    dat <- dat[State == chromatin_state,] 
    eid <- strsplit(bed, "_")[[1]][1]
    bed_subset <- file.path(dirname(output_path), paste0(eid,"_",chromatin_state,"_subset.bed"))
    data.table::fwrite(dat, bed_subset, col.names = F,sep = "\t")
    bed_subset.gz <- paste0(bed_subset,".gz")
    gzip(bed_subset, destname=bed_subset.gz, overwrite = T) 
    return(bed_subset.gz)
  }) %>% unlist() 
  return(output_paths)
}

GS_process_results <- function(RoadMap_ref, 
                               results_path, 
                               annotation_path){
  res <- data.table::fread(file.path(results_path,"GoShifter","GS_results.nperm1000.locusscore"),sep="\t")
  res$Annotation <- basename(annotation_path)
  eid <- strsplit(basename(annotation_path),"_")[[1]][1]
  res$EID <- eid 
  data <- subset(RoadMap_ref, EID == eid)
  res <- data.table:::merge.data.table(res, data.table::data.table(data), by="EID")
  
  # Overal enrichment
  overall_enrichment <- data.table::fread(file.path(results_path,"GoShifter","GS_results.nperm1000.enrich"),sep="\t") 
  overall_enrichment <- overall_enrichment[nperm == max(nperm),]
  res$nSnpOverlap <- overall_enrichment$nSnpOverlap
  res$allSnps <- overall_enrichment$allSnps
  res$enrichment <- overall_enrichment$enrichment 
  return(res)
}
 

GS_run <- function(RoadMap_ref,
                   results_path,
                   permutations = 1000, 
                   goshifter_path = "./echolocatoR/tools/goshifter",
                   chromatin_state = "TssA",
                   R2_filter = 0.8,
                   remove_tmps = T,
                   verbose = T){ 
  bed_names <- GS_bed_names(RoadMap_ref) 
  GS_results_path <- file.path(results_path,"GoShifter")
  GS_results <- lapply(bed_names, function(bed, 
                                           remove_tmps. = remove_tmps,
                                           results_path. = results_path,
                                           RoadMap_ref. = RoadMap_ref,
                                           R2_filter. = R2_filter){  
    
    printer("GoShifter: Testing for enrichment of SNPs within annotation:",basename(bed), v=verbose) 
    output_bed <- GS_download_roadmap_annotations(bed_names = bed, 
                                                  verbose = verbose, 
                                                  chromatin_state = chromatin_state)
    cmd <- paste(
      file.path(goshifter_path,"goshifter.py"),
      "--snpmap",file.path(GS_results_path,"snpmap.txt"),
      "--annotation",output_bed,
      "--permute",permutations,
      "--ld",file.path(GS_results_path,"LD"),
      "--out",file.path(GS_results_path,"GS_results"),
      "--rsquared",R2_filter. # Optional
      # "--window",window, # Optional
      # "--min-shift",min_shift, # Optional
      # "--max-shift",max_shift, # Optional
      # "--ld-extend",ld_extend, # Optional
      # "--nold" # Optional
    ) 
    system(cmd)  
    res <- GS_process_results(RoadMap_ref = RoadMap_ref.,
                              results_path = results_path., 
                              annotation_path = output_bed)
    res$chromatin_state <- chromatin_state
    if(remove_tmps.){ file.remove(output_bed)}
    return(res)
  }) %>% data.table::rbindlist() 
  return(GS_results)
}







####----- Go Shifter: Main Function -----####

GoShifter <- function(goshifter_path = "./echolocatoR/tools/goshifter",
                      results_path, 
                      CredibleSet_only = T,
                      permutations = 1000,
                      Roadmap_chromatinMarks_search = "",
                      chromatin_state = NA,
                      R2_filter = 0.8, # Include LD SNPs at rsquared >= NUM [default: 0.8]
                      remove_tmps = T,
                      verbose = T,
                      save_results = T){ 
  # results_path <-  "./Data/GWAS/Nalls23andMe_2019/LRRK2"; Roadmap_chromatinMarks_search <- "monocytes"; CredibleSet_only = T; permutations = 1000; chromatin_state = NA;  R2_filter = 0.8
  # Cleanup pyc tmp files
  suppressWarnings(file.remove(
    file.path(goshifter_path, 
              c("chromtree.pyc",
                "data.pyc",
                "docopt.pyc",
                "functions.pyc")) ) )
  # Term key for chromatin states
  chromState_key <- data.table::fread(file.path(goshifter_path,
                                                  "annotation_files",
                                                  "ROADMAP_chromatinState_HMM.tsv"))

  # Create GoShifter folder
  # results_path="Data/GWAS/Nalls23andMe_2019/LRRK2"
  GS_results_path <- file.path(results_path,"GoShifter")
  dir.create(GS_results_path, showWarnings = F, recursive = T) 
  
  # snpmap file 
  snp_df <- data.table::fread(file.path(results_path,"Multi-finemap/Multi-finemap_results.txt"))
  if(CredibleSet_only){Support_level <- 1} else {Support_level <- 0}
  snp_df <- subset(snp_df, Support >= Support_level) 
  GS_create_snpmap(snp_df, GS_results_path, verbose = verbose)
  # OR, copy file directly from example
  # file.copy(file.path(goshifter_path,"test_data","bc.snpmappings.hg19.txt"), 
  #           file.path(GS_results_path,"snpmap.txt"),overwrite = T)
  
  # LD file 
  GS_create_LD(results_path, verbose = verbose) 
  # Gather Annotation BED files 
  RoadMap_ref <- GS_construct_reference(fuzzy_search = Roadmap_chromatinMarks_search)
  # Run
  GS_results <- GS_run(goshifter_path = goshifter_path,
                       RoadMap_ref = RoadMap_ref, 
                       results_path = results_path, 
                       chromatin_state = chromatin_state,
                       R2_filter = R2_filter,
                       permutations = permutations, 
                       remove_tmps = remove_tmps,  
                       verbose = verbose) 
  message("GoShifter results data.frame: ",nrow(GS_results), " rows x ", ncol(GS_results)," cols" )
  # Summarise results
  sig_results <- GS_results %>%  
    # Apply bonferonni correction
    # subset(enrichment <= 0.05/nrow(GS_results) ) %>%  
    arrange(desc(enrichment), desc(nSnpOverlap)) %>%
    dplyr::mutate(enrichment = formatC(enrichment,format="e", digits=7) )  
  createDT(sig_results) %>% print()
  
  # Get sig enrichments per tissue per chomatin_state
  sig_count <- sig_results %>% 
    dplyr::group_by(EID) %>% 
    slice(1) %>%
    dplyr::group_by(chromatin_state) %>% 
    count() 
  printer("+++",sig_count$n,"/",nrow(RoadMap_ref),
          "tissues had significant enrichment for the chomatin state(s):",
          paste(chromatin_state, collapse = ", "), v = verbose)
  
  if(save_results){ 
    data.table::fwrite(GS_results, file.path(GS_results_path,"GS_results.txt"),
                       sep="\t", quote = F)
  } 
  if(remove_tmps){
    # Remove GS output files to avoid reading in old files in future runs
    suppressWarnings(
      file.remove(file.path(GS_results_path,paste0("GS_results",c(".nperm1000.enrich",
                                                                  ".nperm1000.locusscore",
                                                                  ".nperm1000.snpoverlap"))) ) )
  }
  return(GS_results)
}

 