# [GitHub](https://github.com/gkichaev/PAINTOR_V3.0)
# [LD Tutorial](https://github.com/gkichaev/PAINTOR_V3.0/wiki/2a.-Computing-1000-genomes-LD)

PAINTOR.install <- function(paintor_path="./echolocatoR/tools/PAINTOR_V3.0/"){
  cmd <- paste("cd", paintor_path,"&& bash install")
  system(cmd)
}

PAINTOR.help_menu <- function(paintor_path="./echolocatoR/tools/PAINTOR_V3.0/"){
  system(file.path(paintor_path,"PAINTOR"))
}

zscore <- function(vec){
  z <- scale(vec, center = T, scale = T)
  return(z)
}


Zscore.get_mean_and_sd <- function(fullSS="./Data/GWAS/Nalls23andMe_2019/nallsEtAl2019_allSamples_allVariants.mod.txt",
                                   Effect_col="beta",
                                   use_saved=T,
                                   output_path="./Data/GWAS/Nalls23andMe_2019/z.info.RDS"){
  if(use_saved & file.exists(output_path)){
    printer("Reading in:",output_path,"...")
    z.info <- readRDS(output_path) 
  } else { 
    printer("Extracting mean and standard deviation from",fullSS,"...")
    sample_x <- data.table::fread(fullSS, nThread = 4, select=c(Effect_col))
    sample.mean <- mean(sample_x$beta, na.rm = T)
    sample.stdv <- sd(sample_x$beta)
    z.info <- list(file.name=fullSS,
                   colname=Effect_col,
                   sample.mean=sample.mean, 
                   sample.stdv=sample.stdv,
                   sample.min=min(sample_x$beta),
                   sample.max=max(sample_x$beta)) 
    saveRDS(z.info, file = output_path )
  } 
  return(z.info)
}


Zscore <- function(x, z.info){ 
  # Need to use the mean and standard deviation of the FULL dataset (i.e. all beta fomr the full summary stats file)
  sample.stdv <- z.info$sample.stdv
  sample.mean <- z.info$sample.mean 
  z <- (x - sample.mean) / sample.stdv 
  return(z)
}

PAINTOR.create_locusFile <- function(results_path,
                                     PT_results_path,
                                     GWAS_dataset_name="Nalls23andMe_2019",
                                     QTL_datasets=c("Fairfax_2014_CD14","Fairfax_2014_IFN",
                                                    "Fairfax_2014_LPS2","Fairfax_2014_LPS24"), #NULL
                                     gene,
                                     locus_name){
  printer("++ Creating Locus File...")
  # The locus file should at the very minimum contain the Z-scores for all the populations,
  # though metadata on each SNP such as chromosome, position, and rsid are recommended.
  # The top line of the locus file must contain a header with the names of the fields.
  # The Z-score of a SNP is the Wald statistic, obtained from standard regression of
  # the phenotype onto the SNP.

  # Import GWAS data 
  GWAS_path <- file.path(results_path,"Multi-finemap/Multi-finemap_results.txt") 
  z.info.gwas <- Zscore.get_mean_and_sd(fullSS="./Data/GWAS/Nalls23andMe_2019/nallsEtAl2019_allSamples_allVariants.mod.txt",
                                        Effect_col="beta", 
                                        use_saved=T,
                                        output_path="./Data/GWAS/Nalls23andMe_2019/z.info.RDS")
  # Get the party started using the GWAS file
  finemap_DT <- data.table::fread(GWAS_path, nThread = 4) %>%
    dplyr::mutate(CHR=paste0("chr",CHR), 
                  RSID=SNP,
                  ZSCORE.P1=Zscore(x = Effect, z.info = z.info.gwas))
  merged_DT <- finemap_DT 
  
  # Merge QTL data (for loop won't run if QTL_datasets=NULL)
  for(i in 1:length(QTL_datasets)){
     printer("PAINTOR:: Merging QTL data - ",QTL_datasets[i])
     z.info.qtl <- Zscore.get_mean_and_sd(fullSS=Directory_info(dataset_name = QTL_datasets[i], variable = "fullSumStats"),
                                          Effect_col="beta", 
                                          use_saved=T,
                                          output_path=file.path("./Data/QTL/Fairfax_2014",paste0("z.info.",QTL_datasets[i],".RDS")) )
    # Merge QTL data together
    merged_DT <- mergeQTL.merge_handler(FM_all = merged_DT, qtl_file = QTL_datasets[i])
    merged_DT <- dplyr::mutate(merged_DT, 
                               QTL.ZSCORE = Zscore(x=QTL.Effect, z.info = z.info.qtl ))
    names(merged_DT)[names(merged_DT) == "QTL.ZSCORE"] <- paste0("ZSCORE.P",i+1)
    QTL.cols <- grep("QTL.",colnames(merged_DT), value = T)
    merged_DT <- dplyr::select(merged_DT, -QTL.cols)
  } 
  # Post-processing
  ## Subset to only necessary columns
  z.cols <- grep("ZSCORE.P",colnames(merged_DT), value=T)
  merged_DT <- subset(merged_DT, select=c("RSID","CHR","POS",z.cols))
  ## Remove NAs (PAINTOR doesn't tolerate missing values)
  merged_DT <- merged_DT[complete.cases(merged_DT),]
  ## Save 
  data.table::fwrite(merged_DT, 
                     file.path(PT_results_path, locus_name),
                     sep=" ", quote = F, na = NA, nThread = 4) 
  return(merged_DT)
}



PAINTOR.prepare_LD <- function(results_path,
                               PT_results_path,
                               locus_name,
                               locus_DT, 
                               gene){
  ### "VERY IMPORTANT!
  # Particular care must be taken when computing LD from a reference panel such the 1000 genomes.
  ## It is imperative that all the reference and alternate alleles for SNPs from which the Z-scores
  ## were computed match the reference and alternate alleles of the reference panel.
  ## The output of PAINTOR will not be correct if there are mismatches of this type in the data.
  ##Please see wiki section 2a for instructions on how to use the LD utility provided with the software."
  printer("++ PAINTOR:: Creating LD Matrix File...")
  finemap_DT <- data.table::fread(file.path(results_path,"Multi-finemap/Multi-finemap_results.txt"), nThread = 4)
  LD_matrix <- LD.load_or_create(results_path=results_path,
                                subset_DT=finemap_DT,
                                gene=gene,
                                download_reference = T,
                                LD_reference="1KG_Phase1",
                                superpopulation="EUR")
  
  load(file.path(results_path,"plink/LD_matrix.RData"))
  ## Make sure SNPs are in the same order as the Locus File
  .LD_file <- file.path(PT_results_path, paste0(locus_name,".ld1"))
  ld.mat <- LD_matrix[locus_DT$RSID, locus_DT$RSID] %>%
    data.table::as.data.table()
  ## Write
  printer("++ PAINTOR:: Writing LD file to ==> ",.LD_file)
  data.table::fwrite(ld.mat,
                     .LD_file,
                     sep = " ", quote = F,
                     col.names = F, row.names = F)
}






PAINTOR.locusName_handler <- function(locus_name=NA, 
                                      gene, 
                                      GWAS_dataset_name=NA, 
                                      QTL_datasets=NA){
  if(is.na(locus_name)){ 
    locus_name <- paste(gene,GWAS_dataset_name,paste(QTL_datasets,collapse = ":"), sep = ".")
    return(locus_name)
  } 
  return(locus_name)
}

PAINTOR.datatype_handler <- function(GWAS_dataset_name=NA, 
                                     QTL_datasets=NA, 
                                     gene){
  if(!is.na(GWAS_dataset_name) & !all(is.na(QTL_datasets))){
    printer("++ PAINTOR:: GWAS and QTL input data detected. Feeding both into PAINTOR...")
    results_path <- file.path("./Data/GWAS",GWAS_dataset_name, gene) 
  } else if(!is.na(GWAS_dataset_name) & all(is.na(QTL_datasets))){
      printer("++ PAINTOR:: Only GWAS input data detected. Feeding into PAINTOR...")
      results_path <- file.path("./Data/GWAS",GWAS_dataset_name, gene)
    } else if(is.na(GWAS_dataset_name) & !all(is.na(QTL_datasets))){
      printer("++ PAINTOR:: Only QTL input data detected. Feeding into PAINTOR...")
      results_path <- file.path("./Data/QTL",paste0(QTL_datasets,collapse=":"),"merged_results", gene)
    } else {
      stop("++ PAINTOR:: Neither GWAS nor QTL data detected. Please enter at least one valid dataset.")
    } 
  PT_results_path <- file.path(results_path,"PAINTOR")
  dir.create(PT_results_path, showWarnings = F, recursive = T)
  printer("++ PAINTOR:: Results will be stored in  ==> ", PT_results_path) 
  return(PT_results_path)
} 



PAINTOR.download_annotations <- function(PT_results_path, 
                                         locus_name,
                                         locus_DT,
                                         XGR_dataset=NA,
                                         ROADMAP_search=NA,
                                         chromatin_state="TssA",
                                         no_annotations=F){
  printer("++ PAINTOR:: Creating Annotation Matrix File...") 
  .annotations_file <- file.path(PT_results_path, paste0(locus_name,".annotations.txt"))
  if(no_annotations){
    printer("+++ PAINTOR:: no_annotations=T; Prior Probability set to 1 for all SNPs.")
    data.table::fwrite(data.frame(Default_Priors = rep(1,nrow(locus_DT))),
                       .annotations_file,
                       sep=" ", quote = F)
    BED_paths.gz <- NA
  } else{
    ## Download GRs and convert to BED
    if(!is.na(XGR_dataset)){
      printer("+++ PAINTOR:: Gathering annotations via XGR for:",paste(XGR_dataset,collapse=", "))
      GR.annotations <- XGR::xRDataLoader(XGR_dataset)
      printer("+++ PAINTOR:: Writing BED file(s) to  ==> ",file.path(PT_results_path,"annotation_files"))
      BED_paths.gz <- GRs.to.BED(GR.annotations = GR.annotations,
                              output_path = file.path(PT_results_path,"annotation_files"),
                              sep=" ")
    }
    if(!is.na(ROADMAP_search)){
      printer("+++ PAINTOR:: Gathering annotations via ROADMAP server for:",paste(ROADMAP_search,collapse=", "))
      # Gather roadmap annotations
      if(ROADMAP_search=="all_cell_types"){ROADMAP_search <- ""}
      # Version 1: Subset locally 
      RoadMap_ref <- GoShifter.search_ROADMAP(fuzzy_search = ROADMAP_search)
      bed_names <- GoShifter.bed_names(RoadMap_ref)
      BED_paths.gz <- GoShifter.get_roadmap_annotations(bed.list = bed_names, 
                                                        chromatin_state = chromatin_state, 
                                                        verbose = T)  
    }
  }
  return(BED_paths.gz)
}


PAINTOR.list_paintor_annotations <- function(annotations_dir=file.path("/sc/orga/projects/pd-omics/brian/PAINTOR_V3.0",
                                                                       "Annotation_directory/Functional_Annotations")){
  types <- basename(list.dirs(annotations_dir, recursive = F))
  printer("PAINTOR:: Available annotations provided by PAINTOR_V3.0...")
  for(atype in types){printer(atype)} 
}


PAINTOR.select_paintor_annotations <- function(annotations_dir=file.path("/sc/orga/projects/pd-omics/brian/PAINTOR_V3.0",
                                                                         "Annotation_directory/Functional_Annotations"),
                                               annotation_types=c("Roadmap_ChromeHMM_15state",
                                                                  "RoadMap_Enhancers",
                                                                  "RoadMap_Promoter",
                                                                  "TFBS")){
  selected_annotations <- list.files(file.path(annotations_dir, annotation_types), full.names = T)
  printer("PAINTOR::",length(selected_annotations),"annotations pulled.")
  return(selected_annotations)
}



PAINTOR.prepare_annotations <- function(paintor_path="./echolocatoR/tools/PAINTOR_V3.0",
                                        BED_paths,
                                        PT_results_path,
                                        locus_name){
  .locus_file <- file.path(PT_results_path, locus_name)
  .annotations_file <- file.path(PT_results_path, paste0(locus_name,".annotations.txt"))
  
  printer("++ PAINTOR:: Merging and formatting BED files using PAINTOR utility.")
  ## Use PAINTOR utility to format merge BED files
  printer("+++ PAINTOR:: Decompressing BED files.")
  gz.files <- grep("*.gz",BED_paths,value = T)
  if(length(gz.files)>0){
     for (gz in gz.files){
      # Unzip but keep original files
      try({gunzip(gz, overwrite=T, remove=F)})
    }
  }
 
  BED_paths <- gsub("*.gz","", BED_paths) 
  annotation_paths <- file.path(PT_results_path,"annotation_paths.txt")
  # Wait until the BED files are decompressed before trying to write them to a file
  while(any(!file.exists(BED_paths))){
    Sys.sleep(.001)
  }
  printer("+++ PAINTOR:: Writing annotations paths file to  ==> ",annotation_paths)
  data.table::fwrite(list(BED_paths),
                     annotation_paths,
                     sep="\n")
  cmd <-  paste("python",file.path(paintor_path,"PAINTOR_Utilities","AnnotateLocus.py"),
                "--input", file.path(PT_results_path,"annotation_paths.txt"),
                "--locus", .locus_file,
                "--out", .annotations_file,
                "--chr","CHR",
                "--pos","POS")
  system(cmd)
  printer("+++ PAINTOR:: Annotation--SNP overlap summaries:")
  colSums( data.table::fread(.annotations_file))
  printer("+++ PAINTOR:: Removing temporary decompressed BED files.")
  file.remove(BED_paths)
}



PAINTOR.survey_annotation <- function(PT_results_path, 
                                      locus_name="Locus1"){
  anno <- data.table::fread(file.path(PT_results_path, paste0(locus_name,".annotations.txt")))
  column_sums <- sort(colSums(anno), decreasing = T)
  signals <- column_sums[column_sums > 0]
  print(signals)
}

PAINTOR.process_results <- function(PT_results_path, 
                                    locus_name="Locus1"){
  results <- data.table::fread(file.path(PT_results_path, paste0(locus_name,".results.txt"))) %>% 
    arrange(desc(Posterior_Prob))
  # subset(results, Posterior_Prob>0) 
  return(results)
}


PAINTOR.run <- function(paintor_path,
                        PT_results_path,
                        n_datasets,
                        n_causal=5){
  printer("+ PAINTOR:: Running PAINTOR...")
  ## RUN
  # https://github.com/gkichaev/PAINTOR_V3.0/wiki/3.-Running-Software-and-Suggested-Pipeline
  cmd <- paste(
    file.path(paintor_path,"PAINTOR"),
    
    #### REQUIRED ####
    # (required) Filename of the input file containing the
    ## list of the fine-mapping loci [default: N/A]
    "-input",file.path(PT_results_path,"input.files"),
    
    #  (required) The name(s) of the Zscore column
    ## in the header of the locus file (comma separated) [default: N/A]
    "-Zhead", paste(paste0("ZSCORE.P",1:n_datasets), collapse=","),
    
    # (required) Suffix(es) for LD files. Must match the order of
    ## Z-scores in which the -Zhead flag is specified (comma separated) [Default:N/A]
    "-LDname", paste(paste0(rep("ld",n_datasets),1), collapse=","),
    
    # specify this flag if you want to enumerate all possible configurations
    ## followed by the max number of causal SNPs (eg. -enumerate 3 considers
    ## up to 3 causals at each locus) [Default: not specified]
    # "-enumerate",n_causal,
    
    #  should the algorithm be run with MCMC? [Default: not specified]
    "-mcmc",
    
    #### OPTIONAL ####
    # The names of the annotations to include in model (comma separated)
    ## [default: N/A]
    # "-annotations",paste(basename(bed),collapse=","),
    
    # Input directory with all run files [default: ./ ]
    "-in", paste0(PT_results_path),
    
    #  Output directory where output will be written [default: ./ ]
    "-out",file.path(PT_results_path),
    
    # Output Filename for enrichment estimates [default: Enrichment.Estimate]
    "-Gname","Enrichment.Estimates.txt",
    
    # Suffix for ouput files of results [Default: results]
    "-RESname","results.txt",
    
    # Suffix for annotation files [Default: annotations]
    "-ANname","annotations.txt",
    
    "-set_seed","2019"
  )
  printer(cmd)
  system(cmd)
  # }
  # EXAMPLE
  # cmd <- paste("cd",paintor_path,"&& ./PAINTOR -input SampleData/input.files -in SampleData/ -out SampleData/ -Zhead Zscore -LDname ld -enumerate 2 -annotations DHS")
  # system(cmd)
}

# ------------------------------#
# ----------- PAINTOR ----------#
# ------------------------------#
PAINTOR <- function(paintor_path = "./echolocatoR/tools/PAINTOR_V3.0",
                    GWAS_dataset_name=NA,
                    QTL_datasets=NA,#c("Fairfax_2014_CD14","Fairfax_2014_IFN", "Fairfax_2014_LPS2","Fairfax_2014_LPS24")
                    gene,
                    locus_name=NA,
                    n_causal=5,
                    XGR_dataset="FANTOM5_Enhancer",
                    ROADMAP_search=NA,
                    chromatin_state="TssA",
                    no_annotations=F){
  # Note: All file formats are assumed to be single space delimited.
  
  ## Quick setup
  # chromatin_state="TssA"; paintor_path = "./echolocatoR/tools/PAINTOR_V3.0"; GWAS_dataset_name = "Nalls23andMe_2019"; QTL_dataset_name = "Fairfax_2014_CD14";locus=NA; gene="LRRK2";  no_annotations=F;  XGR_dataset=NA; ROADMAP_search="monocyte"
  locus_name <- PAINTOR.locusName_handler(locus_name, 
                                          gene, 
                                          GWAS_dataset_name, 
                                          QTL_datasets)
  
  PT_results_path <- PAINTOR.datatype_handler(GWAS_dataset_name=GWAS_dataset_name, 
                                              QTL_datasets=QTL_datasets, 
                                              gene=gene)
  printer("****** Double checking PT_results_path",PT_results_path)
  results_path <- dirname(PT_results_path)


  # 1. Locus File 
  locus_DT <- PAINTOR.create_locusFile(results_path=results_path,
                                       PT_results_path=PT_results_path,
                                       GWAS_dataset_name=GWAS_dataset_name,
                                       QTL_datasets=QTL_datasets,
                                       gene=gene,
                                       locus_name=locus_name) 
  n_datasets <- sum(grepl(colnames(locus_DT), pattern = "ZSCORE"))

  # 2. LD Matrix File
  PAINTOR.prepare_LD(results_path,
                     PT_results_path,
                     locus_name,
                     locus_DT)
  
  # 3. Annotation Matrix File
  # BED_paths.gz <- PAINTOR.download_annotations(PT_results_path, 
  #                                              locus,
  #                                              locus_DT,
  #                                              XGR_dataset,
  #                                              ROADMAP_search,
  #                                              chromatin_state,
  #                                              no_annotations)
 

  if(no_annotations==F){
    PAINTOR.list_paintor_annotations() 
    BED_paths <- PAINTOR.select_paintor_annotations(annotations_dir=file.path("/sc/orga/projects/pd-omics/brian/PAINTOR_V3.0",
                                                                              "Annotation_directory/Functional_Annotations"),
                                                    annotation_types=c("Roadmap_ChromeHMM_15state",
                                                                        "RoadMap_Enhancers",
                                                                        "RoadMap_Promoter",
                                                                        "TFBS"))
    PAINTOR.prepare_annotations(paintor_path,
                                BED_paths.gz, 
                                PT_results_path, 
                                locus)
  }
  
  # 4. Input File
  printer("+ PAINTOR:: Preparing input.files")
  data.table::fwrite(list(locus_name), file.path(PT_results_path, "input.files"), sep="\n")

  # 5. Run PAINTOR!
  PAINTOR.run(paintor_path = paintor_path, 
              PT_results_path = PT_results_path, 
              n_datasets = n_datasets, 
              n_causal = n_causal)

  # Summarise the results 
  paintor.results <- PAINTOR.process_results(PT_results_path=PT_results_path, 
                                             locus_name=locus_name)
  return(paintor.results)
}

