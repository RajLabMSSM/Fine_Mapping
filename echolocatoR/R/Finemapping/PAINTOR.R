# [GitHub](https://github.com/gkichaev/PAINTOR_V3.0)
# [LD Tutorial](https://github.com/gkichaev/PAINTOR_V3.0/wiki/2a.-Computing-1000-genomes-LD)

install_PAINTOR <- function(paintor_path="./echolocatoR/tools/PAINTOR_V3.0/"){
  cmd <- paste("cd", paintor_path,"&& bash install")
  system(cmd)
}

help_PAINTOR <- function(paintor_path="./echolocatoR/tools/PAINTOR_V3.0/"){
  system(file.path(paintor_path,"PAINTOR"))
}

zscore <- function(vec){
  z <- scale(vec, center = T, scale = T)
  return(z)
}


create_locusFile <- function(results_path,
                             PT_results_path,
                             GWAS_dataset_name="Nalls23andMe_2019",
                             eQTL_dataset_name="Fairfax_2014",
                             gene){
  printer("++ Creating Locus File...")
  # The locus file should at the very minimum contain the Z-scores for all the populations,
  # though metadata on each SNP such as chromosome, position, and rsid are recommended.
  # The top line of the locus file must contain a header with the names of the fields.
  # The Z-score of a SNP is the Wald statistic, obtained from standard regression of
  # the phenotype onto the SNP.

  # Import GWAS data 
  GWAS_path <- file.path(results_path,"Multi-finemap/Multi-finemap_results.txt")
  # Get the party started using the GWAS file
  merged_DT <- data.table::fread(GWAS_path, select = c("CHR","POS","SNP","Effect")) %>%
    dplyr::rename(CHR=CHR,
                  POS=POS,
                  RSID=SNP,
                  ZSCORE.P1=Effect) %>% unique() %>%
    dplyr::mutate(CHR=paste0("chr",CHR), 
                  ZSCORE.P1=zscore(ZSCORE.P1))

  if(!is.na(eQTL_dataset_name)){
    # Merge eQTL data together
    if(eQTL_dataset_name=="Fairfax_2014"){
      merged_eQTLs <- merge_eQTL_data(snp_list = unique(merged_DT$RSID),
                                      eQTL_SS_paths = file.path("Data/eQTL/Fairfax_2014",
                                                                c("CD14/LRRK2/LRRK2_Fairfax_CD14.txt",
                                                                  "IFN/LRRK2/LRRK2_Fairfax_IFN.txt",
                                                                  "LPS2/LRRK2/LRRK2_Fairfax_LPS2.txt",
                                                                  "LPS24/LRRK2/LRRK2_Fairfax_LPS24.txt")),
                                      expression_paths = file.path("Data/eQTL/Fairfax_2014",
                                                                   c("CD14/CD14.47231.414.b.txt",
                                                                     "IFN/IFN.47231.367.b.txt",
                                                                     "LPS2/LPS2.47231.261.b.txt",
                                                                     "LPS24/LPS24.47231.322.b.txt")),
                                      ## IMPORTANT: Use the genotype subset that includes all the LRRK2 variants
                                      ### Gathered from minerva using grep
                                      genotype_path = "Data/eQTL/Fairfax_2014/LRRK2.geno.txt",
                                      subset_genotype_file = F,
                                      probe_path = "Data/eQTL/Fairfax_2014/gene.ILMN.map",
                                      .fam_path = "Data/eQTL/Fairfax_2014/volunteers_421.fam",
                                      gene = gene,
                                      save_merged = F)
    }
    # Iterate over eQTL condtions and merge with GWAS file
    eQTL_conditions <-  file.path("./Data/eQTL",eQTL_dataset_name) %>%
      list.dirs(recursive = F) %>%
      basename()
    eQTL_conditions <- eQTL_conditions[!eQTL_conditions=="merged_results"]
    
    for(i in 1:length(eQTL_conditions) ){
      condition <- eQTL_conditions[i]
      DT <- subset(merged_eQTLs, Condition==condition, c("SNP","Effect")) %>%
        unique() %>% dplyr::mutate(Effect = zscore(Effect))
      # DT <- import_table(file_paths[i])[,c("RSID","ZSCORE")]
      colnames(DT) <- c("RSID", paste0("ZSCORE.P",i+1))
      merged_DT <- data.table:::merge.data.table(merged_DT, DT,
                                                 by = "RSID", allow.cartesian = T)
    }
  }

  return(merged_DT)
}

datatype_handler <- function(GWAS_dataset_name=NA, 
                             eQTL_dataset_name=NA, 
                             gene){
  if(!is.na(GWAS_dataset_name) & !is.na(eQTL_dataset_name)){
    printer("++ GWAS and eQTL input data detected. Feeding both into PAINTOR...")
    results_path <- file.path("./Data/GWAS",GWAS_dataset_name, gene) 
  } else if(!is.na(GWAS_dataset_name) & is.na(eQTL_dataset_name)){
      printer("++ Only GWAS input data detected. Feeding into PAINTOR...")
      results_path <- file.path("./Data/GWAS",GWAS_dataset_name, gene)
    } else if(is.na(GWAS_dataset_name) & !is.na(eQTL_dataset_name)){
      printer("++ Only eQTL input data detected. Feeding into PAINTOR...")
      results_path <- file.path("./Data/eQTL",eQTL_dataset_name,"merged_results", gene)
    } else {
      stop("Neither GWAS nor eQTL data detected. Please enter at least one valid dataset.")
    } 
  PT_results_path <- file.path(results_path,"PAINTOR")
  dir.create(PT_results_path, showWarnings = F, recursive = T)
  printer("+++ Results will be stored in  ====> ", PT_results_path) 
  return(PT_results_path)
} 

locusName_handler <- function(locus=NA, gene, GWAS_dataset_name=NA, eQTL_dataset_name=NA){
  if(is.na(locus)){ 
    locus_name <- paste(gene,GWAS_dataset_name,eQTL_dataset_name, sep = ".")
    return(locus_name)
  } else {
    return(locus)
  }
}


# ------------------------------#
# ----------- PAINTOR ----------#
# ------------------------------#
PAINTOR <- function(paintor_path = "./echolocatoR/tools/PAINTOR_V3.0",
                    GWAS_dataset_name=NA,
                    eQTL_dataset_name=NA,
                    gene,
                    locus=NA,
                    n_causal=5,
                    XGR_dataset="FANTOM5_Enhancer",
                    Roadmap_chromatinMarks_search=NA,
                    chromatin_state="TssA",
                    no_annotations=F){
  # Note: All file formats are assumed to be single space delimited.
  ## Quick setup
  # chromatin_state=NA; paintor_path = "./echolocatoR/tools/PAINTOR_V3.0"; GWAS_dataset_name = "Nalls23andMe_2019";eQTL_dataset_name = "Fairfax_2014";locus=NA; gene="LRRK2"
  locus <- locusName_handler(locus, gene, GWAS_dataset_name, eQTL_dataset_name)
  
  PT_results_path <- datatype_handler(GWAS_dataset_name=GWAS_dataset_name, 
                                      eQTL_dataset_name=eQTL_dataset_name, 
                                      gene=gene)
  results_path <- dirname(PT_results_path)


  # 1. Locus File 
  locus_DT <- create_locusFile(results_path=results_path,
                               PT_results_path=PT_results_path,
                               GWAS_dataset_name=GWAS_dataset_name,
                               eQTL_dataset_name=eQTL_dataset_name,
                               gene=gene)
  # Remove NAs (PAINTOR doesn't tolerate missing values)
  locus_DT <- locus_DT[complete.cases(locus_DT),]
  ## Save
  .locus_file <- file.path(PT_results_path, locus)
  data.table::fwrite(locus_DT, .locus_file,
                     sep=" ", quote = F, na = NA)
  n_datasets <- sum(grepl(colnames(locus_DT), pattern = "ZSCORE"))



  # 2. LD Matrix File
  ### "VERY IMPORTANT!
  # Particular care must be taken when computing LD from a reference panel such the 1000 genomes.
  ## It is imperative that all the reference and alternate alleles for SNPs from which the Z-scores
  ## were computed match the reference and alternate alleles of the reference panel.
  ## The output of PAINTOR will not be correct if there are mismatches of this type in the data.
  ##Please see wiki section 2a for instructions on how to use the LD utility provided with the software."
  printer("++ Creating LD Matrix File...")
  load(file.path(results_path,"plink/LD_matrix.RData"))
  ## Make sure SNPs are in the same order as the Locus File
  .LD_file <- file.path(PT_results_path, paste0(locus,".ld1"))
  ld.mat <- LD_matrix[locus_DT$RSID, locus_DT$RSID] %>%
    data.table::as.data.table()
  ## Write
  data.table::fwrite(ld.mat,
                     .LD_file,
                     sep = " ", quote = F,
                     col.names = F, row.names = F)


  # 3. Annotation Matrix File
  printer("++ Creating Annotation Matrix File...")
  .annotations_file <- file.path(PT_results_path, paste0(locus,".annotations.txt"))
  if(no_annotations){
    data.table::fwrite(data.frame(Default_Priors = rep(1,nrow(locus_DT))),
                       .annotations_file,
                       sep=" ", quote = F)
  } else{
    ## Download GRs and convert to BED
    if(!is.na(XGR_dataset)){
      lib.name <- XGR_dataset  # FANTOM5_Enhancer
      GR.annotations <- XGR::xRDataLoader(lib.name)
      BED_paths <- GRs.to.BED(GR.annotations = GR.annotations,
                              output_path = file.path(PT_results_path,"annotation_files"),
                              sep=" ")
    }
    if(!is.na(Roadmap_chromatinMarks_search)){
      if(Roadmap_chromatinMarks_search=="all_cell_types"){Roadmap_chromatinMarks_search <- ""}
      # Version 1: Subset locally 
      RoadMap_ref <- GoShifter.search_ROADMAP(fuzzy_search = Roadmap_chromatinMarks_search)  #
      bed_names <- GoShifter.bed_names(RoadMap_ref)
      BED_paths <- GoShifter.get_roadmap_annotations(bed.list = bed_names, 
                                                     chromatin_state = chromatin_state, 
                                                     verbose = T) 
      try({lapply(BED_paths, gunzip)})
      BED_paths <- gsub(".gz","",BED_paths)
      # Version 2: Subset via tabix
      # BED_paths <- lapply(RoadMap_ref$EID, function(eid){
      #   bed  = Roadmap_tabix(results_path = results_path, 
      #                        chrom = unique(locus_DT$CHR),
      #                        min_pos = min(locus_DT$POS),
      #                        max_pos = max(locus_DT$POS),
      #                        eid = eid, 
      #                        convert_to_GRanges = F)
      #   # Exclude BED files with no overlap with the locus
      #   if(nrow(bed)>0){
      #     bed_name <- gsub(".bgz","",unique(bed$File))
      #     bed_path <- file.path(PT_results_path,"annotation_files",bed_name)
      #     data.table::fwrite(bed, 
      #                        bed_path, 
      #                        sep = "\t", col.names = F)
      #     return(bed_path)
      #   } else {
      #     NULL
      #   } 
      # }) 
      # BED_paths <- BED_paths[!lapply(BED_paths, is.null) %>% unlist()] 

    }
    ## Use PAINTOR utility to format merge BED files 
    gz.files <- list.files("./echolocatoR/tools/Annotations", pattern = ".gz", full.names = T)
    for (g in gz.files){
      try({gunzip(g, overwrite=T)})
    } 
    BED_paths <- list.files("./echolocatoR/tools/Annotations", pattern="mnemonics.bed", full.names = T)
    
    data.table::fwrite(list( BED_paths ),
                       file.path(PT_results_path,"annotation_paths.txt"),
                       sep="\n")
    cmd <-  paste("python",file.path(paintor_path,"PAINTOR_Utilities","AnnotateLocus.py"),
                  "--input", file.path(PT_results_path,"annotation_paths.txt"),
                  "--locus", .locus_file,
                  "--out", .annotations_file,
                  "--chr","CHR",
                  "--pos","POS")
    system(cmd)
    colSums( data.table::fread(.annotations_file))
  }


  # 4. Input File
  data.table::fwrite(list(locus), file.path(PT_results_path,"input.files"), sep="\n")


  # "In order to determine which annotations are relevant to the phenotype being considered,
  ## we recommend running PAINTOR on each annotation independently."
  # for(bed in BED_paths){
    # bed <- BED_paths[1:2]

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

  results <- process_PAINTOR_results(PT_results_path=PT_results_path, locus=locus)
  return(results)
}

survey_annotation <- function(PT_results_path, 
                              locus="Locus1"){
  anno <- data.table::fread(file.path(PT_results_path, paste0(locus,".annotations.txt")))
  column_sums <- sort(colSums(anno), decreasing = T)
  signals <- column_sums[column_sums > 0]
  print(signals)
}

process_PAINTOR_results <- function(PT_results_path, 
                                    locus="Locus1"){
  results <- data.table::fread(file.path(PT_results_path, paste0(locus,".results.txt"))) %>% 
    arrange(desc(Posterior_Prob))
  # subset(results, Posterior_Prob>0) 
  return(results)
}

