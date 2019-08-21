
# ********************************
# ************ fGWAS# ************ 
# ********************************
# """
# The R version of fGWAS appears to be extremely limited. 
# https://github.com/wzhy2000/fGWAS
# 
# Using the command line version instead.
# https://github.com/joepickrell/fgwas
## User Guide:
## https://github.com/joepickrell/fgwas/blob/master/man/fgwas_manual.pdf

# Pre-formatted annotations available here:
# https://github.com/joepickrell/1000-genomes

# """

install_fGWAS <- function(echo_path="./echolocatoR/tools", 
                          download.annotations=F){
  # Install R wrapper
  # devtools::install_github("wzhy2000/fGWAS/pkg")
  # # Install and compile tool
  # system(paste0("git submodule add  https://github.com/wzhy2000/fGWAS.git ", echo_path))
  # system( paste0("cd ",file.path(echo_path,"fGWAS")," & R CMD INSTALL pkg") )w
  tar.gz <- "0.3.6.tar.gz"
  tar <- file.path(echo_path,gsub(".gz","",tar.gz))
  # Download from github
  system( paste0("wget ",
                file.path("https://github.com/joepickrell/fgwas/archive/",tar.gz),
                " ",echo_path) )
  # Decompress
  system(paste0("tar -xf ",tar))
  # Compile
  fgwas.folder <- file.path(echo_path,paste0("fgwas-",gsub(".tar.gz","",tar.gz)) )
  cmd <- paste0("cd ",fgwas.folder," & ./configure & make")
  system(cmd)
  # Create symlink?...
  if(download.annotations){
    fgwas.download_annotations()
  }
}


fgwas.download_annotations <- function(FM_all,
                                       results_path="./Data/GWAS/Nalls23andMe_2019/_genome_wide",
                                       force_new_annot = F,
                                       dataset = ""
                                       ){
  # Path/file names
  Input_path <- file.path(results_path,"fGWAS","Input")
  dir.create(Input_path, showWarnings = F, recursive = T)
  annot.file.name <- file.path(Input_path, paste0("annotations.",dataset,".txt"))
  # Takes a long time: Redo only if it doesn't exist (or its being forced to)
  if(!file.exists(annot.file.name) | force_new_annot==T){
    # https://github.com/joepickrell/1000-genomes
    # system(paste("git submodule add https://github.com/joepickrell/1000-genomes.git",annot.folder))  
    base.url <- "https://github.com/joepickrell/1000-genomes/raw/master"
    FM.snps <- unique(FM_all$SNP) #paste(unique(FM_all$SNP), collapse="|")
    # Download each annotation and merge with data
    annots.filt <- lapply(unique(FM_all$CHR), function(chr){
      printer("fGWAS: Downloading annotations for Chrom",chr)
      file.url <- file.path(base.url, paste0("chr",chr,".annot.wdist.wcoding.gz"))
      chr.dat <- data.table::fread(file.url)
      CHR.dat <- subset(chr.dat, rs %in% FM.snps) 
      return(CHR.dat)
    }) %>% data.table::rbindlist()
    # Save filtered annotations 
    data.table::fwrite(annots.filt, annot.file.name, sep="\t")
  } else {
    printer("+ fGWAS: Existing annotations found.")
    printer("    Importing:",annot.file.name)
    annots.filt <- data.table::fread(annot.file.name)
    } 
  # Merge annotations afterwards
  FM_annot <- data.table:::merge.data.table(FM_all, annots.filt, 
                                            by.x = "SNP",
                                            by.y = "rs", 
                                            all.x = T)
  # Fill NAs with 0s in certain columns
  for (j in names(FM_annot)[24:ncol(FM_annot)]){
    set(FM_annot,which(is.na(FM_annot[[j]])),j,0)
  }  
  return(FM_annot)
}



fgwas.annotation_names <- function(fgwas="./echolocatoR/tools/fgwas-0.3.6/src/fgwas"){ 
  ## Get annotation names from summary file
  annot_files <- data.table::fread(file.path( dirname(dirname(fgwas)),
                                              "annot", "annotation_list_winfo.txt"),
                                   col.names = c('bed.name',"type","description"))
  annot_files$Name <- gsub(".bed.gz","",basename(annot_files$bed.name)) 
  return(annot_files)
}

fgwas.top_annotations <- function(dat.fgwas, annot_files, SNP.Group=""){
  ## Count how many SNPs have hits per annotation
  overlapping.snps <- subset(dat.fgwas, select=as.character(annot_files$Name) ) %>% 
    colSums() 
  printer("fGWAS: Annotations with the most overlapping SNPs:",SNP.Group)
  print(head(overlapping.snps %>% sort(decreasing = T), 5))
  overlap.df <- data.frame(overlapping.snps) %>% `colnames<-`(SNP.Group) 
  # data.table::data.table(overlap.df, keep.rownames = "Annot", key="Annot")
  return(overlap.df)
}
 


fgwas.prepare_input <- function(results_path="./Data/GWAS/Nalls23andMe_2019/_genome_wide",
                                SNP.Groups = c("leadSNP","CS","Consensus"),
                                force_new_annot = F,
                                dataset = "./Data/GWAS/Nalls23andMe_2019" 
                                ){
  # """
  # 1. SNPID: a SNP identifier
  # 2. CHR: the chromosome for the SNP
  # 3. POS: the position of the SNP
  # 4. F: the allele frequency of one of the alleles of the SNP
  # 5. Z: the Z-score for from the association study
  # 6. N: the sample size used in the association study at this SNP (this can vary from SNP to SNP due to, e.g., missing data).
  # """ 
  dir.create(file.path(results_path, "fGWAS"), showWarnings = F, recursive = F)
  
  FM_all <- merge_finemapping_results(minimum_support = 1, 
                                      include_leadSNPs = T) 
  FM_all <- subset(FM_all, Dataset==dataset)
  # Gather annotations data and merge with FM_all
  FM_annot <- fgwas.download_annotations(FM_all,
                                         force_new_annot = force_new_annot, 
                                         dataset = dataset)
  # List annotation names
  annot_files <- fgwas.annotation_names()
  
  
  # Construct input for each SNP.Group
  fgwas.inputs <- data.frame()
  overlap.df.list <- list() 
  for(group in SNP.Groups){
    printer("+ fGWAS: Creating", group, "file")
    # Subset according to which group we're looking at
    if(group=="leadSNP"){
      dat.fgwas <- subset(FM_annot, leadSNP==T)
    }
    if(group=="CS"){
      dat.fgwas <- subset(FM_annot, Support>0)
    }
    if(group=="Consensus"){
      dat.fgwas <- subset(FM_annot, Consensus_SNP==T)
    }
    # Construct data in fGWAS format
    dat.fgwas <- dat.fgwas  %>%
      dplyr::mutate(N=N_cases+N_controls) %>%
      dplyr::rename(SNPID=SNP, CHR=CHR, POS=POS, "F"=Freq, Z=Effect, 
                    N=N, NCASE=N_cases, NCONTROL=N_controls, Gene=Gene) %>% 
      arrange(POS)  
    # Assign each locus a segment ID
    seg.table <- data.table::data.table(Gene = unique(FM_annot$Gene), 
                                        SEGNUMBER = 1:length(unique(FM_annot$Gene)) )
    dat.fgwas <- data.table:::merge.data.table(dat.fgwas, 
                                               seg.table, 
                                               by = "Gene")
    dat.fgwas <- unique(dat.fgwas) %>% arrange(SEGNUMBER, CHR, POS)
    # Save file
    dat.path <- file.path(results_path,"fGWAS","Input",paste0("dat.fgwas.",group,".txt"))
    dir.create(dirname(dat.path), showWarnings = F, recursive = T) 
    data.table::fwrite(dat.fgwas, dat.path, sep = " ")
    gzip(dat.path, overwrite=T) 
    # Add to summary dataframe
    fgwas.inputs <- rbind(fgwas.inputs, data.frame(SNP.Group=group,
                                                   File=paste0(dat.path,".gz"),
                                                   N.SNPs=nrow(dat.fgwas)))
   
    # Report annotations with the most overlap
    overlap.df <- fgwas.top_annotations(dat.fgwas, annot_files, SNP.Group=group)
    overlap.df.list <- append(overlap.df.list, overlap.df)
  } 
  overlap.DF <- cbind.data.frame(overlap.df.list) %>% `row.names<-`(rownames(overlap.df)) 
  return(list(fgwas.inputs=fgwas.inputs, overlap.DF=overlap.DF))
}

 

fgwas.gather_results <- function(results_path, 
                                 output_dir=file.path(results_path,"fGWAS/Output")){
  printer("+ fGWAS: Gathering all results...")
  #### Default output:
  # 1. fgwas.params. The maximum likelihood parameter estimates for each parameter in the model. The columns are the name of the parameter (“pi” is the parameter for the prior probability that any given genomic region contains an association), the lower bound of the 95% confidence interval on the parameter, the maximum likelihood estimate of the parameter, and the upper bound of the 95% confidence interval on the parameter.
  # 2. fgwas.llk. The likelihood of the model. The lines are the log-likelihood of the model under the maximum likelihood parameter estimates, the number of parameters, and the AIC.
  
#### With "-print" flag on:
  # 1.fgwas.ridgeparams. The estimates of the parameters under the penalized likelihood. The first line of this file is the penalty used, then the penalized parameters estimates follow. If you have used the -xv flag, the last line will be the cross-validation likelihood.
  # 2. fgwas.segbfs.gz. The association statistics in each region of the genome defined in the model. The columns of this file are the block number, chromosome, start position, end posi- tion, maximum absolute value of Z-score, log regional Bayes factor, regional prior probability of association, log regional posterior odds for association, and the regional posterior proba- bility of association. The annotations of the region (if any) are in the remaining columns.
  # 3. fgwas.bfs.gz. The association statistics for each SNP in the genome as estimated by the model. The columns of this file are the SNP ID, the chromosome, genomic position, log Bayes factor, Z-score, estimated variance in the effect size, prior probability of association, two columns (pseudologPO and pseudoPPA) for internal use only, the posterior probability of association (conditional on there being an association in the region), and the region number. The annotations in the model (if any) then follow.
  
  # Gather file names
  llk.files <- list.files(output_dir, pattern = ".llk")
  params.files <-  list.files(output_dir, pattern = ".params")
  # Set up preliminary df
  df <- data.frame(llk.file=llk.files, 
                   params.file=params.files)
  df <- tidyr::separate(df, col = "llk.file", into=c("SNP.Group","Annot",NA), 
                        sep="__|\\.",remove=F)
  # Gather results stats for each file
  res_df <- lapply(df$llk.file, function(llk){
    # llk file
    res <- data.table::fread(file.path(output_dir,llk)) %>% data.table::transpose(fill=0)
    colnames(res) <- as.character(res[1,])
    res <- res[-1,]
    res$llk.file <- llk
    # params file
    param <- gsub(".llk",".params",llk)
    res.p <- data.table::fread(file.path(output_dir, param))
    res$parameter <- res.p$parameter
    res$CI_lo <- res.p$CI_lo
    res$CI_hi <- res.p$CI_hi
    res$estimate <- res.p$estimate  
    return(res)
  }) %>% data.table::rbindlist()
  # Merge stats with preliminary df
  DF <- data.table:::merge.data.table(df, res_df, by = "llk.file")
  printer("+ Data for",nrow(DF),"results gathered.")
  return(DF)
}



fgwas.run <- function(results_path="./Data/GWAS/Nalls23andMe_2019/_genome_wide",
                      fgwas.out = file.path(results_path,"fGWAS","Output"),
                      fgwas = "./echolocatoR/tools/fgwas-0.3.6/src/fgwas",
                      remove_tmps = T){
  # Prepare inputs
  prepare_input.list <- fgwas.prepare_input(results_path,
                                            SNP.Groups = c("leadSNP","CS","Consensus")) 
  fgwas.inputs <- prepare_input.list$fgwas.inputs
  overlap.DF <- prepare_input.list$overlap.DF  
  
  # Iterate over each annotation file,
  ## so you can figure out where are significant.
  dir.create(file.path(fgwas.out), showWarnings = F, recursive = T) 
  for(annot in annot_files$Name){
    printer("")
    printer("------------------")
    printer("------------------")
    printer("+ fGWAS: Using annotation =",annot)
    # Iterate over each SNP.group file, 
    ## so you can compare their enrichment results to one another. 
    for(row in 1:nrow(fgwas.inputs)){ 
      # Gather info
      group <- fgwas.inputs[row,"SNP.Group"]
      f <- fgwas.inputs[row,"File"]
      N.SNPs <- fgwas.inputs[row,"N.SNPs"]
      message("++ fGWAS: Running enrichment on:  ",group)
      # Construct command
      cmd <- paste0(#"cd ",fgwas.out," & ",
        fgwas,
        " -i ", f,
        # " -cc ",
        " -o ",file.path(fgwas.out, paste0(group,"__", annot)), 
        " -fine",
        " -w ",annot,
        # " -print" # Produces extra files
        " -k ",N.SNPs
      ) 
      system(cmd) 
    } 
  }
  # Gather results
  results.DF <- fgwas.gather_results(results_path = results_path) 
  data.table::fwrite(results.DF, file.path(results_path,"fGWAS","fGWAS.summary.txt"), sep="\t")
  # Delete tmp files (input/output)
  if(remove_tmps){
    printer("fGWAS: Removing tmp input files...")
    # Remove Input  
    suppressWarnings(file.remove(list.files(file.path(dirname(fgwas.out),"Input"), pattern = ".gz")))
    # Remove Output
    removeDirectory(fgwas.out, recursive = T, mustExist = F) 
  }
  return(results.DF)
}


fgwas.plot <- function(results.DF){
  ggplot(results.DF, aes(x=SNP.Group, y=`AIC:`)) + geom_col()
  ggplot(results.DF, aes(x=SNP.Group, y=estimate)) + geom_col()
}



# dat.path <- file.path(results_path,"fGWAS/fgwas.data.txt")
# data.table::fwrite(dat.fgwas, dat.path, sep="\t", 
#                    row.names = F, col.names = T)
# 
# r <- fGWAS::fg.load.simple(file.simple.snp = dat.path )
# 
# 



