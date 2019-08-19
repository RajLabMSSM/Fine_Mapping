
# ********************************
# ************ fGWAS# ************ 
# ********************************
# """
# The R version of fGWAS appears to be extremely limited. 
# https://github.com/wzhy2000/fGWAS
# 
# Using the command line version instead.
# https://github.com/joepickrell/fgwas

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
  # Create symlink
  # fgwas <- file.path(fgwas.folder,"src/fgwas") # ./echolocatoR/tools/fgwas-0.3.6/src/fgwas
  # system(fgwas)

  if(download.annotations){
    fgwas.download_annotations()
  }
}


fgwas.download_annotations <- function(fgwas = "./echolocatoR/tools/fgwas-0.3.6/src/fgwas"){
  annot.folder <- file.path(fgwas,"annot")
  dir.create(annot.folder, showWarnings = F, recursive = T)
  # https://github.com/joepickrell/1000-genomes
  system(paste("git submodule add https://github.com/joepickrell/1000-genomes.git",annot.folder))  
  # return(annot.folder)
}


 
fgwas.prepare_input <- function(results_path,
                                SNP.Groups = c("leadSNP","CS","Consensus") 
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
  fgwas.inputs <- c()
  for(group in SNP.Groups){
    printer("+ fGWAS: Creating", SNP.Group, "file")
    # Subset according to which group we're looking at
    if(group=="leadSNP"){
      dat.fgwas <- subset(FM_all, leadSNP==T)
    }
    if(group=="CS"){
      dat.fgwas <- subset(FM_all, Support>0)
    }
    if(group=="Consensus"){
      dat.fgwas <- subset(FM_all, Consensus_SNP==T)
    }
    # Construct data in fGWAS format
    dat.fgwas <- dat.fgwas  %>%
      dplyr::mutate(N=N_cases+N_controls) %>%
      dplyr::select(SNPID=SNP, CHR, POS, "F"=Freq, Z=Effect, 
                    N, NCASE=N_cases, NCONTROL=N_controls, Gene) %>% 
      arrange(POS)  
    # Assign each locus a segment ID
    seg.table <- data.table::data.table(Gene = unique(FM_all$Gene), 
                                        SEGNUMBER = 1:length(unique(FM_all$Gene)) )
    dat.fgwas <- data.table:::merge.data.table(dat.fgwas, 
                                               seg.table, 
                                               by = "Gene")
    dat.fgwas <- unique(dat.fgwas) %>% arrange(SEGNUMBER)
    # Write file
    
    dat.path <- file.path(results_path,"fGWAS","Input",paste0("dat.fgwas.",group,".txt"))
    dir.create(dirname(dat.path), showWarnings = F, recursive = T)
    fgwas.inputs <- append(fgwas.inputs,
                          paste0(dat.path,".gz"))
    data.table::fwrite(dat.fgwas, dat.path, sep = " ")
    gzip(dat.path, overwrite=T) 
  } 
  return(unique(fgwas.inputs))
}



fgwas.run <- function(results_path, 
                      fgwas.inputs,
                      fgwas.out = file.path(results_path,"fGWAS","Output"),
                      fgwas = "./echolocatoR/tools/fgwas-0.3.6/src/fgwas",
                      remove_tmps = T){ 
  dir.create(file.path(fgwas.out), showWarnings = F, recursive = T) 
  # Prepare inputs
  fgwas.inputs <- fgwas.prepare_input(results_path,
                                     SNP.Groups = c("leadSNP","CS","Consensus"))
  # Get annotation names
  annot_files <- data.table::fread(file.path( dirname(dirname(fgwas)),
                                              "annot", "annotation_list_winfo.txt"), 
                                   col.names = c('bed.name',"type","description"))
  # Show fGWAS help
  system(fgwas)
  # Iterate over each annotation file,
  ## so you can figure out where are significant.
  for(annot in annot_files$bed.name){
    print("")
    print("------------------")
    print("------------------")
    printer("+ fGWAS: Using annotation =",annot)
    # Iterate over each SNP.group file, 
    ## so you can compare their enrichment results to one another.
    for(f in fgwas.inputs){
      printer("++ fGWAS: Running enrichment on",f)
      group <- strsplit(basename(f),"\\.")[[1]][3]
      cmd <- paste0(#"cd ",fgwas.out," & ",
        fgwas,
        " -i ", f,
        " -cc ",
        " -o ",file.path(fgwas.out, paste0(group,"__", annot)), 
        " -fine ",
        " - w ",annot
        # " -k 10"
      ) 
      system(cmd) 
    } 
  }
  
  if(remove_tmps){
    printer("fGWAS: Removing tmp input files...")
    suppressWarnings(file.remove(fgwas.inputs))
  }
}





# dat.path <- file.path(results_path,"fGWAS/fgwas.data.txt")
# data.table::fwrite(dat.fgwas, dat.path, sep="\t", 
#                    row.names = F, col.names = T)
# 
# r <- fGWAS::fg.load.simple(file.simple.snp = dat.path )
# 
# 



