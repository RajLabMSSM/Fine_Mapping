


TABIX.convert_file <- function(fullSS_path="./Data/GWAS/Nalls23andMe_2019/nallsEtAl2019_allSamples_allVariants.mod.txt",
                               chrom_col="CHR",
                               position_col="POS"){ 
  # Instantaneous header collection without reading in large file!
  # headers <- colnames(data.table::fread(cmd=paste("head -1",fullSS_path)))
  header.path <- file.path(dirname(fullSS_path),"header.txt")
  system(paste("head -1",fullSS_path,">",header.path))
  cDict <-  column_dictionary(file_path = fullSS_path)
  
  print("TABIX:: Converting full summary stats file to tabix format for fast querying...")
  fullSS.gz <- ifelse(endsWith(fullSS_path,".gz"), fullSS_path, paste0(fullSS_path,".gz"))
  gz_cat <- ifelse(endsWith(fullSS_path,".gz"),"gunzip -c","cat")
  cmd <- paste(gz_cat,
               fullSS_path,
               "| tail -n +2", # Get rid of header
               "| sort -k1,1 -k2,2n", #
               # ">",fullSS_path)
               "| bgzip >",fullSS.gz)
  print(cmd)
  system(cmd)
  # Rsamtools:::bgzip(file = fullSS_path, 
  #                   dest = fullSS.gz,
  #                   overwrite = T)  
  # idx <- Rsamtools::indexTabix(file = fullSS.gz,
  #                              format = "bed",
  #                              seq = 1,
  #                              start = 2,
  #                              end = 2,
  #                              skip=1)
  # Rsamtools::headerTabix(tbx)
  # tab <- TabixFile(zipped, idx)
  
   
  # Index
  TABIX.index <- function(fullSS.gz,
                          chrom_i=1,
                          pos_i=2){
    printer("TABIX:: Indexing ")
    cmd2 <- paste("tabix", 
                  "-f",
                  "-S 1",
                  "-s",chrom_i,
                  "-b",pos_i,
                  "-e",pos_i,
                  fullSS.gz)
    print(cmd2)
    system(cmd2)
  }  
  TABIX.index(fullSS.gz=fullSS.gz, 
              chrom_i=cDict[[chrom_col]], 
              pos_i=cDict[[position_col]])  
  return(fullSS.gz)
}


# Query
TABIX.query <- function(fullSS.gz, 
                        chr, 
                        start_pos, 
                        end_pos,
                        subset_path){
  coords <- paste0(chr,":",start_pos,"-",end_pos)
  # cmd4 <- paste("tabix -h",fullSS.gz,coords,">",subset_path)
  printer("TABIX:: Extracting subset of sum stats ==>",subset_path) 
  dat <- data.table::fread(cmd=paste("tabix -h",fullSS.gz,coords))
  printer("++ Returning",paste(dim(dat),collapse=" x "),"data.table")
  return(dat)
}


TABIX <- function(fullSS_path, 
                  subset_path,
                  is_tabix=F,
                  chrom_col="CHR", 
                  position_col="POS",
                  min_POS=NA,
                  max_POS=NA, 
                  chrom=NULL){  
  # Check if it's already an indexed tabix file 
  fullSS.gz <- ifelse(endsWith(fullSS_path,".gz"), fullSS_path, paste0(fullSS_path,".gz"))
  is_tabix <- file.exists(fullSS.gz) & file.exists(paste0(fullSS.gz,".tbi"))
  if(!is_tabix){ 
    fullSS.gz <- TABIX.convert_file(fullSS_path)
  } else { printer("TABIX:: Existing indexed tabix file detected") } 
  # Query
  header.path <- file.path(dirname(fullSS_path),"header.txt")
  cDict <- column_dictionary(file_path = header.path) 
  dat <- TABIX.query(fullSS.gz,
                      chr=chrom, 
                      start_pos=min_POS, 
                      end_pos=max_POS,
                      subset_path=subset_path) 
  colnames(dat) <- colnames(data.table::fread(header.path))
  printer("++ Saving query ==>", subset_path)
  head(dat)
  data.table::fwrite(dat, file = subset_path, nThread = 4, sep="\t")
}




