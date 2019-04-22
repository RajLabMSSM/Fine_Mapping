# %%%%%%%%%%%%%%%%% #
####### QUERY #######
# %%%%%%%%%%%%%%%%% #

## Import Sig GWAS/QTL Summary Statistics
# For each gene, get the position of top SNP (the one with the greatest effect size/Beta)
import_topSNPs <- function(file_path, caption="", sheet = 1,
                           chrom_col="CHR", position_col="POS", snp_col="SNP",
                           pval_col="P", effect_col="Effect", gene_col="Gene"){
  ## Only the significant subset of results
  if(endsWith(file_path, ".xlsx") | endsWith(file_path, ".xlsm")){
    top_SNPs <- read_excel(path = file_path, sheet = sheet)
  } else if (endsWith(file_path, ".csv")){
    top_SNPs <- fread(file=file_path, sep = ",", header = T, stringsAsFactors = F )
  }
  else if (endsWith(file_path, ".txt")){
    top_SNPs <- fread(file=file_path, sep = "\t", header = T, stringsAsFactors = F )
  } else {print("File type must be .xlsx, .cs, or tab-delimited .txt")}

  top_SNPs <- subset(top_SNPs, select=c(chrom_col, position_col, snp_col,
                                        pval_col, effect_col, gene_col ))
  # Standardize names
  colnames(top_SNPs) <- c("CHR","POS","SNP","P","Effect","Gene")
  # Get only the top SNP (sorted by lowest p-val, then highest Effect size) for each gene
  top_SNPs <- top_SNPs %>% arrange(P, desc(Effect)) %>% group_by(Gene) %>% slice(1)
  top_SNPs$CHR <- gsub("chr", "",top_SNPs$CHR)
  top_SNPs <- cbind(Coord= paste(top_SNPs$CHR, top_SNPs$POS, sep=":"),
                    top_SNPs)
  createDT(top_SNPs, caption)
  return(top_SNPs)
}
# top_SNPs <- import_topSNPs(
#   file_path = Data_dirs$Nalls_2019$topSS,
#   sheet="Data",
#   chrom_col = "CHR", position_col = "BP", snp_col="SNP",
#   pval_col="P, all studies", effect_col="Beta, all studies", gene_col="Nearest Gene",
#   caption= "Nalls et al. (2018) PD GWAS Summary Stats")


## Get Flanking SNPs
# Get all genes surrounding the index SNP (default is 500kb upstream + 500kb downstream)
# 1000000 bp
column_dictionary <- function(file_path){
  # Get the index of each column name
  f <- fread(file_path, nrows = 0, header = T)
  cNames <- colnames(f)
  colDict <- setNames(1:length(cNames), cNames  )
  return(colDict)
}




auto_topSNPs_sub <- function(top_SNPs, query, gene){
  # If no top_SNPs dataframe is supplied,
  ## this function will sort by p-value and then effect size,
  ## and use the SNP in the first row.
  if(toString(top_SNPs)=="auto"){
    top_SNPs <- query %>% mutate(Gene=gene) %>%
      arrange(P, desc(Effect)) %>% group_by(Gene) %>% slice(1)
  }
  topSNP_sub <- top_SNPs[top_SNPs$Gene==gene & !is.na(top_SNPs$Gene),][1,]
  return(topSNP_sub)
}

preprocess_subset <- function(gene, top_SNPs, subset_path, file_sep,
                              chrom_col="CHR", position_col="POS", snp_col="SNP",
                              pval_col="P", effect_col="Effect", stderr_col="StdErr",
                              tstat_col="t-stat", return_dt=T){
  cat("\n Preprocessing SNP subset... \n")
  query_check <- data.table::fread(subset_path, sep=file_sep, nrows = 2)

  if(dim(query_check)[1]==0){
    file.remove(subset_path)
    stop("\n Could not find any rows in full data that matched query :(")
  } else{
    query <- data.table::fread(subset_path, header=T, stringsAsFactors = F, sep = file_sep)
    if(stderr_col=="calculate"){
      cat("Calculating Standard Error...\n")
      query <- data.table::fread(subset_path, header=T, stringsAsFactors = F, sep = file_sep)
      # Calculate StdErr
      query$StdErr <- subset(query, select=effect_col) / subset(query, select=tstat_col)
      stderr_col <- "StdErr"
    }
    query <- query %>% subset(select=c(chrom_col,position_col, snp_col, pval_col, effect_col, stderr_col)) %>%
      dplyr::rename(CHR=chrom_col,POS=position_col, SNP=snp_col, P=pval_col, Effect=effect_col, StdErr=stderr_col) %>%
      mutate(Location=paste(CHR,":",POS,sep=""))

    topSNP_sub <- auto_topSNPs_sub(top_SNPs, query, gene)

    ## Remove SNPs with NAs in stats
    query[(query$P<=0)|(query$P>1),"P"] <- 1
    # Get just one SNP per location (just pick the first one)
    query <- query %>% group_by(Location) %>% slice(1)
    # Mark lead SNP
    query$leadSNP <- ifelse(query$SNP==topSNP_sub$SNP, T, F)
    # Only convert to numeric AFTER removing NAs (otherwise as.numeric will turn them into 0s)
    query <- query  %>%
      mutate(Effect=as.numeric(Effect), StdErr=as.numeric(StdErr), P=as.numeric(P))

    data.table::fwrite(query, subset_path, sep = "\t")
    if(return_dt==T){return(query)}
  }
}


query_by_coordinates <- function(top_SNPs, gene, subset_path, fullSS_path,
                                 file_sep, chrom_col, position_col,
                                 minPos, maxPos, bp_distance){
  topSNP_sub <- top_SNPs[top_SNPs$Gene==gene & !is.na(top_SNPs$Gene),]
  if(is.null(minPos)){minPos <- topSNP_sub$POS - bp_distance}
  if(is.null(maxPos)){maxPos <- topSNP_sub$POS + bp_distance}
  cat("---Min snp position:",minPos, "---\n")
  cat("---Max snp position:",maxPos, "---\n")
  colDict <- column_dictionary(fullSS_path)
  awk_cmd <- paste("awk -F '",file_sep,"' 'NR==1 {print $0} NR>1 { if($",colDict[[chrom_col]]," == ",topSNP_sub$CHR,
                   " && ($", colDict[[position_col]]," >= ",minPos," && $",colDict[[position_col]]," <= ",maxPos,")) {print $0} }' ",fullSS_path,
                   " > ",subset_path,sep="")
  cat("\n",awk_cmd,"\n")
  system(awk_cmd)
}

query_by_coordinates_merged <- function(top_SNPs, fullSS_path, subset_path, gene,
                                        chrom_col, file_sep=" ", location_sep=":",
                                        minPos, maxPos, bp_distance){
  topSNP_sub <- top_SNPs[top_SNPs$Gene==gene & !is.na(top_SNPs$Gene),][1,]
  if(is.null(minPos)){minPos <- topSNP_sub$POS - bp_distance}
  if(is.null(maxPos)){maxPos <- topSNP_sub$POS + bp_distance}
  cat("---Min snp position:",minPos, "---\n")
  cat("---Max snp position:",maxPos, "---\n")
  colDict <- column_dictionary(fullSS_path)
  awk_cmd <- paste("cat ",fullSS_path," | tr -s '",location_sep,"' '",file_sep,"'",
                   " | awk -F '",file_sep,"' 'NR==1 {print \"CHR POS \" $2\" \" $3\" \" $4\" \" $5\" \" $6 }",
                   " NR>1 {if($",colDict[[chrom_col]]," == ",topSNP_sub$CHR," && ",
                   "($",colDict[[chrom_col]]+1," >=",minPos,"&& $",colDict[[chrom_col]]+1," <=",maxPos,")) {print $0}}'",
                   " | tr -s '",location_sep,"' '\t' > ",subset_path, sep="")
  # awk_cmd <- paste("cat ",fullSS_path," | tr -s '",location_sep,"' '",file_sep,"'",
  #                  " | awk -F '",file_sep,"' 'NR==1 {print \"CHR POS \" $2\" \" $3\" \" $4\" \" $5\" \" $6 }",
  #                  " NR>1 {if($1 == 1 && ($2 >=",minPos,"&& $2 <=",maxPos,")) {print $0}}'",
  #                  " | tr -s '",location_sep,"' '\t' > ",subset_path, sep="")
  cat("\n",awk_cmd,"\n")
  system(awk_cmd)
}

query_by_gene <- function(fullSS_path, subset_path, gene, gene_col, file_sep){
  colDict <- column_dictionary(fullSS_path)
  awk_cmd <- paste("awk -F '",file_sep,"' 'NR==1{print $0} $",colDict[[gene_col]]," == \"",gene,"\" {print $0}' ",fullSS_path,
                   "| tr -s '",file_sep,"' '\t'  > ",subset_path, sep="")
  cat("\n",awk_cmd,"\n")
  system(awk_cmd)
}

query_by_probe <- function(fullSS_path, subset_path, gene, gene_col, file_sep, probe_path){
  probe_path="Data/eQTL/Fairfax/gene.ILMN.map"
  subset_path="Data/eQTL/Fairfax/LRRK2_Fairfax_CD14.txt"
  file_sep="\t"
  gene="LRRK2"
  probes <- find_probes(probe_path, gene)
  PROBES <- paste("\"",paste(probes[1], collapse="\" || \""),"\"",sep="")
   
  
  colDict <- column_dictionary(fullSS_path)
  awk_cmd <- paste("awk -F '",file_sep,"' 'NR==1 {print $0} NR>1 $",colDict[[gene_col]]," == \"",probes[1],"\" {print $0}}' ",fullSS_path,
                   " > ",subset_path, sep="")
  cat("\n",awk_cmd,"\n")
  system(awk_cmd)
  
  data.table::fread(fullSS_path, nrows = 5)
  
}



query_handler <- function(gene, fullSS_path, top_SNPs, subset_path, #top_SNPs="auto",
                          gene_col="Gene", minPos=NULL, maxPos=NULL, bp_distance=500000,
                          chrom_col="CHR", position_col="POS", file_sep="\t",
                          query_by="coordinates"){
  cat("\n ++ Query Method: '",query_by,"'\n", sep="")
  if(query_by=="coordinates"){
    query_by_coordinates(top_SNPs=top_SNPs, gene=gene,
                         subset_path=subset_path, fullSS_path=fullSS_path,
                         chrom_col=chrom_col, position_col=position_col,
                         minPos=minPos, maxPos=maxPos, bp_distance=bp_distance,
                         file_sep=file_sep)
  }
  if(query_by=="coordinates_merged"){
    query_by_coordinates_merged(top_SNPs=top_SNPs, fullSS_path=fullSS_path, subset_path=subset_path,
                                file_sep=file_sep,  gene=gene, location_sep= ":",
                                minPos=minPos, maxPos=maxPos, bp_distance=bp_distance,
                                chrom_col=chrom_col)
  }
  if(query_by=="gene"){
    query_by_gene(fullSS_path=fullSS_path, subset_path=subset_path,
                  gene=gene, gene_col=gene_col, file_sep=file_sep)
  }

}


check_if_empty <- function(file_path, file_sep="\t"){
  rowCheck <- dim(data.table::fread(file_path, nrows = 2, sep = file_sep))[1]
  if(rowCheck==0){
    stop("\n No SNPs identified within the summary stats file that met your criterion. :o \n")
  } else {cat("\n Subset file looks good! :)")}
}




extract_SNP_subset <- function(gene, fullSS_path, subset_path="auto", top_SNPs="auto", bp_distance=500000,
                               chrom_col="CHR", position_col="POS", snp_col="SNP", gene_col="Gene",
                               pval_col="P", effect_col="Effect", stderr_col="StdErr",
                               tstat_col="t-stat", superpopulation="",
                               minPos=NULL, maxPos=NULL, file_sep="\t",
                               query_by="coordinates"){
  # topSNP_sub <- top_SNPs[top_SNPs$Gene==gene & !is.na(top_SNPs$Gene),][1,]
  # if(is.null(minPos)){minPos <- topSNP_sub$POS - bp_distance}
  # if(is.null(maxPos)){maxPos <- topSNP_sub$POS + bp_distance}
  # cat("---Min snp position:",minPos, "---\n")
  # cat("---Max snp position:",maxPos, "---\n")
  if(file.exists(subset_path)){
    cat("\ Subset file already exists. Importing",subset_path,"...\n")
    check_if_empty(subset_path, file_sep=file_sep)
    query <- data.table::fread(subset_path, header=T, stringsAsFactors=F, sep="\t")
  } else {
    # Extract subset with awk
    cat("Extracting relevant variants from fullSS...\n")
    start_query <- Sys.time()
    # Function selects different methods of querying your SNPs
    query_handler(gene=gene, fullSS_path=fullSS_path, subset_path=subset_path, top_SNPs=top_SNPs,
                  gene_col=gene_col, chrom_col=chrom_col, position_col=position_col,
                  file_sep=file_sep, minPos=minPos, maxPos=maxPos, bp_distance=bp_distance,
                  query_by=query_by)
    # query_by_coordinates(topSNP_sub, subset_path,fullSS_path, file_sep, chrom_col, position_col, minPos, maxPos)
    # Clean file
    query <- preprocess_subset(gene=gene, top_SNPs=top_SNPs, subset_path=subset_path,
                               file_sep=file_sep, chrom_col=chrom_col,
                               position_col=position_col, snp_col=snp_col,pval_col=pval_col,
                               effect_col=effect_col, stderr_col=stderr_col,
                               tstat_col=tstat_col)
    end_query <- Sys.time()
    cat("\nExtraction completed in", round(end_query-start_query, 2),"seconds \n")
  }
  return(query)
}
# flankingSNPs <- extract_SNP_subset(gene="PTK2B", top_SNPs,
#                                   file_path = Data_dirs$MESA_CAU$topSS,
#                                   chrom_col = "chr",  pval_col = "pvalue", snp_col = "snps",
#                                   effect_col = "beta", position_col = "pos_snps",
#                                   stderr_col = "calculate")
#


find_probes <- function(map_file, genes){
  df  = data.table::fread(map_file)
  return(subset(df, GENE %in% genes)$PROBE_ID) 
}
 


# *********** PREPROCESSING ***********
# ## 23andME
# add_snpIDs_to_fullSS <-function(snpInfo_path, fullSS_path){
#   snp_ids <- fread(snpInfo_path,
#                    select=c("all.data.id", "assay.name", "scaffold", "position"),
#                    key="all.data.id")
#   snp_ids$scaffold <- gsub("chr", "",snp_ids$scaffold)
#   fullSS <- fread(fullSS_path,  key="all.data.id")
#   dim(fullSS)
#   SS_merged <- snp_ids[fullSS] %>% rename(SNP="assay.name", CHR="scaffold")
#   dim(SS_merged)
#   fwrite(SS_merged, "Data/Parkinsons/23andMe/PD_all_post30APRIL2015_5.2_extended.txt",
#          sep = "\t", quote = F, row.names = F)
#   remove(snp_ids, fullSS, SS_merged)
# }
# # add_snpIDs_to_fullSS(snpInfo_path = "Data/Parkinsons/23andMe/all_snp_info-5.2.txt",
# #                      fullSS_path = "Data/Parkinsons/23andMe/PD_all_post30APRIL2015_5.2.dat")



# Preprocessing
# add_snp_info <-function(snpInfo_path, fullSS_path, newSS_path){
#   # sqldf::read.csv.sql(file=snpInfo_path, header = T, sep = "\t",
#   #                     sql = paste("select * from file where location IN (", paste("'",snp_locs,"'", collapse=",",sep=""),")",sep="")
#   #                     )
#   cat("\nLoading SNP Info file...\n")
#   snp_ids <- fread(snpInfo_path,header = T, stringsAsFactors = F, key="location", colClasses = rep("character",2))
#
#   cat("\nLoading full SS file...\n")
#   fullSS <- fread(fullSS_path, stringsAsFactors = F,
#                   colClasses = c(rep("character",3), rep("numeric",7), "character", rep("numeric",4)))
#   fullSS$MarkerName <- gsub("chr","",fullSS$MarkerName)
#   fullSS <- fullSS %>% rename(location="MarkerName")
#   fullSS <- data.table(fullSS, key="location")
#   dim(fullSS)
#   cat("\nMerging files...\n")
#   SS_merged <- snp_ids[fullSS]
#   dim(SS_merged)
#   cat("\nWriting new file...\n")
#   fwrite(SS_merged, newSS_path, sep = "\t", quote = F, row.names = F)
#   remove(snp_ids, fullSS, SS_merged)
# }
# add_snp_info(snpInfo_path = "Data/HRC.RSID.txt",
#              fullSS_path = "Data/Parkinsons/Nalls_2019/nallsEtal2019_no23andMe.tab.txt",
#              newSS_path="Data/Parkinsons/Nalls_2019/nallsEtal2019_no23andMe_extended.txt")

#
# split_location_col <- function(input_path, output_path){
#   awk_cmd <- paste("awk 'BEGIN {print \"CHR\tPOS\tRSID\tA1\tA2\tfreq\tbeta\tse\tp\tN_cases\tN_controls\"}; FNR>1 {split($1, c, \":\"); print c[1] \"\t\" c[2] \"\t\" $2 \"\t\" $4 \"\t\" $5 \"\t\" $6 \"\t\" $7 \"\t\" $8 \"\t\" $9 \"\t\" $10 \"\t\" $11}' "
#          input_path," > ", output_path, sep="")
#   cat(awk_cmd)
#   system(awk_cmd)
# }



############ QTL ############
# Subset eQTL data to just the parts you need
subset_eQTL_SS <- function(fullSS_path, subset_path, gene,
                           gene_col="Gene", chrom_col="CHR", position_col="POS",
                           snp_col="SNP", pval_col="P", effect_col="Effect",
                           bp_distance=500000, file_sep="\t",
                           force_new_subset=T){
  start_eQTL_sub <- Sys.time()

  if(file.exists(subset_path) & force_new_subset==F){
    cat("Subset file already exists. Importing",subset_path,"...\n")
    check_if_empty(subset_path)
    # Extract lead SNPs per gene from subset data, using the subset SS as a reference file
    top_SNPs <- data.table::fread(subset_path) %>% mutate(Gene=gene) %>%
      arrange(P, desc(Effect)) %>% group_by(Gene) %>% slice(1)
  } else {
    colDict <- column_dictionary(fullSS_path)
    # Extract subset with awk
    cat("Extracting relevant variants from fullSS...\n")
    if(chrom_col==position_col){
      # Automatically detect when the chrom and position cols are merged into one and parse accordingly
      query_by_location(fullSS_path = fullSS_path, file_sep = file_sep, location_sep = ":", chr = top_SNPs$CHR,
                        minPos = (top_SNPs$POS - bp_distance), maxPos =  (top_SNPs$POS + bp_distance))
    } else {
      # Otherwise, just subset by gene
      query_by_gene(fullSS_path = fullSS_path, subset_path = subset_path,
                    gene = gene, gene_col = gene_col)
    }
    check_if_empty(subset_path)
    # Extract lead SNPs per gene from subset data, using the subset SS as a reference file
    top_SNPs <- import_topSNPs(
      file_path = subset_path,
      chrom_col = chrom_col, position_col= position_col, snp_col=snp_col,
      pval_col=pval_col, effect_col=effect_col, gene_col=gene_col,
      caption= paste(superpop,": eQTL Summary Stats"))
    # Make format/col names consistent
    preprocess_subset(topSNP_sub = top_SNPs, subset_path = subset_path,
                      chrom_col = chrom_col, position_col = position_col, snp_col = snp_col,
                      pval_col = pval_col, effect_col = effect_col, stderr_col = stderr_col, file_sep = "\t",
                      return_dt=F)
  }
  end_eQTL_sub <- Sys.time()
  cat("\nExtraction completed in", round(end_eQTL_sub-start_eQTL_sub, 2),"seconds \n")
  return(top_SNPs)
}



translate_population <- function(superpopulation){
  pop_dict <- list("AFA"="AFR", "CAU"="EUR", "HIS"="AMR",
                   "AFR"="AFR","EUR"="EUR", "AMR"="AMR")
  return(pop_dict[superpopulation][[1]])
}






# import info from FUMA instead

import_FUMA <- function(topSS_path, geneList, subset_path){
  # topSS_path = Data_dirs$Kunkle_2019$topSS
  # risk_loci <- fread(file.path(dirname(topSS_path), "FUMA/GenomicRiskLoci.txt"))
  # mapped_genes <- fread(file.path(dirname(topSS_path), "FUMA/magma.genes.out"))
  annovar <- data.table::fread(file.path(dirname(topSS_path), "FUMA/annov.txt"))
  annovar_sub <- subset(annovar, symbol %in% geneList)

  mapped_genes_sub <- subset(mapped_genes, SYMBOL %in% geneList)

  # Get gene range directly from biomart, then query for snps in that range
  geneSNPs <- function(gene, file_path,
                       chrom_col="CHR", position_col="POS", snp_col="SNP",
                       pval_col="P", effect_col="Effect", stderr_col="StdErr", sep="\t"){
    # bm <- biomart_geneInfo(gene)
    gene_start <- min(annovar_sub$pos)
    gene_end <- max(annovar_sub$pos)
    gene_chr <- unique(annovar_sub$chr)[1]

    colDict <- column_dictionary(file_path)
    superpop <- translate_population(superpopulation)
    dataset_name <- get_dataset_name(file_path)
    subset_path <- paste(dirname(file_path),"/",gene,"_",dataset_name,"_subset.txt",sep="")
    awk_cmd <- paste("awk -F \"",sep,"\" 'NR==1{print $0}{ if(($",colDict[chrom_col]," == ",gene_chr,")",
                     " && ($",colDict[position_col]," >= ",gene_start," && $",colDict[position_col]," <= ",gene_end,")) { print } }' ",file_path,
                     " > ",subset_path,sep="")
    system(awk_cmd)
  }
}




# BIOMART
biomart_snps_to_geneInfo <- function(snp_list){
  # listMarts()
  snp_mart = useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp")
  # View(listFilters(snp_mart))
  # View(listAttributes(snp_mart))
  snp_results <- biomaRt::getBM(snp_mart, filters="snp_filter", values=snp_list,
                                attributes=c("refsnp_id","snp","chr_name", "chrom_start","chrom_end",
                                             "associated_gene","ensembl_gene_stable_id" ) )
  # # Split ensembl IDs
  gene_mart = useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
  gene_results <- biomaRt::getBM(mart = gene_mart, filters = "ensembl_gene_id",
                                 # values = unlist(strsplit(snp_results$ensembl, ";")),
                                 values = snp_results$ensembl_gene_stable_id,
                                 attributes = c("hgnc_symbol","external_gene_name","ensembl_gene_id",
                                                "chromosome_name", "start_position", "end_position") )
  snp_results <-snp_results %>%
    mutate(ensembl = strsplit(as.character(ensembl_gene_stable_id), ";")) %>%
    tidyr::unnest(ensembl)
  merged_df <- data.table(gene_results, key = "ensembl_gene_id")[data.table(snp_results, key = "ensembl")]
  return(merged_df)
}
# biomart_snps_to_geneInfo(c("rs114360492"))

biomart_geneInfo <- function(geneList){
  # listDatasets(useMart("ENSEMBL_MART_ENSEMBL") )
  gene_mart = biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
  # View(listFilters(gene_mart))
  # View(listAttributes(gene_mart))
  gene_results <- biomaRt::getBM(mart = gene_mart, filters = "hgnc_symbol",
                                 # values = unlist(strsplit(snp_results$ensembl, ";")),
                                 values = geneList,
                                 attributes = c("hgnc_symbol","external_gene_name","ensembl_gene_id",
                                                "chromosome_name", "start_position", "end_position") )
  return(gene_results)
}
# biomart_geneInfo(c("PTK2B","CLU","APOE"))

