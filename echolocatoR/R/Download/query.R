# %%%%%%%%%%%%%%%%% #
####### QUERY #######
# %%%%%%%%%%%%%%%%% #

## Import Sig GWAS/QTL Summary Statistics
# For each gene, get the position of top SNP (the one with the greatest effect size/Beta)
import_topSNPs <- function(topSS_path, 
                           caption="", 
                           sheet = 1,
                           chrom_col="CHR", 
                           position_col="POS", 
                           snp_col="SNP",
                           pval_col="P", 
                           effect_col="Effect", 
                           gene_col="Gene",
                           group_by_locus=F,
                           locus_col="Locus",
                           remove_variants=F
                           ){
  # Import top SNPs
  topSNPs_reader <-function(topSS_path, sheet = 1){
    if(endsWith(topSS_path, ".xlsx") | endsWith(topSS_path, ".xlsm")){
      top_SNPs <- readxl::read_excel(path = topSS_path, sheet = sheet) %>% data.table::data.table()
    } else if (endsWith(topSS_path, ".csv")){
      top_SNPs <- data.table::fread(file=topSS_path, sep = ",", header = T, stringsAsFactors = F )
    }
    else if (endsWith(topSS_path, ".txt")){
      top_SNPs <-  data.table::fread(file=topSS_path, sep = "\t", header = T, stringsAsFactors = F )
    } else {printer("File type must be .xlsx, .cs, or tab-delimited .txt")}
    return(top_SNPs)
  }
  top_SNPs <- topSNPs_reader(topSS_path, sheet)
  orig_top_SNPs <- top_SNPs

  # Standardize col names 
  top_SNPs <- top_SNPs %>%
    dplyr::select(Gene=gene_col,
                  CHR=chrom_col, 
                  POS=position_col, 
                  SNP=snp_col,
                  P=pval_col, 
                  Effect=effect_col) 
    # Add Locus column
    if(locus_col %in% colnames(orig_top_SNPs)){
      LOCUS_vector <- dplyr::select(orig_top_SNPs, Locus=locus_col)
      top_SNPs <- cbind(LOCUS_vector, top_SNPs)
    } else {
      top_SNPs <- cbind(Locus=top_SNPs$Gene, top_SNPs)
    }
  
    # Remove specific variants
    if(remove_variants != F){
      top_SNPs <- subset(top_SNPs, !(SNP %in% remove_variants))
    }
    
    # Get only the top SNP (sorted by lowest p-val, then highest Effect size) for each gene
    top_SNPs <- top_SNPs %>% 
      arrange(P, desc(Effect)) %>% 
      group_by(Gene) %>% 
      dplyr::slice(1)
    top_SNPs$CHR <- gsub("chr", "",top_SNPs$CHR) 
    top_SNPs$CHR <- as.numeric(top_SNPs$CHR) 
    # Get the top representative SNP and Gene per locus (by lowest p-value)
    if(group_by_locus){
      top_SNPs <- top_SNPs %>%
        arrange(P) %>% 
        dplyr::group_by(Locus) %>% dplyr::slice(1) %>% 
        replace(., .=="NA", NA) %>% 
        subset(!is.na(Locus)) 
    }
    
    createDT(top_SNPs, caption)
    return(data.table::data.table(top_SNPs))  
}
  
column_dictionary <- function(file_path){
  # Get the index of each column name
  f <- data.table::fread(file_path, nrows = 0, header = T)
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
      arrange(P, desc(Effect)) %>% group_by(Gene) %>% dplyr::slice(1)
  }
  topSNP_sub <- top_SNPs[top_SNPs$Gene==gene & !is.na(top_SNPs$Gene),][1,]
  return(topSNP_sub)
}


 
query_by_coordinates <- function(top_SNPs, 
                                 gene, 
                                 subset_path, 
                                 fullSS_path,
                                 file_sep, 
                                 chrom_col, 
                                 position_col,
                                 min_POS, 
                                 max_POS, 
                                 bp_distance){
  gz.reader <- ifelse(endsWith(fullSS_path,".gz"), " gzcat ","")
  topSNP_sub <- top_SNPs[top_SNPs$Gene==gene & !is.na(top_SNPs$Gene),]
  if(is.na(min_POS)){min_POS <- topSNP_sub$POS - bp_distance}
  if(is.na(max_POS)){max_POS <- topSNP_sub$POS + bp_distance}
  printer("---Min snp position:",min_POS, "---")
  printer("---Max snp position:",max_POS, "---")
  colDict <- column_dictionary(fullSS_path)
  awk_cmd <- paste0(gz.reader,fullSS_path," | awk -F '",file_sep,"' 'NR==1 {print $0} NR>1 { if($",colDict[[chrom_col]]," == ",topSNP_sub$CHR,
                   " && ($", colDict[[position_col]]," >= ",min_POS," && $",colDict[[position_col]]," <= ",max_POS,")) {print $0} }'"," > ",subset_path)
  printer(awk_cmd)
  system(awk_cmd)
}

query_by_coordinates_merged <- function(top_SNPs, fullSS_path, subset_path, gene,
                                        chrom_col, file_sep=" ", location_sep=":",
                                        min_POS, max_POS, bp_distance){
  topSNP_sub <- top_SNPs[top_SNPs$Gene==gene & !is.na(top_SNPs$Gene),][1,]
  if(is.na(min_POS)){min_POS <- topSNP_sub$POS - bp_distance}
  if(is.na(max_POS)){max_POS <- topSNP_sub$POS + bp_distance}
  printer("---Min snp position:",min_POS, "---")
  printer("---Max snp position:",max_POS, "---")
  colDict <- column_dictionary(fullSS_path)
  awk_cmd <- paste("cat ",fullSS_path," | tr -s '",location_sep,"' '",file_sep,"'",
                   " | awk -F '",file_sep,"' 'NR==1 {print \"CHR POS \" $2\" \" $3\" \" $4\" \" $5\" \" $6 }",
                   " NR>1 {if($",colDict[[chrom_col]]," == ",topSNP_sub$CHR," && ",
                   "($",colDict[[chrom_col]]+1," >=",min_POS,"&& $",colDict[[chrom_col]]+1," <=",max_POS,")) {print $0}}'",
                   " | tr -s '",file_sep,"' '\t' > ",subset_path, sep="")
  # awk_cmd <- paste("cat ",fullSS_path," | tr -s '",location_sep,"' '",file_sep,"'",
  #                  " | awk -F '",file_sep,"' 'NR==1 {print \"CHR POS \" $2\" \" $3\" \" $4\" \" $5\" \" $6 }",
  #                  " NR>1 {if($1 == 1 && ($2 >=",min_POS,"&& $2 <=",max_POS,")) {print $0}}'",
  #                  " | tr -s '",location_sep,"' '\t' > ",subset_path, sep="")
  printer(awk_cmd)
  system(awk_cmd)
}

query_by_gene <- function(fullSS_path, subset_path, gene, gene_col, file_sep){
  colDict <- column_dictionary(fullSS_path)
  # if(endsWith(fullSS_path,".gz")){fullSS_path <- paste("<(gzcat",fullSS_path,")")} 
  awk_cmd <- paste("awk -F '",file_sep,"' 'NR==1{print $0} NR>1{if($",colDict[[gene_col]]," == \"",gene,"\"){print $0}}' ",fullSS_path,
                   "| tr -s '",file_sep,"' '\t'  > ",subset_path, sep="")
  printer(awk_cmd)
  system(awk_cmd)
}

query_by_probe <- function(fullSS_path, subset_path, gene, gene_col, chrom_col,
                           file_sep, probe_path, coordinates_merged=T, location_sep=":"){
  ## EXAMPLE COMMAND LINE
  # awk -F ' ' 'NR==1{print $0} NR>1{if($2 == "ILMN_2226015" || $2 == "ILMN_1776649") {print $0}}' cis.eqtls.fairfax.all.chr.IFN.47231.367.b.qced.f.txt > LRRK2_Fairfax_IFN.txt
  # probe_path="Data/eQTL/Fairfax/gene.ILMN.map"
  # subset_path="Data/eQTL/Fairfax/LRRK2_Fairfax_CD14.txt"
  # file_sep="\t"
  # gene="LRRK2"
  colDict <- column_dictionary(fullSS_path)
  probes <- find_probes(probe_path, gene)
  probe_string <- paste(paste(paste("$",colDict[[gene_col]], sep="")," == \"" ,probes,"\"", sep=""), collapse=" || ")
  
  awk_cmd <- paste("awk -F '",file_sep,"' 'NR==1 {print $0} NR>1 if(" ,probe_string,") {print}' ",fullSS_path,
                   " > ",subset_path, sep="")
  printer(awk_cmd)
  system(awk_cmd)
  
  if(coordinates_merged){
    ##  EXAMPLE COMMAND LINE 
    
    # awk -F ' ' 'NR==1{print "Coord","CHR","POS",$2,$3,$4,$5,$6 } NR>1{split($1,a,":"); print $1, a[1], a[2], $2, $3, $4, $5, $6}' LRRK2_Fairfax_CD14.txt | tr -s " " "\t" > tmp.txt && mv tmp.txt LRRK2_Fairfax_CD14.txt
    awk_cmd <- paste("awk -F '",file_sep,"' 'NR==1{print \"Coord\",\"CHR\",\"POS\",$2,$3,$4,$5,$6 }",
                     "NR>1{split($",colDict[[chrom_col]],",a,\":\"); print $1, a[1], a[2], $2, $3, $4, $5, $6}' ",
                     subset_path, " | tr -s ' ' '\t' > tmp.txt && mv tmp.txt ",subset_path, sep="")
    printer(awk_cmd)
    system(awk_cmd)
  }
  
}

query_fullSS <- function(fullSS_path, subset_path){
  file.copy(fullSS_path, subset_path)
}
   

coordinates_to_SNPs <- function(){
  # head /sc/orga/projects/ad-omics/satesh/HRC.RSID
  # SNP	RSID
  # 1:13380	rs571093408
  # 1:16071 rs541172944 
  array = paste('"',paste(data.table::fread("Data/eQTL/Fairfax/eQTL_effect_sizes.csv")$Coord, collapse='","'),'"',sep="")
  printer("awk -F '\t' 'NR==FNR{arr=[",array,"];next} NR==1{print $0} NR>1{if($1 in arr){print $0}}' /sc/orga/projects/ad-omics/satesh/HRC.RSID", sep="")
  # awk -F '\t' 'NR==1{print $0} NR>1{if($1 in ["12:40614434","12:40614434","12:40922572","12:40922572","12:40614434","12:40614434","12:40922572","12:40922572","12:40614434","12:40614434","12:40922572","12:40922572","12:40614434","12:40614434","12:40922572","12:40922572"]){print $0}}' /sc/orga/projects/ad-omics/satesh/HRC.RSID
  
}

# 
# manual_query <- function(fullSS_path, subset_path){
#   # Query via the Minerva R interface
#   gene_col="gene"
#   gene="LRRK2"
#   file_sep="\t"
#   column_dictionary <- function(file_path){
#     # Get the index of each column name
#     f <- data.table::fread(file_path, nrows = 0, header = T)
#     cNames <- colnames(f)
#     colDict <- setNames(1:length(cNames), cNames  )
#     return(colDict)
#   }
# 
#   query_by_probe <- function(fullSS_path, subset_path, gene, gene_col, file_sep){
#     colDict <- column_dictionary(fullSS_path)
#     # probe_path="Data/eQTL/Fairfax/gene.ILMN.map" 
#     file_sep="\t"
#     gene="LRRK2"
#     probes <- c("ILMN_2226015","ILMN_1776649")#find_probes(probe_path, gene)
#     probe_string <- paste(paste(paste("($",colDict[[gene_col]], sep="")," == \"" ,probes,"\")", sep=""), collapse=" || ")
#     
#     awk_cmd <- paste("awk -F '",file_sep,"' 'NR==1 {print $0} NR>1{if(" ,probe_string,") {print $0}}' ",fullSS_path,
#                      " > ",subset_path, sep="")
#     # awk_cmd <- paste("awk '/ILMN_2226015/'",fullSS_path,">",subset_path)
#     printer(awk_cmd)
#     system(awk_cmd)
#   }
#   query_by_probe(fullSS_path, subset_path, gene, gene_col, file_sep)
# }
# # manual_query("/hpc/users/rajt01/ad-omics/data/fairfax/sumstats/cis.eqtls.fairfax.all.chr.CD14.47231.414.b.qced.f.txt",
# #              "/hpc/users/rajt01/ad-omics/data/fairfax/sumstats/LRRK2_Fairfax_CD14.txt")


query_handler <- function(gene, 
                          fullSS_path, 
                          top_SNPs=NULL, 
                          subset_path, #top_SNPs="auto", 
                          min_POS=NA, 
                          max_POS=NA, 
                          bp_distance=500000,
                          gene_col="Gene",
                          chrom_col="CHR", 
                          position_col="POS",
                          file_sep="\t",
                          query_by="coordinates", 
                          probe_path = "./Data/eQTL/gene.ILMN.map"){ 
  printer("++ Query Method: '",query_by, sep="")
  if(query_by=="tabix"){
    topSNP_sub <- top_SNPs[top_SNPs$Gene==gene & !is.na(top_SNPs$Gene),]
    if(is.na(min_POS)){min_POS <- topSNP_sub$POS - bp_distance}
    if(is.na(max_POS)){max_POS <- topSNP_sub$POS + bp_distance}
    printer("---Min snp position:",min_POS, "---")
    printer("---Max snp position:",max_POS, "---")
    TABIX(fullSS_path=fullSS_path, 
          subset_path=subset_path,
          # is_tabix=F, 
          chrom_col=chrom_col, 
          position_col=position_col,
          min_POS=min_POS,
          max_POS=max_POS,
          chrom= gsub("chr","",topSNP_sub$CHR[1])
          )
  }
  if(query_by=="coordinates"){
    query_by_coordinates(top_SNPs=top_SNPs, gene=gene,
                         subset_path=subset_path, fullSS_path=fullSS_path,
                         chrom_col=chrom_col, position_col=position_col,
                         min_POS=min_POS, max_POS=max_POS, bp_distance=bp_distance,
                         file_sep=file_sep)
  }
  if(query_by=="coordinates_merged"){
    query_by_coordinates_merged(top_SNPs=top_SNPs, fullSS_path=fullSS_path, subset_path=subset_path,
                                file_sep=file_sep,  gene=gene, location_sep= ":",
                                min_POS=min_POS, max_POS=max_POS, bp_distance=bp_distance,
                                chrom_col=chrom_col)
  }
  if(query_by=="gene"){
    query_by_gene(fullSS_path=fullSS_path, subset_path=subset_path,
                  gene=gene, gene_col=gene_col, file_sep=file_sep)
  }
  if(query_by=="probes"){
    query_by_probe(fullSS_path=fullSS_path, subset_path=subset_path, 
                   gene=gene, gene_col=gene_col, chrom_col=chrom_col,
                   file_sep=file_sep, probe_path=probe_path)
  }
  if(query_by=="fullSS"){
    query_fullSS(fullSS_path=fullSS_path, 
                 subset_path = subset_path)
  }

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


 
# add_snp_info <-function(snpInfo_path, fullSS_path, newSS_path){
#   # sqldf::read.csv.sql(file=snpInfo_path, header = T, sep = "\t",
#   #                     sql = paste("select * from file where location IN (", paste("'",snp_locs,"'", collapse=",",sep=""),")",sep="")
#   #                     )
#   printer("\nLoading SNP Info file...\n")
#   snp_ids <- fread(snpInfo_path,header = T, stringsAsFactors = F, key="location", colClasses = rep("character",2))
#
#   printer("\nLoading full SS file...\n")
#   fullSS <- fread(fullSS_path, stringsAsFactors = F,
#                   colClasses = c(rep("character",3), rep("numeric",7), "character", rep("numeric",4)))
#   fullSS$MarkerName <- gsub("chr","",fullSS$MarkerName)
#   fullSS <- fullSS %>% rename(location="MarkerName")
#   fullSS <- data.table(fullSS, key="location")
#   dim(fullSS)
#   printer("\nMerging files...\n")
#   SS_merged <- snp_ids[fullSS]
#   dim(SS_merged)
#   printer("\nWriting new file...\n")
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
#   printer(awk_cmd)
#   system(awk_cmd)
# }

# 
# 
# ############ QTL ############
# # Subset eQTL data to just the parts you need
# subset_eQTL_SS <- function(fullSS_path, subset_path, gene,
#                            gene_col="Gene", chrom_col="CHR", position_col="POS",
#                            snp_col="SNP", pval_col="P", effect_col="Effect",
#                            bp_distance=500000, file_sep="\t",
#                            force_new_subset=T){
#   start_eQTL_sub <- Sys.time()
# 
#   if(file.exists(subset_path) & force_new_subset==F){
#     printer("Subset file already exists. Importing",subset_path,"...\n")
#     check_if_empty(subset_path)
#     # Extract lead SNPs per gene from subset data, using the subset SS as a reference file
#     top_SNPs <- data.table::fread(subset_path) %>% mutate(Gene=gene) %>%
#       arrange(P, desc(Effect)) %>% group_by(Gene) %>% slice(1)
#   } else {
#     colDict <- column_dictionary(fullSS_path)
#     # Extract subset with awk
#     printer("Extracting relevant variants from fullSS...\n")
#     if(chrom_col==position_col){
#       # Automatically detect when the chrom and position cols are merged into one and parse accordingly
#       query_by_location(fullSS_path = fullSS_path, file_sep = file_sep, location_sep = ":", chr = top_SNPs$CHR,
#                         min_POS = (top_SNPs$POS - bp_distance), max_POS =  (top_SNPs$POS + bp_distance))
#     } else {
#       # Otherwise, just subset by gene
#       query_by_gene(fullSS_path = fullSS_path, subset_path = subset_path,
#                     gene = gene, gene_col = gene_col)
#     }
#     check_if_empty(subset_path)
#     # Extract lead SNPs per gene from subset data, using the subset SS as a reference file
#     top_SNPs <- import_topSNPs(
#       file_path = subset_path,
#       chrom_col = chrom_col, position_col= position_col, snp_col=snp_col,
#       pval_col=pval_col, effect_col=effect_col, gene_col=gene_col,
#       caption= paste(superpop,": eQTL Summary Stats"))
#     # Make format/col names consistent
#     preprocess_subset(topSNP_sub = top_SNPs, subset_path = subset_path,
#                       chrom_col = chrom_col, position_col = position_col, snp_col = snp_col,
#                       pval_col = pval_col, effect_col = effect_col, stderr_col = stderr_col, file_sep = "\t",
#                       return_dt=F)
#   }
#   end_eQTL_sub <- Sys.time()
#   printer("\nExtraction completed in", round(end_eQTL_sub-start_eQTL_sub, 2),"seconds \n")
#   return(top_SNPs)
# }
 



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



