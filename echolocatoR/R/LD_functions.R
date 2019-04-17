# %%%%%%%%%%%%%%%%% #
####### LD ####### 
# %%%%%%%%%%%%%%%%% #


# * Nalls et al. 2018: Imputation Panel Notes
# + _"One of the limitations of this study is the use of multiple imputation panels, due to logistic constraints.
# Adding datasets from non-European populations would be helpful to further improve our granularity in association
# testing and ability to fine-map loci through integration of more variable LD signatures."_
# + _"Post-Chang 23andMe samples were imputed using a combination of Finch for phasing (an in-house developed fork of Beagle)
# and miniMac2 for imputation with all-ethnicity samples from the __September 2013 release of 1000 Genomes Phase1__
# as reference haplotypes."_
# + _"The Nalls et al . 2014 and Chang et al . 2017 samples were imputed with Minimac2 using
# __1000 Genomes phase 1 haplotypes__.
# All additional sample series except for the post-Chang et al . 2017 samples from 23andMe were imputed using the
# __Haplotype Reference Consortium (HRC)__  on the University of Michigan imputation server under default settings
# with Eagle v2.3 phasing based on reference panel HRC r1.1 2016"_

download_all_vcfs <- function(vcf_folder="../1000_Genomes_VCFs"){
  # PHASE 3 DATA
  path3 <- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
  for(chrom in c(1:22)){
    cat("\nDownloading Chromosome",chrom,"\n")
    URL <- paste("ALL.chr",chrom,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",sep = "")
    system(paste("wget -P",file.path(vcf_folder,"Phase3"), file.path(path3, URL) ))
  }
  X_chrom <-"ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz"
  system(paste("wget -P",file.path(vcf_folder,"Phase3"), file.path(path3, X_chrom)))
  Y_chrom <- "ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz"
  system(paste("wget -P",file.path(vcf_folder,"Phase3"), file.path(path3, Y_chrom) ))
  
  popDat_URL = file.path(path3, "integrated_call_samples_v3.20130502.ALL.panel")
  popDat <- read.delim(popDat_URL, header = F, row.names = NULL)
  write.table(popDat,file=file.path(vcf_folder,"Phase3","integrated_call_samples_v3.20130502.ALL.panel"), row.names = F, sep="\t", quote = F, col.names = F)
  
  # PHASE 1 DATA
  path1 <- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521"
  for(chrom in c(1:22)){
    cat("\nDownloading Chromosome",chrom,"\n")
    URL <- paste("ALL.chr",chrom, ".phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz", sep="")
    system(paste( "wget -P",file.path(vcf_folder,"Phase1"), file.path(path1, URL) ))
  }
  X_chrom <- "ALL.chrX.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz"
  system( paste("wget -P",file.path(vcf_folder,"Phase1"), file.path(path1, X_chrom)) )
  
  popDat_URL = file.path(path1, "phase1_integrated_calls.20101123.ALL.panel")
  popDat <- read.delim(popDat_URL, header = F, row.names = NULL)
  write.table(popDat,file=file.path(vcf_folder,"Phase1","phase1_integrated_calls.20101123.ALL.panel"),  row.names = F, sep="\t", quote = F, col.names = F)
}



gaston_LD <- function(flankingSNPs, gene, reference="1KG_Phase1", superpopulation="EUR", 
                      vcf_folder=F, min_r2=0, LD_block=F, block_size=.7, min_Dprime=F){
  # Download portion of vcf from 1KG website
  region <- paste(flankingSNPs$CHR[1],":",min(flankingSNPs$POS),"-",max(flankingSNPs$POS), sep="")
  chrom <- flankingSNPs$CHR[1]
  # PHASE 3 DATA
  if(reference=="1KG_Phase3"){
    cat("LD Reference Panel = 1KG_Phase3 \n")
    if(vcf_folder==F){## With internet
      vcf_URL <- paste("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr",chrom,
                       ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",sep="")
      popDat_URL = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
    }else{## WithOUT internet
      vcf_URL <- paste(vcf_folder, "/ALL.chr",chrom,
                       ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",sep="")
      popDat_URL = file.path(vcf_folder,"integrated_call_samples_v3.20130502.ALL.panel")
    }
    
    # PHASE 1 DATA
  } else if (reference=="1KG_Phase1") {
    cat("LD Reference Panel = 1KG_Phase1 \n")
    if(vcf_folder==F){## With internet
      vcf_URL <- paste("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr",chrom,
                       ".phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz", sep="")
      popDat_URL = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/phase1_integrated_calls.20101123.ALL.panel"
    }else{## WithOUT internet
      vcf_URL <- paste(vcf_folder,"/ALL.chr",chrom,
                       ".phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz", sep="")
      popDat_URL = file.path(vcf_folder, "phase1_integrated_calls.20101123.ALL.panel")
    }
  }
  
  phase <- gsub("1KG_","",reference) 
  popDat <- read.delim(popDat_URL, header = F, row.names = NULL)
  colnames(popDat) <- c("sample","population","superpop","gender")
  
  # library(Rsamtools); #BiocManager::install("Rsamtools")
  subset_vcf <- file.path("../1000_Genomes_VCFs", phase, paste(gene,"subset.vcf",sep="_"))
  # Create directory if it doesn't exist
  if(!dir.exists(dirname(dirname(subset_vcf))) ) {
    dir.create(path = dirname(subset_vcf),recursive =  T)
  }else(cat("Creating '../1000_Genomes_VCFs' directory.\n"))
  # Download and subset vcf if the subset doesn't exist already
  if(!file.exists(subset_vcf)){
    tabix_cmd <- paste("tabix -fh",vcf_URL, region, ">", subset_vcf)
    cat(tabix_cmd)
    system(tabix_cmd)
    vcf_name <- paste(basename(vcf_URL), ".tbi", sep="")
    file.remove(vcf_name)
  }else{cat("Identified matching VCF subset file. Importing...", subset_vcf,"\n")} 
  # Import w/ gaston and further subset
  cat("\n------ Importing VCF as bed file... ------\n")
  bed.file <- gaston::read.vcf(subset_vcf, verbose = F) 
  ## Subset rsIDs
  bed <- gaston::select.snps(bed.file, id %in% flankingSNPs$SNP & id !=".")
  
  gaston::write.bed.matrix(bed, "./plink_tmp/matrix")
  # Subset Individuals
  selectedInds <- subset(popDat, superpop == superpopulation)
  bed <- gaston::select.inds(bed, id %in% selectedInds$sample)
  # Cleanup extra files
  remove(bed.file)
  # file.remove("subset.vcf")
  
  # Calculate pairwise LD for all SNP combinations
  #### "Caution that the LD matrix has to be correlation matrix" -SuSiER documentation
  ### https://stephenslab.github.io/susieR/articles/finemapping_summary_statistics.html  
  # Gaston LD method
  # LD_matrix <- gaston::LD(bed, lim = c(1,ncol(bed)), measure ="r") #"D"
  # LD_matrix[!is.finite(LD_matrix)] <- 0
  
  # Get lead SNP rsid
  leadSNP = subset(flankingSNPs, leadSNP==T)$SNP #rs76904798
  # Plink LD method
  LD_matrix <- calculate_LD(leadSNP, min_r2 = min_r2, min_Dprime = min_Dprime) 
  # Make sure only includes SNPs in flankingSNPs
  # LD_matrix <- LD_matrix[row.names(LD_matrix) %in% flankingSNPs$SNP, colnames(LD_matrix) %in% flankingSNPs$SNP]
  
  # Filter by r2 with lead SNP
  # filter_by_LD <- function(LD_matrix, leadSNP, min_r2=0.2){
  #   if(min_r2>0){
  #     l = LD_matrix[leadSNP,]
  #     ld_filt <- l[lapply(l, function(x){x^2}) >= min_r2] %>% names()
  #     # LD_matrix[LD_matrix^2>=.2,]
  #     return(LD_matrix[ld_filt, ld_filt])
  #   } else{return(LD_matrix)}
  # }
  # LD_matrix <- filter_by_LD(LD_matrix, leadSNP, min_r2)
  
  # Filter out SNPs not in the same LD block as the lead SNP
  if(LD_block==T){
    block_snps <- leadSNP_block(leadSNP, "./plink_tmp", block_size)
    LD_matrix <- LD_matrix[row.names(LD_matrix) %in% block_snps, colnames(LD_matrix) %in% block_snps]
  }
  
  # LD plot
  try({ 
    lead_index = match(leadSNP, row.names(LD_matrix)) 
    if(dim(LD_matrix)[1]<10){
      start = lead_index - dim(LD_matrix)[1]
      end = lead_index + dim(LD_matrix)[1]
    } else{
      start = lead_index - 10
      end = lead_index + 10
    } 
    gaston::LD.plot( LD_matrix[start:end, start:end], snp.positions = bed@snps$pos[start:end] )
  })
  # Double check subsetting
  # LD_matrix <- ld.x[row.names(ld.x) %in% flankingSNPs$SNP, colnames(ld.x) %in% flankingSNPs$SNP]
  return(LD_matrix)
}
# LD_matrix <- gaston_LD(flankingSNPs)

Dprime_table <- function(SNP_list, plink_folder="./plink_tmp"){
  system( paste("./echolocatoR/tools/plink1.9 --bfile",file.path(plink_folder,"matrix"),"--ld-snps", paste(SNP_list, collapse=" "),
                "--r dprime-signed --ld-window 10000000 --ld-window-kb 10000000 --out",file.path(plink_folder,"plink")) ) 
  #--ld-window-r2 0
  plink.ld <- data.table::fread(file.path(plink_folder, "plink.ld"), select = c("SNP_A", "SNP_B","DP","R"))
  plink.ld$R2 <- plink.ld$R^2
  # plink.ld <- plink.ld[complete.cases(plink.ld) ] 
  return(plink.ld)
}

complex_LD <- function(bim, plink_folder="./plink_tmp", min_Dprime=F){
  # METHOD 1
  plink.ld <- Dprime_table(bim)
  # Filter NaNs
  plink.ld <- plink.ld[!is.nan(plink.ld$DP),]
  plink.ld <- plink.ld[!is.nan(plink.ld$R),]
  
  if(min_Dprime!=F){
    cat("\n++++++++++ Filtering by DPrime ++++++++++\n")
    plink.ld <- subset(plink.ld, DP>=min_Dprime)
  }
  ld.matrix <- data.table::dcast.data.table(plink.ld, formula = SNP_B ~ SNP_A, value.var="R",
                                            fill=0, drop=F, fun.aggregate = mean)
  ld.matrix <-  data.frame(ld.matrix, row.names = ld.matrix$SNP_B) %>% subset(select = -SNP_B) %>% 
    data.table() %>% as.matrix()
  return(ld.matrix)
}

simple_LD <-function(bim, plink_folder="./plink_tmp"){
  # METHOD 2 (faster, but less control over parameters. Most importantly, can't get Dprime)
  system( paste("./echolocatoR/tools/plink1.9 --bfile",file.path(plink_folder,"matrix"),
                "--extract",file.path(plink_folder,"SNPs.txt"),
                "--r square bin --out", file.path(plink_folder,"plink")) )
  bin.vector <- readBin(file.path(plink_folder, "plink.ld.bin"), what = "numeric", n=length(bim$SNP)^2)
  ld.matrix <- matrix(bin.vector, nrow = length(bim$SNP), dimnames = list(bim$SNP, bim$SNP))
  return(ld.matrix)
}


calculate_LD <-function(leadSNP, plink_folder="./plink_tmp", min_r2=0, min_Dprime=F){
  # Dprime ranges from -1 to 1
  start <- Sys.time()
  # Calculate LD 
  cat("\n++++++++++ Reading in BIM file... ++++++++++\n")
  bim <- data.table::fread(file.path(plink_folder, "matrix.bim"), col.names = c("CHR","SNP","V3","POS","A1","A2")) 
  data.table::fwrite(subset(bim, select="SNP"), file.path(plink_folder,"SNPs.txt"), col.names = F)
  
  cat("\n++++++++++ Calculating LD ++++++++++\n")
  ld.matrix <- simple_LD(bim, plink_folder)
  
  if((min_Dprime != F) | (min_r2!=0) ){ 
    plink.ld <- Dprime_table(SNP_list = row.names(ld.matrix), plink_folder) 
    keep.pairs <- plink.ld
    # DPrime filter
    if(min_Dprime != F){
      keep.pairs <- subset(removed.pairs, (SNP_A==leadSNP & DP>=min_Dprime) | (SNP_B==leadSNP & DP>=min_Dprime))
    }
    # R2 filter
    if(min_r2!=0){
      keep.pairs <- subset(plink.ld, (SNP_A==leadSNP & R2>=min_r2) | (SNP_B==leadSNP & R2>=min_r2))  
    snp_list <-  unique(keep.pairs$SNP_A, keep.pairs$SNP_A)
    ld.matrix <- ld.matrix[row.names(ld.matrix) %in% snp_list, colnames(ld.matrix) %in% snp_list] 
  } 
    # Apply filters
    ## Manually remove rare variant
    # ld.matrix <- ld.matrix[rownames(ld.matrix)!="rs34637584", colnames(ld.matrix)!="rs34637584"]
    
    # !IMPORTANT!: Fill NAs (otherwise susieR will break)
    ld.matrix[is.na(ld.matrix)] <- 0
    end <- Sys.time()
    cat("\n++++++++++ LD matrix calculated in",round(as.numeric(end-start),2),"seconds. ++++++++++\n")
    return(ld.matrix)
  }
}

LD_blocks <- function(plink_folder="./plink_tmp", block_size=.7){
  cat("\n++++++++++ Calculating LD blocks... ++++++++++\n")
  # PLINK 1.07 LD: http://zzz.bwh.harvard.edu/plink/ld.shtml
  # PLINK 1.9 LD: https://www.cog-genomics.org/plink/1.9/ld
  # system("./echolocatoR/tools/plink1.9 -h")
  # Identify duplicate snps
  # system("./echolocatoR/tools/plink1.9 --vcf subset.vcf --list-duplicate-vars")
  # Convert vcf to plink format
  # system("./echolocatoR/tools/plink1.9 --vcf subset.vcf --exclude ./plink_tmp/plink.dupvar --make-bed --out PTK2B")
  
  # Estimate LD blocks
  # Defaults: --blocks-strong-lowci = 0.70, --blocks-strong-highci .90
  
  # Reucing "--blocks-inform-frac" is the only parameter that seems to make the block sizes larger
  system( paste("./echolocatoR/tools/plink1.9 --bfile",file.path(plink_folder,"matrix"),
                "--blocks no-pheno-req no-small-max-span --blocks-max-kb 100000",
                # "--blocks-strong-lowci .52 --blocks-strong-highci 1",
                "--blocks-inform-frac",block_size," --blocks-min-maf 0 --out",file.path(plink_folder,"plink")) )
  # system( paste("./echolocatoR/tools/plink1.9 --bfile plink --ld-snp-list snp_list.txt --r") )
  blocks <- data.table::fread("./plink_tmp/plink.blocks.det") 
  return(blocks)
}


leadSNP_block <- function(leadSNP, plink_folder="./plink_tmp", block_size=.7){
  cat("\n Returning lead SNP's block...\n")
  blocks <- LD_blocks(plink_folder, block_size)
  splitLists <- strsplit(blocks$SNPS,split = "[|]")
  block_snps <- lapply(splitLists, function(l, leadSNP){if(leadSNP %in% l){return(l)} }, leadSNP=leadSNP) %>% unlist()
  cat("\n Number of SNPs in LD block =", length(block_snps), "\n")
  return(block_snps)
}


# LD_clumping <- function(subset_vcf, subset_SS){ 
#   # PLINK clumping: http://zzz.bwh.harvard.edu/plink/clump.shtml
#   # Convert vcf to .map (beagle)
#   ## https://www.cog-genomics.org/plink/1.9/data
#   system(paste("./echolocatoR/tools/plink1.9 --vcf",subset_vcf,"--recode beagle --out ./plink_tmp/plink"))
#   # Clumping
#   system("./echolocatoR/tools/plink1.9 --file ./plink_tmp/plink.chr-8 --clump",subset_SS,"--out ./plink_tmp")
# }