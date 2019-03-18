#### ------- echolocatoR ------- ####

# Author: Brian M. Schilder
  # Bioinformatician II
  # Icahn School of Medicine at Mount Sinai
  # New York City, New York, USA
  # https://bschilder.github.io/BMSchilder

     # =/\                  /\=
    #  / \'._   (\_/)   _.'/ \
   #  / .''._'--(o.o)--'_.''. \
  #  /.' _/ |`'=/ " \='`| \_ `.\
 #  /` .' `\;-,'\___/',-;/` '. '\
#  /.-'       `\(-V-)/`       `-.\
# `              "   "           `


# You can learn more about package authoring with RStudio at:
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

# Load libraries
.libPaths()

library(readxl)
library(DT)
library(data.table)
library(dplyr)
library(ggplot2)
library(plotly)
library(cowplot)
library(ggrepel)
library(curl)
library(biomaRt)
library(sqldf)
# Ensembl LD API
library(httr)
library(jsonlite)
library(xml2)
library(gaston)
library(RCurl)
library(tidyr)
library(biomaRt)

# *** susieR ****
# library(knitrBootstrap) #install_github('jimhester/knitrBootstrap')
library(susieR) # devtools::install_github("stephenslab/susieR")

# *** finemapr ****
## finemapr contains: finemap, CAVIAR, and PAINTOR
# library(finemapr) # devtools::install_github("variani/finemapr")
# library(roxygen2) #roxygenize()

# *** locuscomparer ****
# https://github.com/boxiangliu/locuscomparer
# library(locuscomparer); #devtools::install_github("boxiangliu/locuscomparer")

# thm <- knitr::knit_theme$get("bipolar")
# knitr::knit_theme$set(thm)



## General Functions
createDT <- function(DF, caption="", scrollY=400){
  data <- DT::datatable(DF, caption=caption,
                        extensions = 'Buttons',
                        options = list( dom = 'Bfrtip',
                                        buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                                        scrollY = scrollY, scrollX=T, scrollCollapse = T, paging = F,
                                        columnDefs = list(list(className = 'dt-center', targets = "_all"))
                        )
  )
  return(data)
}
createDT_html <- function(DF, caption="", scrollY=400){
  htmltools::tagList( createDT(DF, caption, scrollY))
}

# https://stackoverflow.com/questions/42866055/rshiny-download-button-within-rmarkdown
# downloadHandler(filename = function() {
#     return(paste('Example', input$SS, '.csv', sep=''))
#
#  }, content = function(file) {
#    write.csv(RandomSample(), file)
#  })


# Data Preprocessing
snps_to_genes <- function(snp_list){
  # listMarts()
  snp_mart = useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp")
  # View(listFilters(snp_mart))
  # View(listAttributes(snp_mart))
  snp_results <- getBM(snp_mart, filters="snp_filter", values=snp_list,
                       attributes=c("refsnp_id","snp","chr_name", "chrom_start",
                                    "associated_gene","ensembl_gene_stable_id" ) )
  # Split ensembl IDs
  snp_results <-snp_results %>%
    mutate(ensembl = strsplit(as.character(ensembl_gene_stable_id), ";")) %>%
    tidyr::unnest(ensembl)
  # listDatasets(useMart("ENSEMBL_MART_ENSEMBL") )
  gene_mart = useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
  # View(listFilters(gene_mart))
  # View(listAttributes(gene_mart))
  gene_results <- getBM(mart = gene_mart, filters = "ensembl_gene_id",
                        values = unlist(strsplit(snp_results$ensembl, ";")),
                        attributes = c("external_gene_name","ensembl_gene_id") )
  merged_df <- data.table(gene_results, key = "ensembl_gene_id")[data.table(snp_results, key = "ensembl")]
  return(merged_df)
}


#
# locations_to_snps <- function(locations){
#   # listMarts()
#   locations = gsub(pattern = "chr", replacement = "", locations)
#   # CHR = lapply(locations, function(x){strsplit(x, ":")[[1]][1]}) %>% unlist()
#   # POS = lapply(locations, function(x){strsplit(x, ":")[[1]][2]}) %>% unlist()
#
#   snp_mart = useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp")
#   # View(listFilters(snp_mart))
#   # View(listAttributes(snp_mart))
#   snp_results <- getBM(snp_mart, filters=c("chr_name", "start"), values=list("1", "60320992"),
#                        attributes=c("refsnp_id","snp","chr_name", "chrom_start",
#                                     "associated_gene","ensembl_gene_stable_id" ) )
#   # Split ensembl IDs
#   snp_results <-snp_results %>%
#     mutate(ensembl = strsplit(as.character(ensembl_gene_stable_id), ";")) %>%
#     tidyr::unnest(ensembl)
#   # listDatasets(useMart("ENSEMBL_MART_ENSEMBL") )
#   gene_mart = useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
#   # View(listFilters(gene_mart))
#   # View(listAttributes(gene_mart))
#   gene_results <- getBM(mart = gene_mart, filters = "ensembl_gene_id",
#                         values = unlist(strsplit(snp_results$ensembl, ";")),
#                         attributes = c("external_gene_name","ensembl_gene_id") )
#   merged_df <- data.table(gene_results, key = "ensembl_gene_id")[data.table(snp_results, key = "ensembl")]
#   return(merged_df)
# }


## TSS Data
# Use bioMart to get TSS positions for each gene
get_TSS_position <- function(gene){
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  # att <- listAttributes(mart)
  # grep("transcription", att$name, value=TRUE)
  TSS <- getBM(mart=mart,
               attributes=c("hgnc_symbol","transcription_start_site", "version"),
               filters="hgnc_symbol", values=gene)
}



## Import Sig GWAS Summary Statistics
# For each gene, get the position of top SNP (the one with the greatest effect size/Beta)
import_sig_GWAS <- function(file_path, caption="", sheet = 1,
                            chrom_col="CHR", position_col="POS", snp_col="SNP",
                            pval_col="P", effect_col="Effect", gene_col="Gene"){
  ## Only the significant subset of results
  top_SNPs <- read_excel(path = file_path, sheet = sheet)[,c(chrom_col, position_col, snp_col,
                                                                 pval_col, effect_col, gene_col )]
  # Standardize names
  colnames(top_SNPs) <- c("CHR","POS","SNP","P","Effect","Gene")
  top_SNPs <- top_SNPs %>% group_by(Gene) %>% top_n(-1, "P") #`Beta, all studies`
  top_SNPs$CHR <- gsub("chr", "",top_SNPs$CHR)
  top_SNPs <- cbind(Coord= paste(top_SNPs$CHR, top_SNPs$POS, sep=":"),
                    top_SNPs)
  createDT(top_SNPs, caption)
  return(top_SNPs)
}
# top_SNPs <- import_sig_GWAS(
#   file_path = Data_dirs$Nalls_2019$topSS,
#   sheet="Data",
#   chrom_col = "CHR", position_col = "BP", snp_col="SNP",
#   pval_col="P, all studies", effect_col="Beta, all studies", gene_col="Nearest Gene",
#   caption= "Nalls et al. (2018) PD GWAS Summary Stats")
#
#
# finemapped_PD <- finemap_geneList(top_SNPs, geneList=c("LRRK2","GBAP1","SNCA","VPS13C","GCH1"),
#                                   file_path=Data_dirs$Nalls_2018$fullSS,
#                                   snp_col = "MarkerName", pval_col = "P.value")

## Get Flanking SNPs

# Get all genes surrounding the index SNP (default is 500kb upstream + 500kb downstream)
# 1000000 bp
get_flanking_SNPs <- function(gene, top_SNPs, bp_distance=500000, file_path,
                              chrom_col="CHR", position_col="POS", snp_col="SNP",
                              pval_col="P", effect_col="Effect", stderr_col="StdErr"){
  read.delim(file_path, nrows = 2)

  if(gene %in% top_SNPs$Gene==F){
    cat("\n ----- Could not find gene '",gene,"' in topSNPs dataframe. Try a different gene name. ----- \n" )
  }else{
    topSNP_sub <- top_SNPs[top_SNPs$Gene==gene & !is.na(top_SNPs$Gene),][1,]
    minPos <- as.numeric(topSNP_sub$POS) - bp_distance
    maxPos <- as.numeric(topSNP_sub$POS) + bp_distance
    query <- sqldf::read.csv.sql(file_path, sep="\t", stringsAsFactors=F,
                                 sql = paste('select * from file where',chrom_col,'=',topSNP_sub$CHR,
                                             'AND',position_col,'BETWEEN', minPos, 'AND', maxPos))
    query <- query[,c(chrom_col,position_col, snp_col, pval_col, effect_col, stderr_col)]
    geneSubset <- query %>% dplyr::rename(CHR=chrom_col,POS=position_col, SNP=snp_col,
                                          P=pval_col, Effect=effect_col, StdErr=stderr_col) %>%
      mutate(Effect=as.numeric(Effect), StdErr=as.numeric(StdErr), P=as.numeric(P))
    ## Remove SNPs with NAs in stats
    geneSubset <- geneSubset[complete.cases(geneSubset),]
    geneSubset$leadSNP <- ifelse(geneSubset$SNP==topSNP_sub$SNP, T, F)
    return(geneSubset)
    # Close the connection
    base::close(query)
    unlink(file_path)
  }
}

# flankingSNPs <- get_flanking_SNPs("LRRK2", top_SNPs,
#                                   file_path = "Data/Parkinsons/META.PD.NALLS2014.PRS.tsv",
#                                   snp_col = "MarkerName", pval_col = "P.value")
# flankingSNPs <- get_flanking_SNPs("CLU/PTK2B", top_SNPs,
#                                   file_path = "Data/Alzheimers/Posthuma/phase3.beta.se.hrc.txt",
#                                   effect_col = "BETA", stderr_col = "SE", position_col = "BP")
# flankingSNPs <- get_flanking_SNPs("LRRK2", top_SNPs,
#                                   file_path = "Data/Parkinsons/23andMe/PD_all_post30APRIL2015_5.2_extended.txt",
#                                   snp_col = "SNP", pval_col = "pvalue", chrom_col = "CHR",
#                                   position_col = "position", effect_col = "effect",stderr_col = "stderr")
#
#



# susieR

## Gaston LD
#
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

gaston_LD <- function(flankingSNPs, reference="1KG_Phase1", superpopulation="EUR", vcf_folder=F){
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

  popDat <- read.delim(popDat_URL, header = F, row.names = NULL)
  colnames(popDat) <- c("sample","population","superpop","gender")

  # library(Rsamtools); #BiocManager::install("Rsamtools")
  system(paste("tabix -fh ",vcf_URL ,region, "> subset.vcf"))
  vcf_name <- paste(basename(vcf_URL), ".tbi",sep="")
  file.remove(vcf_name)
  # Import w/ gaston and further subset
  bed.file <- read.vcf("subset.vcf")
  ## Subset rsIDs
  bed <- select.snps(bed.file, id %in% flankingSNPs$SNP)
  # Subset Individuals
  selectedInds <- subset(popDat, superpop == superpopulation)
  bed <- select.inds(bed, id %in% selectedInds$sample)
  # Cleanup extra files
  remove(bed.file)
  file.remove("subset.vcf")

  # Calculate pairwise LD for all SNP combinations
  #### "Caution that the LD matrix has to be correlation matrix" -SuSiER documentation
  ### https://stephenslab.github.io/susieR/articles/finemapping_summary_statistics.html
  ld.x <- gaston::LD(bed, lim = c(1,ncol(bed)), measure ="r" )
  LD_matrix <- ld.x
  LD_matrix[!is.finite(LD_matrix)] <- 0
  # LD plot
  try({
    leadSNP = subset(flankingSNPs, leadSNP==T)$SNP
    lead_index = match(leadSNP, rownames(LD_matrix))
    start = lead_index-10
    end = lead_index+10
    LD.plot( LD_matrix[start:end, start:end], snp.positions = bed@snps$pos[start:end] )
  })
  # Double check subsetting
  # LD_matrix <- ld.x[row.names(ld.x) %in% flankingSNPs$SNP, colnames(ld.x) %in% flankingSNPs$SNP]
  return(LD_matrix)
}
# LD_matrix <- gaston_LD(flankingSNPs)


## susieR Function
#
# * Notes on L parameter
# + L is the expected number of causal variants
# + Increasing L increases computational time
# + L=1: Gives a good amount of variation in PIP.
# + L=2: Warns "IBSS algorithm did not converge in 100 iterations!", but gives good variation in PIP.
# + L=3: Warns "IBSS algorithm did not converge in 100 iterations!". All PIPs 1s and 0s.
# + These results seem to be at least partially dependent on whether the ethnic composition of the LD matrix.
# * Notes on variance:
#   + If 'estimate_residual_variance' = TRUE _without_ providing 'var_y' _and_ L>1, susieR will throw error:
#   __"Estimating residual variance failed: the estimated value is negative"__
# + Running susieR with 'var_y = var(b)' provides _exactly_ the same results.
# * Statistical Terms:
#   + posterior inclusion probability (PIP)
# + coefficient estimate (Beta)
# + Effect allele frequency (EAF)
# + The I^2 statistic describes the percentage of variation across studies that seems not to be due to chance.
susie_on_gene <- function(gene, top_SNPs,
                          bp_distance=500000, file_path, num_causal=1,
                          chrom_col="CHR", position_col="POS", snp_col="SNP",
                          pval_col="P", effect_col="Effect", stderr_col="StdErr",
                          LD_reference="1KG_Phase1", superpopulation="EUR", vcf_folder=F){
  cat("\n + Extracting SNPs flanking lead SNP... \n")
  flankingSNPs <- get_flanking_SNPs(gene, top_SNPs, bp_distance=bp_distance, file_path=file_path,
                                    chrom_col=chrom_col, position_col=position_col, snp_col=snp_col,
                                    pval_col=pval_col, effect_col=effect_col, stderr_col=stderr_col)
  ### Get LD matrix
  cat("\n + Creating LD matrix... \n")
  LD_matrix <- gaston_LD(flankingSNPs, LD_reference, superpopulation, vcf_folder = vcf_folder)
  ## Turn LD matrix into positive semi-definite matrix
  # LD_matrix2 <- ifelse(matrixcalc::is.positive.semi.definite(LD_matrix),
  #        LDmatrix,
  #        Matrix::nearPD(LD_matrix)$mat %>% as.matrix() )

  ## Subset summary stats to only include SNPs found in query
  geneSubset <- flankingSNPs %>% subset(SNP %in% unique(row.names(LD_matrix), colnames(LD_matrix) ) )
  LD_matrix <- LD_matrix[geneSubset$SNP,  geneSubset$SNP]

  b <- geneSubset$Effect
  se <- geneSubset$StdErr

  # Run Susie
  cat("\n + Fine mapping with SusieR... \n")
  fitted_bhat <- susie_bhat(bhat = b, shat = se,
                            R = LD_matrix,
                            n = nrow(LD_matrix),

                            L = num_causal, # 1
                            # scaled_prior_variance = 0.1,
                            estimate_residual_variance = T, # TRUE
                            estimate_prior_variance = T, # FALSE
                            verbose = T,

                            # var_y = var(b),
                            standardize = T
  )
  # Format results
  credible_set <- row.names(LD_matrix)[as.numeric(summary(fitted_bhat)$cs$variable)]

  geneSubset$Coord <- paste(geneSubset$CHR, geneSubset$POS, sep=":")
  susieDF <- data.frame(SNP=names(fitted_bhat$X_column_scale_factors), PIP=fitted_bhat$pip) %>%
    base::merge(subset(geneSubset, select=c("CHR","POS","SNP","Effect","P","Coord","leadSNP")), by="SNP") %>%
    mutate(POS=as.numeric(POS))
  susieDF$credible_set <- ifelse(susieDF$SNP==credible_set, T, F)

  return(susieDF)
}
# susieDF <- susie_on_gene("LRRK2", top_SNPs,
#                          file_path = "Data/Parkinsons/META.PD.NALLS2014.PRS.tsv",
#                           snp_col = "MarkerName", pval_col = "P.value")
# susieDF <- susie_on_gene(gene="CLU/PTK2B", top_SNPs,
#                          file_path="Data/Alzheimers/Posthuma/phase3.beta.se.hrc.txt",
#                          effect_col = "BETA", stderr_col = "SE", position_col = "BP")

## Before-After Plots
before_after_plots <- function(gene, susieDF, before_var="P"){
  roundBreaks <- seq(plyr::round_any(min(susieDF$POS),10000), max(susieDF$POS),250000)
  ## Before fine-mapping
  if(before_var=="Effect"){
    geneSelect <- susieDF %>% arrange(desc(abs(Effect)))
  }else{
    # Sort by pval and then absolute Effect size
    geneSelect <- susieDF %>% arrange()
  }
  # Label original and fine mapped top snps
  leadSNP_before <- subset(geneSelect, leadSNP==T)[1,]
  leadSNP_after <- susieDF %>% arrange(desc(PIP))
  leadSNP_after <- leadSNP_after[1,]
  labelSNPs <- rbind(leadSNP_before, leadSNP_after)
  labelSNPs$type <- c("before", "after")
  labelSNPs$color <- c("red", "green")

  yLimits <- c(min(-log10(susieDF[before_var])), max(-log10(susieDF[before_var]))*1.1)
  yLimits <- c(0,-log10(leadSNP_before$P)+5)

  before_plot <- ggplot(geneSelect, aes(x=POS, y=-log10(eval(parse(text=before_var))), label=SNP, color= -log10(P) )) +
    ylim(yLimits) +
    geom_hline(yintercept=0, alpha=.5, linetype=1, size=.5) +
    geom_point(alpha=.5) +
    geom_point(data=labelSNPs[1,], pch=21, fill=NA, size=4, colour=labelSNPs[1,"color"], stroke=1) +
    geom_point(data=labelSNPs[2,], pch=21, fill=NA, size=4, colour=labelSNPs[2,"color"], stroke=1) +
    geom_segment(aes(xend=POS, yend= yLimits[1], color= -log10(P)), alpha=.5) +
    geom_text_repel(data=labelSNPs, aes(label=SNP), color=labelSNPs$color, segment.alpha = .2, nudge_x = .5) +
    labs(title=paste(gene," (",length(susieDF$PIP)," variants)","\nBefore Fine Mapping",sep=""),
         y=paste("-log10(p-value)",sep=""), x="Position", color="-log10(p-value)") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(breaks = roundBreaks)

  ## After fine-mapping
  susieDF <- susieDF %>% arrange(desc(PIP))
  yLimits <- c(min(susieDF$PIP), max(susieDF$PIP)+.1)

  after_plot <- ggplot(susieDF, aes(x=POS, y=PIP, label=SNP, color= -log10(P) )) +
    ylim(yLimits) +
    geom_hline(yintercept=0,alpha=.5, linetype=1, size=.5) +
    geom_point(alpha=.5) +
    geom_point(data=labelSNPs[1,], pch=21, fill=NA, size=4, colour=labelSNPs[1,"color"], stroke=1) +
    geom_point(data=labelSNPs[2,], pch=21, fill=NA, size=4, colour=labelSNPs[2,"color"], stroke=1) +
    geom_segment(aes(xend=POS, yend=yLimits[1], color= -log10(P)), alpha=.5) +
    geom_text_repel(data=labelSNPs, aes(label=SNP), color=labelSNPs$color, segment.alpha = .2, nudge_x = .5) +
    labs(title=paste(gene," (",length(susieDF$PIP)," variants)","\nAfter Fine Mapping",sep=""), y="PIP", x="Position",
         color="-log10(p-value)") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(breaks = roundBreaks)

  plot_grid(before_plot, after_plot, nrow = 2) %>% print()
  # susie_plot(fitted_bhat, y="PIP", b=b, add_bar = T)

  createDT_html(susieDF, paste("susieR Results: ", gene), scrollY = 200)
}
# before_after_plots(gene = "LRRK2", susieDF, topVariants = 3)
# before_after_plots(gene = "CLU/PTK2B", susieDF, topVariants = 3)

## Report SNP Overlap
before_after_consensus <- function(gene, top_SNPs, susieDF, max_SNPs=10){
  # Get top SNPs from Nalls et all fine mapping
  SS_dat  <- subset(top_SNPs, Gene==gene) %>%
    arrange(P) %>% mutate(SNP=as.character(SNP))
  SumStats_SNPs <-if(max_SNPs > length(SS_dat$SNP)){ SS_dat$SNP}else{SS_dat$SNP[1:max_SNPs]}
  # Get top SNPs from susieR fine mapping
  susieR_dat <- subset(susieDF, PIP!=0) %>%
    arrange(desc(PIP)) %>%  mutate(SNP=as.character(SNP))
  susieR_SNPs <-if(max_SNPs > length(susieR_dat$SNP)){ susieR_dat$SNP}else{susieR_dat$SNP[1:max_SNPs]}
  # Calculate percent
  overlap <-  length(intersect(SumStats_SNPs, susieR_SNPs))
  percentOverlap <- overlap / length(susieR_SNPs) * 100
  cat("\n",overlap," / ",length(susieR_SNPs), " (",round(percentOverlap,2),
      "%) of SNPs of the SNPs in the summary stats were confirmed after fine-mapping.","\n",sep="")
}
# before_after_consensus(top_SNPs, susieDF, max_SNPs=10)


# Fine Map Iteratively

#Nalls_SS %>% group_by(`Nearest Gene`) %>% tally() %>% subset(n>2)
finemap_geneList <- function(top_SNPs, geneList, file_path,
                             bp_distance=500000, num_causal=1,
                             chrom_col="CHR", position_col="POS", snp_col="SNP",
                             pval_col="P", effect_col="Effect", stderr_col="StdErr",
                             LD_reference="1KG_Phase1", superpopulation="EUR",
                             topVariants=3,vcf_folder=F){

  fineMapped_topSNPs <- data.table()
  for (gene in geneList){
    cat('\n')
    cat("###", gene, "\n")
    susieDF <- susie_on_gene(gene=gene, top_SNPs=top_SNPs, num_causal = 1,
                             file_path=file_path, bp_distance=bp_distance,
                             chrom_col=chrom_col, position_col=position_col, snp_col=snp_col,
                             pval_col=pval_col, effect_col=effect_col, stderr_col=stderr_col,
                             LD_reference=LD_reference, superpopulation=superpopulation)
    try({
      before_after_plots(gene, susieDF)
    })
    try({
      before_after_consensus(gene, top_SNPs, susieDF, max_SNPs=10)
    })
    # Create summary table for all genes
    newEntry <- cbind(data.table(Gene=gene), subset(susieDF, PIP==max(PIP)) %>% as.data.table())
    fineMapped_topSNPs <- rbind(fineMapped_topSNPs, newEntry)
    cat('\n')
    createDT_html(susieDF)
  }
  createDT_html(fineMapped_topSNPs, "Potential Causal SNPs Identified by susieR", scrollY = 200)
  return(susieDF)
}



