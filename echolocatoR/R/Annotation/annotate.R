# %%%%%%%%%%%%%%%%% #
####### Annotate ####### 
# %%%%%%%%%%%%%%%%% # 


lead.SNP.coords <- function(){
  annot <- readxl::read_excel("./Data/annotated_finemapping_results.xlsx")
  annot <- find_consensus_SNPs(annot, support_thresh = 2)
  annot[is.na(annot$mean.PP),"mean.PP"] <-0
  
  annot.sub <- subset(annot, leadSNP==T, select=c(Gene, SNP, CHR, POS)) %>% 
    dplyr::rename(lead.SNP=SNP) %>% 
    dplyr::mutate(min.POS=POS - 1e+06, max.POS=POS + 1e+06)
  data.table::fwrite(annot.sub, "./Data/lead.SNP.coords.csv", sep=",")
  annot[is.na(annot)] <-0
  
  SNPgroup.summary <- function(DF, group.name=""){
    n.total = nrow(DF)
    n.per.locus <- (DF %>% group_by(Gene) %>% tally())$n %>% mean() 
    means <- DF[,c("P","Effect","mean.PP","MAF")] %>% dplyr::summarise_all(function(x){abs(mean(x, na.rm = T))})
    means <- cbind(`SNP Group`=group.name, means, `SNPs / locus` = n.per.locus, `Total SNPs` = n.total) 
    # Count distribution
    mean.size.statCS <- table((statCS %>% tally())$n )
    print(mean.size.statCS)
    return(means)
  }
  
  # All SNPs 
  allSNPs.means <- SNPgroup.summary(annot, group.name="All GWAS")
  # Nom sig GWAS SNPs
  nomSigSNPs <- subset(annot, P<0.05) 
  nomSigSNPs.means <- SNPgroup.summary(nomSigSNPs, group.name="nom. sig. GWAS")
  # FDR sig GWAS SNPs
  fdrSigSNPs <- subset(annot, P<5e-8) 
  fdrSigSNPs.means <- SNPgroup.summary(fdrSigSNPs, group.name="FDR sig. GWAS")
  
  #CS stats
  
  # Lead GWAS SNPs
  leadSNPs <- subset(annot, leadSNP) 
  leadSNP.means <- SNPgroup.summary(leadSNPs, group.name="Lead GWAS SNP")
  #CS stats
  ## Statistical FM
  statCS <- annot %>% 
    group_by(Gene) %>% 
    subset(SUSIE.Credible_Set>0 | FINEMAP.Credible_Set > 0| PAINTOR.Credible_Set >0)  
  statCS.means <- SNPgroup.summary(statCS, group.name="Statistical.CS")
  # Functional FM
  funcCS <- annot %>% group_by(Gene) %>% subset(PAINTOR.Credible_Set>0)
  funcCS.means <- SNPgroup.summary(funcCS, group.name="Functional.CS")
  # Consensus stats
  consensus <- annot %>% group_by(Gene) %>% subset(Consensus_SNP)
  consensus.means <- SNPgroup.summary(consensus, group.name="Consensus")
  percent.loci.w.consensus <- length(unique(consensus$Gene)) / length(unique(annot$Gene))
  
  # Merged
  merged.means <- rbind(allSNPs.means, nomSigSNPs.means, fdrSigSNPs.means, leadSNP.means,statCS.means, funcCS.means, consensus.means)
  merged.means[,-c(1:2)] <-round(merged.means[,-c(1:2)],3)
  merged.means$P <- formatC(merged.means$P , format = "e", digits = 3)
  data.table::fwrite(merged.means, "./Data/GWAS/Nalls23andMe_2019/_genome_wide/SNP.summary.csv")
  
  library(patchwork)
  bins=150
  x.var="Effect"
  alpha=1
  ggplot() + 
    geom_histogram(data=leadSNPs, aes(x=eval(parse(text=x.var)), fill="Lead GWAS SNPs"), alpha=alpha,   fill="red", bins = bins) +  
    labs(title = "Lead GWAS SNPs", x=x.var) +
    theme_classic() + 
    xlim(min(annot[[x.var]]), max(annot[[x.var]])) +
    
    ggplot() + 
    geom_histogram(data = statCS, aes(x=eval(parse(text=x.var)), fill="Statistical Credible Set"), alpha=alpha,  fill="green", bins = bins) + 
    labs(title = "Statistical Credible Set", x=x.var) +
    theme_classic() + 
    xlim(min(annot[[x.var]]), max(annot[[x.var]])) +
    
    ggplot() + 
    geom_histogram(data = funcCS, aes(x=eval(parse(text=x.var)), fill="Functional Credible Set"), alpha=alpha,   fill="green4", bins = bins) + 
    labs(title = "Functional Credible Set", x=x.var) +
    theme_classic() + 
    xlim(min(annot[[x.var]]), max(annot[[x.var]])) +
    
    ggplot() + 
    geom_histogram(data = consensus, aes(x=eval(parse(text=x.var)), fill="Consensus SNPs"), alpha=alpha,  fill="goldenrod2", bins = bins) + 
    labs(title = "Consensus SNPs", x=x.var) +
    theme_classic() + 
    xlim(min(annot[[x.var]]), max(annot[[x.var]])) +
    
    patchwork::plot_layout(ncol = 1)
  
  return(merged.means)
}





compare_finemapping_methods <- function(dataset="./Data/GWAS/Nalls23andMe_2019"){
  FM_orig <- merge_finemapping_results(minimum_support = 0, 
                                       dataset = dataset, 
                                       exclude_methods = NULL)   
  # counts <- (FM_orig %>% group_by(Gene) %>% count()) 
  # plotrix::std.error(counts$n)
  # min(counts$n)
  # max(counts$n)
  
  # Remove loci that were manually added
  ## Also remove TRIM40 for now bc I'm having issues getting its LD
  FM_all <- subset(FM_orig, !(Gene %in%c("ATG14","SP1","LMNB1","ATP6V0A1","TRIM40")) )  
  # Identify any un-finemapped loci
  FM_all %>% group_by(Gene) %>% summarise(leadGWAS=sum(leadSNP==T)) %>% arrange(leadGWAS)
  
  # Proportion of CS SNPs that are the leadSNP (by tool)
  dat <- FM_all %>%  
    summarise_at(vars(ends_with(".Credible_Set")), 
                          funs(CS.sum=sum(.>0, na.rm = T),
                               leadGWAS.sum=sum(leadSNP==T, na.rm = T),
                               leadGWAS.CS.sum=sum(leadSNP==T & .>0, na.rm = T),
                               PROP.SNPS=sum(leadSNP==T & .>0, na.rm = T)/sum(.>0, na.rm = T)) )
  colnames(dat) <- gsub("\\.Credible_Set","", colnames(dat))
  dat
  
  # Proportion of loci in which the CS contains the leadSNP (by tool)
  # dat3.1 <- FM_all %>% group_by(Gene) %>% 
  #   summarise_at(vars(ends_with(".Credible_Set")), 
  #                funs(leadGWAS.CS.size=sum(leadSNP==T & .>0, na.rm = T)))
  # dat3.2 <- FM_all %>% group_by(Gene) %>% 
  #   summarise_at(vars(ends_with(".Credible_Set")), 
  #                funs(leadGWAS.CS.size=sum(.>0, na.rm = T)))
  # colSums(dat3.1[,-1]) / colSums(dat3.2[,-1])
  
  # Number of CS SNPs (by tool)
  dat2 <- FM_all %>% group_by(Gene) %>% summarise_at(vars(ends_with(".Credible_Set")), 
                          funs(CS=sum(.,na.rm=T)))
  colnames(dat2) <- gsub("\\.Credible_Set","", colnames(dat2))
  # Proportion of loci with at least one CS SNP (by tool)
  colSums(dat2[,-1]>0) / nrow(dat2) 
  
  
   
  
  # subset(dat, ABF._CS.size>0) %>%
  
  # dat %>% dplyr::select(ends_with("_SUM")) 
  # 
  # dat %>%  summarise_at(vars(ends_with("_SUM")), 
  #                       funs(PROP_LOCI=sum(.,na.rm = T)/ sum(.>0) ) )
  #  
}


plot_snpGroupPP_by_tool <- function(FM_all){ 
  
  dat <- FM_all %>% summarise_at(vars(ends_with(".PP")), 
                            funs(Overall=mean(.,na.rm=T), 
                                 GWAS.nom.sig=mean(subset(.,P<.05),na.rm=T),
                                 GWAS.sig=mean(subset(.,P<5e-8),na.rm=T), 
                                 GWAS.lead=mean(subset(.,leadSNP==T),na.rm=T),
                                 Credible.Set=mean(subset(.,Support>0),na.rm=T),
                                 Consensus=mean(subset(.,Consensus_SNP==T),na.rm=T)
                                 )) 
  pdat <- data.frame(PP=t(dat))
  pdat$Method <- gsub("\\.PP_.*", "", rownames(pdat))
  pdat$SNP.Group <- gsub(".*\\.PP_", "", rownames(pdat))
  pdat$SNP.Group <- factor(pdat$SNP.Group, levels = unique(pdat$SNP.Group), ordered = T)
  
  ggplot(pdat, aes(x=SNP.Group, y=PP, fill=Method)) + 
    geom_col(position = "dodge") + 
    labs(y="mean PP") +
    theme_bw() 
}


top_finemapped_loci <- function(dataset="./Data/GWAS/Nalls23andMe_2019",
                                save_results=T, 
                                biomart=T){ 
  FM_orig <- merge_finemapping_results(minimum_support = 0, 
                                      dataset = dataset, 
                                      exclude_methods = NULL)   
  FM_all <- FM_orig 
  # missing.loci <- unique(subset(FM_all, !is.na(FINEMAP.Probability))$Gene)
  # FM_all <- subset(FM_all, !(Gene %in% missing.loci))  
    # dplyr::group_by(FM_all, Gene) %>% 
    # summarise(n_consensus = sum(Consensus_SNP == T), 
    #           n_CS = sum(Support)) %>% 
    #   subset(n_consensus>0) %>% arrange(n_consensus, n_CS)  
  FM_all <- subset(FM_all, !(Gene %in%c("ATG14","SP1","LMNB1","ATP6V0A1")) )
  
  # List all available QTL datasets
  list_Data_dirs() %>% dplyr::filter(grepl("QTL",type)) %>% dplyr::select(Dataset, type)
  qlt.list <- c("psychENCODE_eQTL",
                "Fairfax_2014_CD14",
                "Fairfax_2014_IFN",
                "Fairfax_2014_LPS2",
                "Fairfax_2014_LPS24",
                "Cardiogenics_macrophages",
                "Cardiogenics_monocytes",
                "MESA_CAU"
                # "Brain.xQTL.Serve_eQTL",
                # "Brain.xQTL.Serve_haQTL",
                # "Brain.xQTL.Serve_mQTL"
                ) 
  FM_tmp <- FM_all
  for(qtl in qlt.list){
    FM_tmp <- mergeQTL.merge_handler(FM_all = FM_tmp,  
                                     qtl_file = qtl, 
                                     force_new_subset = F)
    FM_tmp$QTL.sig <- ifelse(data.frame(FM_tmp)[,"QTL.P"]<=5e-8,"Y","N") 
    colnames(FM_tmp) <- gsub("^QTL\\.",paste0(qtl,"."), colnames(FM_tmp))  
  } 
  FM_all <- FM_tmp
  qtl.sig.cols <- grep("\\.sig$",colnames(FM_all),value = T)
  FM_all$QTL.count <- rowSums(data.frame(FM_all)[,qtl.sig.cols]=="Y", na.rm = T)
  
  # FM_all %>% dplyr::mutate()
  # Biomart Annotations
  if(biomart){
    query_snps <- unique(subset(FM_all, Consensus_SNP)$SNP)
    printer("+ BIOMART:: Gathering annotations for",length(query_snps),"SNPs...")
    SNP.info <- biomart_snp_info(snp_list = query_snps) 
    SNP.info.collapse <- SNP.info %>% 
      dplyr::rename(SNP=refsnp_id) %>% 
      dplyr::select(SNP, consequence_type_tv, reg_consequence_types) %>% 
      dplyr::group_by(SNP) %>% 
      dplyr::summarise(consequence_type_tv= paste0(unique(consequence_type_tv), collapse = "/"),
                       reg_consequence_types= paste0(unique(reg_consequence_types), collapse = "/") ) %>% 
      dplyr::mutate(consequence_type_tv=gsub(", ,|NA, ", "",consequence_type_tv),
                    reg_consequence_types=gsub(", ,|NA, ", "",reg_consequence_types)) 
    FM_all <- data.table:::merge.data.table(FM_all, 
                                            SNP.info.collapse, 
                                            by = "SNP")
  }
  
  # Create summary data.frame 
  library(tidyverse) 
  FM_all$GWAS.lead <- ifelse(FM_all$leadSNP==T,"Y","N")
  grouped.dat <- FM_all %>% group_by(Gene) 
  cols <- list( 
    ## SNP Group counts
    grouped.dat %>% 
      summarise(Consensus.RSID=paste(SNP[Consensus_SNP==T], collapse=", "),
                Consensus.ID=paste(SNP_id[Consensus_SNP==T], collapse=", "),
                CredSet.RSID=paste(SNP[Support>0], collapse=", "),
                Total.size=n(), 
                GWAS.nom.sig.size=sum(P<0.05),
                GWAS.sig.size=sum(P<5e-8), 
                CredSet.size=sum(Support>0), 
                Consensus.size=sum(Consensus_SNP==T)) %>%
      dplyr::select(-Gene),
    ## Text cols
    ### Is the Consensus SNP the GWAS lead?
    grouped.dat %>%
      summarise_at(vars(GWAS.lead, ends_with(".sig")), 
                     funs(Consensus=paste(replace_na(subset(., Consensus_SNP), "N"),collapse=", "),
                          CredSet=paste(replace_na(subset(., Support>0), "N"),collapse=", ")) ),
     
    # Numeric cols
    ## As separated text
    grouped.dat %>%  
        summarise_at(vars(ends_with(".PP"), ends_with("QTL.count"),ends_with("Effect"), MAF),
                     funs(paste(round(subset(., Consensus_SNP),3),collapse=", ")) ) %>%
      dplyr::select(-Gene),
    # As means 
    grouped.dat %>%
      summarise_at(vars(mean.PP, ends_with("QTL.count")),
                   funs(avg=mean(subset(., Consensus_SNP), na.rm=T)) ) %>%
      dplyr::select(-Gene)  
  )
  top_loci <- do.call(cbind, cols)
    
  # Check whether the locus is novel according to the most recent Nalls et al (2019) PD GWAS
  ## `Known GWAS locus within 1MB` (locus-level)***
  ## `Locus within 250KB` (SNP-level?)
  Nalls <- readxl::read_excel("./Data/GWAS/Nalls23andMe_2019/Nalls2019_TableS2.xlsx")
  Nalls.novel <- (Nalls %>% dplyr::mutate(Novel.Locus=ifelse(`Known GWAS locus within 1MB`==1,"N","Y")))[c("Nearest Gene","Novel.Locus") ]%>% dplyr::rename(Gene=`Nearest Gene`) %>% unique()
  # Nalls.novel <-Nalls.novel %>% 
  #   group_by(Gene) %>% 
  #   summarise_each(funs(paste(., collapse = ", ")))
  top_loci <- data.table:::merge.data.table(data.table::data.table(top_loci), 
                                            data.table::data.table(Nalls.novel), 
                                            by="Gene") %>% dplyr::rename(Locus=Gene)
  top_loci <- top_loci %>% 
    dplyr::mutate(GWAS.lead_Consensus.any = grepl("Y",GWAS.lead_Consensus) ,
                  GWAS.lead_CredSet.any = grepl("Y",GWAS.lead_CredSet))
 
  # Sort/filter by criterion
  top_loci_filt <- top_loci %>% 
    ## [0] It's one of the PD GWAS loci
    subset(GWAS.sig.size>0) %>%
    ## [1] There is at least one consensus SNP
    subset(Consensus.size==1) %>% 
    ## [2] None of the consensus SNPs are the lead GWAS SNP
    # dplyr::filter(grepl("N",GWAS.lead) & !grepl("Y",GWAS.lead)) %>%
    subset(GWAS.lead_Consensus.any==F) %>%
    ## [3] There's at least one QTL
    subset(QTL.count_avg>0) %>%
    ## [3] Just one consensus SNP and small Credible Set
    arrange(Consensus.size, desc(QTL.count), CredSet.size)  
  top_loci_filt <- rbind(top_loci_filt, subset(top_loci, Locus %in% c("LRRK2")))
  
  top_loci_sort <- subset(top_loci, GWAS.sig.size>0 & Consensus.size>0) %>% 
    arrange(Consensus.size, GWAS.lead_Consensus.any, desc(QTL.count_avg))
  
  if(save_results){
    topLoci.path <- file.path(dataset,"_genome_wide/top_loci.csv")
    printer("+ Saving top loci ==>",topLoci.path)
    data.table::fwrite(top_loci, topLoci.path)
  } 
  
  return(top_loci) 
}




# multi_finemap_results_table <- function(multi_finemap_DT, 
#                                         finemap_method_list, 
#                                         fancy_table=F, 
#                                         minimum_support=0,
#                                         include_leadSNPs=T){ 
#   # finemap_DT <- data.table::fread(file.path(results_path,"Multi-finemap/Multi-finemap_results.txt"),stringsAsFactors = F)
#   CS_cols <- colnames(multi_finemap_DT)[endsWith(colnames(multi_finemap_DT), ".Credible_Set")]
#   if(include_leadSNPs){
#     support_DT <- subset(multi_finemap_DT, Support >= minimum_support | leadSNP==T)
#   } else {
#     support_DT <- subset(multi_finemap_DT, Support >= minimum_support)
#   } 
#   # support_DT <- subset(support_DT, select=c("Gene","SNP","CHR","POS","P","leadSNP",CS_cols,"Support"))  %>%
#   #   arrange(desc(Support))
#   support_DT <- dplyr::select(support_DT, -dplyr::one_of(c("Dataset"))) %>%
#       arrange(desc(Support))
#   # Plot table 
#   if(fancy_table){
#     customGreen0 = "#DeF7E9" 
#     customGreen = "#71CA97" 
#     customRed = "#ff7f7f" 
#     CS_formatter <- 
#       formattable::formatter("span", 
#                              style = x ~ style( 
#                                color = ifelse(x > 0, customGreen, ifelse(x == 0, "black", "black")))) 
#     formattable::formattable(support_DT, 
#                              align =c("l","c","c","c","c", "c", "c", "c", "r"), 
#                              list( P = formattable::color_tile(customGreen, customGreen0),
#                                    SUSIE.Credible_Set = CS_formatter,
#                                    ABF.Credible_Set = CS_formatter,
#                                    FINEMAP.Credible_Set = CS_formatter,
#                                    COJO.Credible_Set = CS_formatter,
#                                    Support = formattable::color_tile("white", "green")) 
#     )
#     
#   }
#   return(support_DT)
# }


leadSNP_comparison <- function(top_SNPs, merged_results){
  leadSNP_summary_table <- data.table:::merge.data.table(
    top_SNPs %>% dplyr::select(leadSNP=SNP, Gene),
    merged_results %>% dplyr::select(finemappedSNP=SNP, Gene),
    by="Gene", all=T) %>% 
    dplyr::group_by(Gene, leadSNP) %>% 
    dplyr::mutate(Overlap = leadSNP %in% finemappedSNP) %>% 
    dplyr::group_by(Gene) %>% 
    dplyr::summarise(Overlap = sum(Overlap)) %>% dplyr::mutate(leadSNP_in_CredSet = Overlap > 0 )
  
  percent_leadSNPs <- round( sum(leadSNP_summary_table$leadSNP_in_CredSet) /
                               length(leadSNP_summary_table$leadSNP_in_CredSet) * 100,2) 
  return(leadSNP_summary_table[,c("Gene","leadSNP_in_CredSet")])
}


merge_finemapping_results <- function(minimum_support=0, 
                                      include_leadSNPs=T,
                                      xlsx_path="./Data/annotated_finemapping_results.xlsx",
                                      from_storage=T,
                                      haploreg_annotation=F,
                                      regulomeDB_annotation=F,
                                      biomart_annotation=F,
                                      verbose=T,
                                      dataset="./Data/GWAS",
                                      PP_threshold=.95, 
                                      consensus_thresh=2, 
                                      exclude_methods=NULL){ 
  if(from_storage){
    printer("+ Gathering all fine-mapping results from storage...", v=verbose)
    # Find all multi-finemap_results files
    multi_dirs <- list.files(dataset, pattern = "Multi-finemap_results.txt", 
                             recursive = T, full.names = T)
    dataset_names <- dirname(dirname(dirname(multi_dirs))) %>% unique() 
    # Loop through each GENE
    finemap_results <- lapply(dataset_names, function(dn, multi_dirs.=multi_dirs){
      gene_dirs <- dirname(dirname(multi_dirs.))
      # Loop through each gene folder
      all_results <- lapply(gene_dirs, function(gd){
        gene <- basename(gd)
        dirname(gd) 
        printer("+ Importing results...",gene, v=verbose)
        multi_data <- data.table::fread(file.path(gd,"Multi-finemap/Multi-finemap_results.txt"), 
                                        nThread = 4)
        multi_data <- cbind(data.table::data.table(Dataset=dn, Gene=gene), multi_data)
        return(multi_data)
      }) %>% data.table::rbindlist(fill=TRUE) # Bind genes
    }) %>% data.table::rbindlist(fill=TRUE) # Bind datasets    
  }
  
  
  # Add/Update Support/Consensus cols 
  merged_results <- find_consensus_SNPs(finemap_DT = finemap_results, 
                                        credset_thresh = PP_threshold, 
                                        consensus_thresh = consensus_thresh, 
                                        exclude_methods = exclude_methods)
  
  # Loop through each DATASET
  # merged_results <- lapply(unique(finemap_results$Dataset), function(dname, include_leadSNPs.=include_leadSNPs){ 
  #   multi_finemap_DT <- subset(finemap_results, Dataset==dname) 
  #   CS_cols <- colnames(multi_finemap_DT)[endsWith(colnames(multi_finemap_DT), ".Credible_Set")]
  #   finemap_method_list <- gsub(".Credible_Set","",CS_cols) 
  #   # Create support table
  #   support_DT <- multi_finemap_results_table(multi_finemap_DT,
  #                                             finemap_method_list,
  #                                             fancy_table = F,
  #                                             minimum_support = minimum_support,
  #                                             include_leadSNPs = include_leadSNPs.)
  #   support_DT <- cbind(Dataset = dname, support_DT) %>% arrange(Gene, desc(Support)) 
  # }) %>% data.table::rbindlist() 
  
  
  
  
  # Annotate with haplorR
  if(haploreg_annotation){
    HR_query <- haploR.HaploReg(snp_list = unique(merged_results$SNP), verbose = verbose)
    merged_results <- data.table:::merge.data.table(merged_results, HR_query, 
                                                    by.x = "SNP", 
                                                    by.y = "rsID", 
                                                    all = T,
                                                    allow.cartesian=TRUE)
  }
  if(regulomeDB_annotation){
    regDB_query <- haploR.regulomeDB(snp_list = unique(merged_results$SNP), verbose = verbose)
    merged_results <- data.table:::merge.data.table(merged_results, regDB_query, 
                                                    by.x = "SNP", 
                                                    by.y = "rsID", 
                                                    all = T,
                                                    allow.cartesian=TRUE)
  }
  if(biomart_annotation){
    biomart_query <- biomart_snp_info(snp_list = merged_results$SNP, verbose = verbose) 
    merged_results <- data.table:::merge.data.table(merged_results, biomart_query, 
                                                      by.x = "SNP",
                                                      by.y = "refsnp_id", 
                                                      all = T, 
                                                      allow.cartesian=TRUE)  
  }
  
  if(xlsx_path!=F){
    # data.table::fwrite(merged_results, file = csv_path, quote = F, sep = ",")
    openxlsx::write.xlsx(merged_results, xlsx_path)
  } 
  # createDT_html(merged_results) %>% print() 
  return(merged_results)
}


counts_summary <- function(top_SNPs, merged_results, verbose=T){
  # Get total # of SNPs per gene per dataset  
  candidate_counts <- merged_results %>% dplyr::group_by(Dataset, Gene) %>% 
    count(name = "Total_SNPs")
  max_consensus <- sum(endsWith(colnames(merged_results),".Credible_Set"))
  candidate_counts <- merge.data.frame(candidate_counts, 
                                       dplyr::group_by(merged_results, Gene) %>% 
                                         count(name = "Credible_Set"), 
                                       by = "Gene", all = T)
  candidate_counts <- merge.data.frame(candidate_counts,  
                                       subset(merged_results, Support == max_consensus) %>%
                                         dplyr::group_by(Gene) %>% 
                                         count(name = "Consensus_SNPs"),
                                       all = T)
  # Add lead rsid column
  candidate_counts <-  merge.data.frame(candidate_counts,  
                                        top_SNPs[,c("Gene","SNP")] %>% dplyr::rename(leadSNP=SNP), 
                                        on = "Gene", all = T) 
  # Add "is lead SNP in Credible Set" column
  candidate_counts <-  merge.data.frame(candidate_counts,  
                                        leadSNP_comparison(top_SNPs, merged_results), 
                                        on = "Gene", all = T) 
  # Gather credible set rsids
  CredSet_rsids <- merged_results %>% dplyr::group_by(Dataset, Gene) %>%
    subset(Support==max_consensus) %>%
    dplyr::summarise(CredSet_rsids = paste(SNP,collapse="; ")) %>%
    data.table::data.table()
  candidate_counts <- merge.data.frame(candidate_counts,  
                                       CredSet_rsids, 
                                       on = "Gene", all = T) 
  # Gather consensus rsids
  consensus_rsids <- merged_results %>% dplyr::group_by(Dataset, Gene) %>%
    subset(Support==T) %>%
    dplyr::summarise(Consensus_rsids = paste(SNP,collapse="; ")) %>%
    data.table::data.table()
  candidate_counts <- merge.data.frame(candidate_counts,  
                                       consensus_rsids, 
                                       on = "Gene", all = T) 
  # Fill 0s
  candidate_counts$Consensus_SNPs[is.na(candidate_counts$Consensus_SNPs)] <- 0 
  means <- c(Gene=" ",
             Dataset="[Average]",
             candidate_counts[,c("Total_SNPs","Credible_Set","Consensus_SNPs")] %>% colMeans() %>% round(1),
             leadSNP_in_CredSet = paste0(round(sum(candidate_counts$leadSNP_in_CredSet) / dim(candidate_counts)[1]*100,2),"%"),
             CredSet_rsids = "",
             Conensus_rsids = ""
  )
  # Add averages 
  candidate_counts <- suppressWarnings(rbind(candidate_counts,  means))
  percent_model_convergence <- round(sum(candidate_counts$Consensus_SNPs>0)  / length(candidate_counts$Consensus_SNPs) *100, 2)
  max_consensus_set <- max(candidate_counts$Consensus_SNPs)
  # Check if lead SNP is in the credible sets for each locus
  printer("\n + All",max_consensus,"models converged upon 1 to",
          max_consensus_set,"SNPs in",percent_model_convergence,"% of loci.", 
          v=verbose)
  # createDT_html(candidate_counts) %>% print()
  return(candidate_counts)
}









# HaploR
# https://cran.r-project.org/web/packages/haploR/vignettes/haplor-vignette.html
haploR.HaploReg <- function(snp_list, verbose=T, chunk_size=NA){
  printer("+ Gathering annotation data from HaploReg...", v=verbose)
  # Break into smaller chunks
  snp_list <- unique(snp_list)
  if(is.na(chunk_size)){chunk_size <- length(snp_list)}
  chunked_list <- split(snp_list, ceiling(seq_along(snp_list)/chunk_size)) 
  
  HR_query <- lapply(names(chunked_list), function(i){ 
    printer("++ Submitting chunk",i,"/",length(chunked_list)) 
    chunk <- chunked_list[[i]]
    HR_query <-  haploR::queryHaploreg(query = chunk, file = NULL, study = NULL, ldThresh = NA,
                                       ldPop = "EUR", epi = "vanilla", cons = "siphy", genetypes = "gencode",
                                       url = "https://pubs.broadinstitute.org/mammals/haploreg/haploreg.php",
                                       timeout = 500, encoding = "UTF-8", verbose = FALSE) 
    
    return(data.table::as.data.table(HR_query))
  }) %>% data.table::rbindlist()
  
  return(HR_query)
}

haploR.regulomeDB <- function(snp_list, verbose=T, chunk_size=NA){
  printer("+ Gathering annotation data from HaploReg...", v=verbose)
  # Break into smaller chunks
  snp_list <- unique(snp_list)
  if(is.na(chunk_size)){chunk_size <- length(snp_list)}
  chunked_list <- split(snp_list, ceiling(seq_along(snp_list)/chunk_size)) 
  
  rDB_query <- lapply(names(chunked_list), function(i){ 
    printer("++ Submitting chunk",i,"/",length(chunked_list)) 
    chunk <- chunked_list[[i]]
    rdb_query <-  haploR::queryRegulome(query = chunk, 
                                       timeout = 500, 
                                       verbose = F)  
    return(data.table::as.data.table(rdb_query))
  }) %>% data.table::rbindlist()
  
  return(rDB_query)
}


biomart_snp_info <- function(snp_list, reference_genome="grch37", verbose=T){
  printer("+ Gathering annotation data from Biomart...", v=verbose) 
  mart = biomaRt::useMart(biomart="ENSEMBL_MART_SNP", 
                 host=paste0(reference_genome,".ensembl.org"), 
                 path="/biomart/martservice", 
                 dataset="hsapiens_snp")
  # View(biomaRt::listFilters(mart))
  # View(biomaRt::listAttributes(mart))
  biomart_query = biomaRt::getBM(attributes =  c('refsnp_id',
                                 'allele',
                                 'chr_name',
                                 'chrom_start',
                                 'chrom_end',
                                 'chrom_strand',
                                 'ensembl_gene_stable_id',
                                 'consequence_type_tv',
                                 'polyphen_prediction',
                                 'polyphen_score',
                                 'sift_prediction',
                                 'sift_score',
                                 'reg_consequence_types',
                                 'validated'
                                 ),
                  filters = c('snp_filter'),
                  values = unique(snp_list),
                  mart = mart) 
  biomart_query <- data.table::as.data.table(biomart_query)
  biomart_query[biomart_query$consequence_type_tv=="",]$consequence_type_tv <- NA
  # Only take the first annotation per variant
  # annotated_results %>% dplyr::group_by(Dataset, Gene, SNP) %>% slice(1)
  return(biomart_query)
}

# BIOMART
biomart_snps_to_geneInfo <- function(snp_list, reference_genome="grch37"){
  # listMarts()
  snp_mart = useMart("ENSEMBL_MART_SNP", 
                     dataset="hsapiens_snp", 
                     host =  paste0(reference_genome,".ensembl.org"))
  # View(listFilters(snp_mart))
  # View(listAttributes(snp_mart))
  snp_results <- biomaRt::getBM(snp_mart, filters="snp_filter", 
                                values=snp_list,
                                attributes=c("refsnp_id","snp","chr_name", "chrom_start","chrom_end",
                                             "associated_gene","ensembl_gene_stable_id" ) )
  # # Split ensembl IDs
  gene_mart = useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
  gene_results <- biomaRt::getBM(mart = gene_mart, 
                                 filters = "ensembl_gene_id",
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

biomart_geneInfo <- function(geneList, reference_genome="grch37"){
  # listDatasets(useMart("ENSEMBL_MART_ENSEMBL") )
  gene_mart = biomaRt::useMart("ENSEMBL_MART_ENSEMBL", 
                               dataset="hsapiens_gene_ensembl",  
                               host = paste0(reference_genome,".ensembl.org"))
  # View(listFilters(gene_mart))
  # View(listAttributes(gene_mart))
  gene_results <- biomaRt::getBM(mart = gene_mart, 
                                 filters = "hgnc_symbol",
                                 # values = unlist(strsplit(snp_results$ensembl, ";")),
                                 values = geneList,
                                 attributes = c("hgnc_symbol","external_gene_name","ensembl_gene_id",
                                                "chromosome_name", "start_position", "end_position") )
  return(gene_results)
}
# biomart_geneInfo(c("PTK2B","CLU","APOE"))

SNPs_by_mutation_type <- function(merged_results, 
                                  mutation_type="missense_variant"){
  potential_missense <- subset(merged_results, consequence_type_tv == mutation_type) %>% 
    dplyr::group_by(Dataset, Gene, SNP) %>% 
    dplyr::select(Dataset, Gene, SNP, consequence_type_tv) %>% 
    unique()
  potential_missense_full <- subset(merged_results, SNP %in% potential_missense$SNP) %>% 
    dplyr::group_by(Dataset, Gene, SNP) %>% 
    dplyr::select(Dataset, Gene, SNP, consequence_type_tv) %>% 
    unique()
  return(potential_missense_full)
}



epigenetics_summary <- function(merged_results, 
                                tissue_list = c("BRN","BLD"),
                                epigenetic_variables = c("Promoter_histone_marks","Enhancer_histone_marks") # Chromatin_Marks 
                                ){  
  merged_results <- data.table(merged_results)
  summary_func <- function(ev){
    boolean <- lapply(ev, function(x){ intersect(tissue_list, strsplit(as.character(x), ", ")[[1]]) %>% length() > 0 }) %>% unlist()  
    n_hits <- dim(merged_results[boolean,])[1]
    Total <- dim(merged_results)[1]
    Percent_Total <- round(n_hits / Total*100,2)
    return(list(Hits = n_hits,
                Total = Total,
                Percent_Total = Percent_Total))
  } 
  epi_summary <- merged_results[, lapply(.SD, summary_func), .SDcols = epigenetic_variables] %>% t()
  colnames(epi_summary) <- c("Hits","Total_SNPs","Percent_Total")
  print(epi_summary)
}
  

epigenetics_enrichment <- function(snp_list1, 
                                   snp_list2, 
                                   chisq=T, 
                                   fisher=T,
                                   epigenetic_variables = c("Promoter_histone_marks","Enhancer_histone_marks"),
                                   tissue_list = c("BRN","BLD")){
  printer("Conducting SNP epigenomic annotation enrichment tests...")
  # Foreground
  printer("+++ SNP list 1 :")
  HR1 <- haploR.HaploReg(snp_list1)
  summ1 <- epigenetics_summary(HR1, tissue_list = tissue_list)
  # Background
  printer("+++ SNP list 2 :")
  HR2 <- haploR.HaploReg(snp_list2)
  summ2 <- epigenetics_summary(HR2, tissue_list = tissue_list)
  
  for(epi_var in epigenetic_variables){
    printer("++ Testing for enrichment of '",epi_var,
            "' in the tissues '",paste(tissue_list, collapse=" & "),"'") 
    # Create contingency table
    cont_tab <- rbind(list1 = summ1[epi_var, c("Hits","Total_SNPs")] %>% unlist, 
                      list2 = summ2[epi_var, c("Hits","Total_SNPs")] %>% unlist() ) %>% as.table()
    # Conduct tests
    if(chisq){
      chisq.results <- chisq.test(cont_tab, simulate.p.value = TRUE) 
      print(chisq.results)
    }
    if(fisher){
      fisher.results <- fisher.test(cont_tab)
      print(fisher.results)
    } 
  }
}

 
 