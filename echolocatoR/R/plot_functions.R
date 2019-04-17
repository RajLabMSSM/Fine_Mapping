# %%%%%%%%%%%%%%%%% #
####### PLOT ####### 
# %%%%%%%%%%%%%%%%% #


## Before-After Plots
before_after_plots <- function(gene, susieDF, before_var="P"){
  roundBreaks <- seq(plyr::round_any(min(susieDF$POS),10000), max(susieDF$POS),250000)
  
  # Label top snps
  ## BEFORE fine-mapping
  if(before_var=="Effect"){
    leadSNP_before <- subset( susieDF %>% arrange(desc(abs(Effect))), leadSNP==T)[1,]
    yLimits1 <- c(min(susieDF[before_var]), max(susieDF[before_var])*1.1 )
    y_var = "Effect"
  }else{
    # Sort by pval and then absolute Effect size
    leadSNP_before <- subset( susieDF %>% arrange(P), leadSNP==T)[1,]
    yLimits1 <- c(0, max(-log10(leadSNP_before[before_var]))*1.1)
    y_var = "-log10(p-value)"
  }
  leadSNP_before$type <- "before"
  leadSNP_before$color <- "red"
  # AFTER fine-mapping
  leadSNP_after = subset(susieDF, credible_set==T) %>% arrange(desc(PIP))
  leadSNP_after$type <- "after"
  leadSNP_after$color <- "darkgreen"
  leadSNP_after[leadSNP_after$PIP==max(leadSNP_after$PIP), "color"] <- "green3"
  
  labelSNPs <- rbind(leadSNP_before, leadSNP_after)
  
  # -log10(eval(parse(text=before_var)))
  before_plot <- ggplot(susieDF, aes(x=POS, y=-log10(eval(parse(text=before_var))), label=SNP, color= -log10(P) )) +
    # ylim(yLimits1) +
    geom_hline(yintercept=0, alpha=.5, linetype=1, size=.5) +
    geom_point(alpha=.5) +
    geom_segment(aes(xend=POS, yend=yLimits1[1], color= -log10(P) ), alpha=.5) +
    geom_point(data=labelSNPs[labelSNPs$type=="before",], pch=21, fill=NA, size=4, colour=labelSNPs[labelSNPs$type=="before","color"], stroke=1) +
    geom_point(data=labelSNPs[labelSNPs$type=="after",][1,], pch=21, fill=NA, size=4, colour=labelSNPs[labelSNPs$type=="after","color"][1], stroke=1) +
    geom_text_repel(data=labelSNPs, aes(label=SNP), color=labelSNPs$color, segment.alpha = .5, nudge_x = .5) +
    labs(title=paste(gene," (",length(susieDF$PIP)," variants)","\nBefore Fine Mapping",sep=""),
         y=y_var, x="Position", color="-log10(p-value)") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(breaks = roundBreaks)
  
  ## After fine-mapping
  yLimits2 <- c(min(susieDF$PIP), max(susieDF$PIP)*1.1)
  
  after_plot <- ggplot(susieDF, aes(x=POS, y=PIP, label=SNP, color= -log10(P) )) +
    # ylim(yLimits) +
    geom_hline(yintercept=0,alpha=.5, linetype=1, size=.5) +
    geom_point(alpha=.5) +
    geom_segment(aes(xend=POS, yend=yLimits2[1], color= -log10(P)), alpha=.5) +
    geom_point(data=labelSNPs[labelSNPs$type=="before",], pch=21, fill=NA, size=4, colour=labelSNPs[labelSNPs$type=="before","color"], stroke=1) +
    geom_point(data=labelSNPs[labelSNPs$type=="after",][1,], pch=21, fill=NA, size=4, colour=labelSNPs[labelSNPs$type=="after","color"][1], stroke=1) +
    geom_text_repel(data=labelSNPs, aes(label=SNP), color=labelSNPs$color, segment.alpha = .5, nudge_x = .5) +
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
