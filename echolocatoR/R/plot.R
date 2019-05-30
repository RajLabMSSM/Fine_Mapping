# %%%%%%%%%%%%%%%%% #
####### PLOT ####### 
# %%%%%%%%%%%%%%%%% # 

 


# SNP_list = c("rs76904798","rs34637584","rs117073808")
COJO_plot <- function(cojo_DT, results_path, conditioned_snps, show_plot=T, save_plot=T){
  gene <- basename(results_path)
  # Label independent SNPs
  ind_SNPs <- subset(cojo_DT, COJO.Credible_Set==T)
  ind_SNPs <- cbind(ind_SNPs, type="Independent",color="purple") 
  labelSNPs <- ind_SNPs#rbind(labelSNPs, ind_SNPs) 
  
  effect_SNPs <- cojo_DT %>% arrange(desc(abs(COJO.Conditioned_Effect)))
  effect_SNPs <- effect_SNPs[1:5,]
  effect_SNPs$type <- "effect_SNPs"
  effect_SNPs$color <- "turquoise3"
  labelSNPs <- rbind(labelSNPs, effect_SNPs)
 
  spacing <- if(length(cojo_DT$SNP)>1000){250000}else{50000}
  roundBreaks <- seq(plyr::round_any(min(cojo_DT$POS),10000), max(cojo_DT$POS),spacing) 
  
  cp <- ggplot(cojo_DT, aes(x=POS, y = -log10(P), label=SNP, color= -log10(P))) +
    # ylim(yLimits1) +
    geom_hline(yintercept=0, alpha=.5, linetype=1, size=.5) +
    geom_point(alpha=.5) +
    geom_segment(aes(xend=POS, yend=0, color= -log10(P) ), alpha=.5) +
    geom_point(data=labelSNPs, pch=21, size=4, colour=labelSNPs$color, stroke=1) + 
    geom_label_repel(data=labelSNPs, aes(label=SNP), color=NA, nudge_x = .5, box.padding = .5,
                     label.size=NA, alpha=.8, seed = 1) +
    geom_label_repel(data=labelSNPs, aes(label=SNP), color=labelSNPs$color, segment.alpha = .5, 
                     nudge_x = .5, box.padding = .5, fill = NA, alpha=1, seed = 1) +
    labs(title=paste("Conditional & Stepwise Results: COJO"),
         subtitle = paste("Conditioned on:",conditioned_snps,"(Independent SNPs Highlighted)"),
         y="-log10(p-value)", x="Position", color="-log10(p-value)") +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5) ) +
    scale_x_continuous(breaks = roundBreaks)
  
  if(save_plot){
    png(filename=file.path(results_path,"COJO/COJO_plot.png"), width = 1000, 800)
    print(cp)
    dev.off()
  }  
  if(show_plot){print(cp)}else{return(cp)} 
}


manhattan_plot <- function(subset_DT, 
                           gene="", title="", subtitle="",
                           SNP_list=c(), alt_color_SNPs=c(), 
                           show_plot=T){
  # Read in data subset
  DT <- subset_DT #data.table::fread(subset_path) %>% data.frame() 
  # Color selected list (assumes first one should be colored differently)
  if(length(SNP_list)==0){
    DT = DT %>% arrange(P, desc(Effect))
    SNP_list = DT$SNP[1:5]
  } 
  labelSNPs = subset(DT, SNP %in% SNP_list)
  labelSNPs$type = "after"
  labelSNPs$color = "red"
  if(length(alt_color_SNPs)>0){ 
    labelSNPs[labelSNPs$SNP %in% alt_color_SNPs, "color"] <- "darkred"
  }
  # write.csv(labelSNPs, "table.csv")
  roundBreaks <- seq(plyr::round_any(min(DT$POS),10000), max(DT$POS),250000)
  yLimits1 <- c(0, max(-log10(max(DT$P)))*1.1)
  title <- if(title==""){
    paste0(gene," : Before fine-mapping\n",length(DT$P)," variants")
  } else {title }
  
  mplot <- ggplot(DT, aes(x=POS, y=-log10(P), label=eval(SNP), color= -log10(P) )) +
    # ylim(yLimits1) +
    geom_hline(yintercept=0, alpha=.5, linetype=1, size=.5) +
    geom_point(alpha=.5) +
    geom_segment(aes(xend=POS, yend=yLimits1[1], color= -log10(P) ), alpha=.5) +
    geom_point(data=labelSNPs, pch=21, fill=NA, size=4, colour=labelSNPs$color, stroke=1) +
    # geom_point(data=labelSNPs[labelSNPs$type=="after",][1,], pch=21, fill=NA, size=4, colour=labelSNPs[labelSNPs$type=="after","color"][1], stroke=1) +
    geom_text_repel(data=labelSNPs, aes(label=SNP), color=labelSNPs$color, segment.alpha = .5, nudge_x = .5) +
    labs(title=title,
         subtitle = subtitle,
         y="-log10(p-value)", x="Position", color="-log10(p-value)") +
    theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5)) +
    scale_x_continuous(breaks = roundBreaks)
  if(show_plot){print(mplot)}else{return(mplot)} 
}






construct_SNPs_labels <- function(DT, lead=T, method=T, consensus=T){ 
  cat("\n + Constructing SNP labels...")
  labelSNPs <- data.table::data.table()
  ## BEFORE fine-mapping  
  if(lead){
    before <- subset( DT %>% arrange(P), leadSNP == T)[1,] 
    y_var = "-log10(p-value)" 
    before$type <- "before"
    before$color <- "red"
    labelSNPs <- rbind(labelSNPs, before)
  }
  if(method){
    # AFTER fine-mapping
    after = subset(DT, Credible_Set>0) %>% arrange(desc(Probability)) 
    if(dim(after)[1]>0){
      after$type <- "after"  
      after$color<- "green3"
      labelSNPs <- rbind(labelSNPs, after) 
    } 
  } 
  if(consensus & "Consensus_SNP" %in% colnames(DT)){
    # Conensus across all fine-mapping tools
    cons_SNPs <- subset(DT, Consensus_SNP==T)
    if(dim(cons_SNPs)[1]>0){
      cons_SNPs$type <- "consensus"
      cons_SNPs$color <- "darkgoldenrod1"
      labelSNPs <- rbind(labelSNPs, cons_SNPs)
    } 
  } 
  return(labelSNPs)
}


snp_plot <- function(finemap_DT,
                     gene,
                     method="original",
                     show_plot=T,
                     subtitle=NA,
                     multi = T){
  score_dict <- setNames(c("-log10(P-value)", "PIP","PP", "PIP", "Conditional Probability","-log10(P-value)"),
                         c("original","SUSIE", "ABF", "FINEMAP", "COJO","COLOC"))
  # X-tick spacing
  spacing <- if(length(finemap_DT$SNP)>1000){250000}else{50000}
  roundBreaks <- seq(plyr::round_any(min(finemap_DT$POS),10000), max(finemap_DT$POS), spacing) 
  
  # yLimits2 <- c(0,1.1)
  
  if(method=="original"){
    DT <- finemap_DT
    p <- ggplot(data = DT, aes(x=POS, y= -log10(P), label=SNP, color= -log10(P) ))
    title <- paste0(gene," : Before fine-mapping")
    labelSNPs <- construct_SNPs_labels(DT, lead=T, method = F, consensus = T)
    subtitle <- if(is.na(subtitle)){paste0(length(DT$SNP)," SNPs")}else{subtitle}
  } else {  
    # COLOC
    if(method=="COLOC"){
      DT <- finemap_DT %>% dplyr::rename(Credible_Set = "Colocalized" )
      title <- paste0(gene," : Colocalization (",method,")") 
      p <- ggplot(data = DT, aes(x=POS, y= -log10(P), label=SNP, color= -log10(P) ))  
    # Fine-mapping methods
    } else if (multi){
      title <- paste0(gene," : After fine-mapping (",method,")")
      DT <- finemap_DT %>% dplyr::rename(Probability = paste0(method,".Probability"),
                                         Credible_Set = paste0(method,".Credible_Set"))
      p <- ggplot(data = DT, aes(x=POS, y=Probability, label=SNP, color= -log10(P) )) + 
        ylim(c(0,1.1)) 
      subtitle <- if(is.na(subtitle)){
        paste0(length(subset(DT, Credible_Set>0)$SNP), " Candidate SNPs")
        }else{subtitle}
    } else {
      title <- paste0(gene," : After fine-mapping (",method,")")
      DT <- finemap_DT
      p <- ggplot(data = DT, aes(x=POS, y=Probability, label=SNP, color= -log10(P) )) + 
        ylim(c(0,1.1)) 
    } 
    labelSNPs <- construct_SNPs_labels(DT, lead=T, method = T, consensus = F)
  } 
  
 

  p <- p + geom_hline(yintercept=0,alpha=.5, linetype=1, size=.5) +
    geom_point(alpha=.5) +
    geom_segment(aes(xend=POS, yend=0, color= -log10(P)), alpha=.5) +
    geom_point(data=labelSNPs, pch=21, fill=NA, size=4, colour=labelSNPs$color, stroke=1) +
    # Labels (one for background, one for text)
    geom_label_repel(data=labelSNPs, aes(label=SNP), color=NA, nudge_x = .5, box.padding = 1,
                     label.size=NA, alpha=.5, seed = 1) +
    geom_label_repel(data=labelSNPs, aes(label=SNP), color=labelSNPs$color, segment.alpha = .5, 
                     nudge_x = .5, box.padding = .5, fill = NA, alpha=1, seed = 1) +
    labs(title=title,
         subtitle=subtitle,
         y=score_dict[[method]], x="Position",
         color="-log10(p-value)") +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
    scale_x_continuous(breaks = roundBreaks) 
  if(show_plot){print(p)} else{ return(p) }
}





multi_finemap_plot <- function(finemap_DT,
                               results_path,
                               finemap_method_list, 
                               gene, 
                               conditioned_snps,
                               original=T, 
                               save_plot=T,
                               ncol=1,
                               width=500,
                               height=1000){ 
  method_list <- if(original){c("original", finemap_method_list)}else{finemap_method_list} 
  
  # Assemble plots in list
  plot_list <- lapply(method_list, function(method){
    cat("\n Plotting...",method)
    if(method=="COJO"){
      p <- COJO_plot(cojo_DT = finemap_DT, 
                results_path = results_path, 
                show_plot = F, 
                save_plot = T,
                conditioned_snps = conditioned_snps)
    } else{
      p <- snp_plot(finemap_DT = finemap_DT, 
                    gene = gene, 
                    method = method, 
                    show_plot = F, 
                    multi = T)
    } 
    return(p) 
  }) 
  # Multi-plot
  ## Adjust plotting parameters
  if(length(finemap_method_list)>3){ncol <- 2; width <- width*2}
  cp <- cowplot::plot_grid(plotlist = plot_list, ncol = ncol)
  # Save plot
  if(save_plot){
    png(filename=file.path(results_path,"Multi-finemap/multi_finemap_plot.png"),width = width, height)
    print(cp)
    dev.off()
  } 
  return(cp)
}




eQTL_barplots <- function(subset_path_list, group_list, SNP_list,
                          x_lab="Population", y_lab="Effect",
                          snp_col="snps", effect_col="beta",
                          title="", writeCSV=F){
  library(ggplot2) 
  dt = data.table::data.table()
  i=1
  for (s in subset_path_list){
    dt_sub <- subset(data.table::fread(s), eval(parse(text = snp_col)) %in% SNP_list)
    dt_sub$Group = group_list[i] #strsplit(strsplit(s,"/")[[1]][4], "_")[[1]][2]
    dt = rbind(dt, dt_sub)
    i=i+1
  } 
  g <- ggplot(data=dt, aes(x=Group, y=eval(parse(text = effect_col)), fill=Group)) +
    geom_col() + facet_grid(~eval(parse(text = snp_col))) + 
    labs(title = title, x = x_lab, y = y_lab) 
  print(g)
  if(writeCSV!=F){
    data.table::fwrite(dt, writeCSV)
  }
} 



# 
# # MESA
# eQTL_barplots(subset_path_list = c("Data/eQTL/MESA/LRRK2_AFA_MESA.txt",
#                                  "Data/eQTL/MESA/LRRK2_CAU_MESA.txt",
#                                  "Data/eQTL/MESA/LRRK2_HIS_MESA.txt"), 
#               group_list = c("AFA","CAU","HIS"),
#               SNP_list = c("rs76904798","rs34637584","rs117073808","rs140722239"), 
#               title = "MESA eQTL Effect Sizes:\nNominated Variants",
#               writeCSV = "Data/eQTL/MESA/eQTL_effect_sizes.csv")
# # Fairfax
# eQTL_barplots(subset_path_list = c("Data/eQTL/Fairfax/LRRK2_Fairfax_CD14.txt",
#                                    "Data/eQTL/Fairfax/LRRK2_Fairfax_IFN.txt",
#                                    "Data/eQTL/Fairfax/LRRK2_Fairfax_LPS2.txt",
#                                    "Data/eQTL/Fairfax/LRRK2_Fairfax_LPS24.txt"),
#               group_list = c("CD14","IFN","LPS2","LP24"), 
#               snp_col="Coord",
#               SNP_list = c("12:40614434","12:40734202","12:40922572") , 
#               x_lab = "Stimulation Condition",
#               title = "Fairfax eQTL Effect Sizes:\nNominated Variants",
#               writeCSV = "Data/eQTL/Fairfax/eQTL_effect_sizes.csv")

# 
# ## Before-After Plots
# before_after_plots <- function(results_path, 
#                                gene, 
#                                finemap_DT, 
#                                before_var="P", 
#                                finemap_method="", 
#                                plot_COJO=F){
#   score_dict <- setNames(c("Posterior Inclusion Probability","Posterior Probability",
#                            "Posterior Probability", "Conditional Probability","Probability"),
#                          c("SUSIE", "ABF", "FINEMAP", "COJO"))
#   roundBreaks <- seq(plyr::round_any(min(finemap_DT$POS),10000), max(finemap_DT$POS),250000)
#   labelSNPs <- label_SNPS(finemap_DT)
#   
#   # -log10(eval(parse(text=before_var)))
#   before_plot <- ggplot(finemap_DT, aes(x=POS, y=-log10(eval(parse(text=before_var))), label=SNP, color= -log10(P) )) +
#     # ylim(yLimits1) +
#     geom_hline(yintercept=0, alpha=.5, linetype=1, size=.5) +
#     geom_point(alpha=.5) +
#     geom_segment(aes(xend=POS, yend=yLimits1[1], color= -log10(P) ), alpha=.5) +
#     geom_point(data=labelSNPs[labelSNPs$type=="before",], pch=21, fill=NA, size=4, colour=labelSNPs[labelSNPs$type=="before","color"], stroke=1) +
#     geom_point(data=labelSNPs[labelSNPs$type=="after",][1,], pch=21, fill=NA, size=4, colour=labelSNPs[labelSNPs$type=="after","color"][1], stroke=1) +
#     geom_text_repel(data=labelSNPs, aes(label=SNP), color=labelSNPs$color, segment.alpha = .5, nudge_x = .5, box.padding = .5) +
#     labs(title="Before Fine-mapping",
#          subtitle = paste(gene," (",length(finemap_DT$Probability)," variants)",sep=""),
#          y=y_var, x="Position", color="-log10(p-value)") +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
#     scale_x_continuous(breaks = roundBreaks)
#   
#   ## After fine-mapping
#   yLimits2 <- c(0, max(finemap_DT$Probability)*1.1)
#   
#   after_plot <- ggplot(finemap_DT, aes(x=POS, y=Probability, label=SNP, color= -log10(P) )) +
#     # ylim(yLimits) +
#     geom_hline(yintercept=0,alpha=.5, linetype=1, size=.5) +
#     geom_point(alpha=.5) +
#     geom_segment(aes(xend=POS, yend=yLimits2[1], color= -log10(P)), alpha=.5) +
#     geom_point(data=labelSNPs[labelSNPs$type=="before",], pch=21, fill=NA, size=4, colour=labelSNPs[labelSNPs$type=="before","color"], stroke=1) +
#     geom_point(data=labelSNPs[labelSNPs$type=="after",][1,], pch=21, fill=NA, size=4, colour=labelSNPs[labelSNPs$type=="after","color"][1], stroke=1) +
#     geom_text_repel(data=labelSNPs, aes(label=SNP), color=labelSNPs$color, segment.alpha = .5, nudge_x = .5, box.padding = .5) +
#     labs(title=paste("After Fine-mapping: ",finemap_method,sep=""),
#          subtitle=paste(gene," (",length(finemap_DT$Probability)," variants)", sep=""),
#          y=score_dict[finemap_method][[1]], x="Position",
#          color="-log10(p-value)") +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
#     scale_x_continuous(breaks = roundBreaks)
#   
#   # Combine plots 
#   # Add COJO plot
#   if (plot_COJO){
#     cojo_plot <- COJO_plot(results_path, gene)
#     plot_grid(before_plot,cojo_plot, after_plot, ncol = 1) %>% print() 
#     # print(cojo_plot)
#   } else {
#     plot_grid(before_plot, after_plot, ncol = 1) %>% print() 
#   }
#   createDT_html(finemap_DT, paste("SUSIE Results: ", gene), scrollY = 200)
# } 
