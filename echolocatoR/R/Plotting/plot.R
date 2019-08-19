# %%%%%%%%%%%%%%%%% #
####### PLOT ####### 
# %%%%%%%%%%%%%%%%% # 

 
# cojo_DT = data.table::fread(file.path("Data/GWAS/Nalls23andMe_2019/LRRK2/COJO/COJO_results.txt"), sep="\t")

# SNP_list = c("rs76904798","rs34637584","rs117073808")
COJO_plot <- function(gene, 
                      cojo_DT, 
                      results_path, 
                      conditioned_snps, 
                      show_plot=T, 
                      save_plot=T){
  gene <- basename(results_path)
  # Label independent SNPs
  ind_SNPs <- subset(cojo_DT, COJO.Credible_Set==T)
  ind_SNPs <- cbind(ind_SNPs, type="Independent",color="purple") 
  labelSNPs <- ind_SNPs#rbind(labelSNPs, ind_SNPs) 
  
  effect_SNPs <- cojo_DT %>% arrange(desc(abs(COJO.Conditioned_Effect)))
  effect_SNPs <- effect_SNPs[1:5,]
  effect_SNPs$type <- "effect_SNPs"
  effect_SNPs$color <- "blue3"
  labelSNPs <- rbind(labelSNPs, effect_SNPs, fill=T)
  
  # topEffect_snps <- cojo_DT %>% arrange(desc(COJO.Conditioned_Effect))
  # topEffect_snps <- topEffect_snps$SNP[1:5]
  # cojo_DT <- cojo_DT %>% dplyr::mutate(color = ifelse(COJO.Credible_Set | COJO.Conditioned_Effect %in% topEffect_snps,
  #                                                     ifelse(COJO.Credible_Set & COJO.Conditioned_Effect %in% topEffect_snps, "purple", 
  #                                                            ifelse(COJO.Credible_Set & !(COJO.Conditioned_Effect %in% topEffect_snps), "lightpurple", "lightblue") ), NA)
  # )
 
  spacing <- if(length(cojo_DT$SNP)>1000){250000}else{50000}
  # roundBreaks <- seq(plyr::round_any(min(cojo_DT$POS),10000), max(cojo_DT$POS),spacing) 
  
  cp <- ggplot(cojo_DT, aes(x=POS, y = -log10(P), label=SNP, color= -log10(P))) +
    # ylim(yLimits1) +
    geom_hline(yintercept=0, alpha=.5, linetype=1, size=.5) +
    geom_point(alpha=.5) +
    geom_segment(aes(xend=POS, yend=0, color= -log10(P) ), alpha=.5) +
    geom_point(data=labelSNPs, pch=21, size=4, colour=labelSNPs$color, stroke=1) + 
    geom_label_repel(data=labelSNPs, aes(label=SNP), color=NA, nudge_x = .5, box.padding = .5,
                     label.size=NA, alpha=.8, seed = 1) +
    geom_label_repel(data=labelSNPs, aes(label=SNP), color=labelSNPs$color, segment.alpha = .5, 
                     nudge_x = .5, box.padding = .5, fill = NA, alpha=1, seed = 1, show.legend = T ) +
    labs(title=paste(gene,": Conditional & Stepwise Results (COJO)"),
         subtitle = paste("Purple = Independent signals from stepwise procedure\n",
                          "Blue = Residual effects conditioned on:",conditioned_snps),
         y="-log10(p-value)", x="Position", color="-log10(p-value)") +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5) )
    # scale_x_continuous(breaks = roundBreaks)
  
  if(save_plot){
    png(filename=file.path(results_path,"COJO/COJO_plot.png"), width = 1000, 800)
    print(cp)
    dev.off()
  }  
  if(show_plot){print(cp)}else{return(cp)} 
}





construct_SNPs_labels <- function(DT, lead=T, method=T, consensus=T){ 
  printer("\n + Constructing SNP labels...")
  labelSNPs <- data.table::data.table() 
  DT <- data.table::as.data.table(DT)
  
  ## BEFORE fine-mapping  
  if(lead){
    before <- subset( DT %>% arrange(P), leadSNP == T)[1,]  
    before$type <- "Lead SNP"
    before$color <- "red"
    labelSNPs <- rbind(labelSNPs, before, fill=T)
  }
  if(method){
    # AFTER fine-mapping
    after = subset(DT, Support>0) 
    if(dim(after)[1]>0){
      after$type <- "Credible Set"  
      after$color<- "green3"
      labelSNPs <- rbind(labelSNPs, after, fill=T) 
    } 
  } 
  if(consensus & "Consensus_SNP" %in% colnames(DT)){
    # Conensus across all fine-mapping tools
    cons_SNPs <- subset(DT, Consensus_SNP==T)
    if(dim(cons_SNPs)[1]>0){
      cons_SNPs$type <- "Consensus SNP"
      cons_SNPs$color <- "darkgoldenrod1"
      labelSNPs <- rbind(labelSNPs, cons_SNPs, fill=T)
    } 
  } 
  # Convert to GRanges object
  # labelGR <- transformDfToGr(data=labelSNPs, seqnames = "CHR", start = "POS", end = "POS")
  # names(labelGR) <- labelGR$SNP
  # plotGrandLinear(gr.snp, aes(y = P, x=POS), highlight.gr = labelGR)
  return(as.data.frame(labelSNPs))
}
 
LD_with_leadSNP <- function(LD_matrix, LD_SNP){  
  # Get R2 values from LD matrix
  printer("LD Matrix dimensions:", paste(dim(LD_matrix), collapse=" x "))
  printer("Extracting LD subset for lead SNP:",LD_SNP)
  LD_sub <- subset(LD_matrix, select=LD_SNP) %>% 
    as.data.table(keep.rownames = T) %>% 
    `colnames<-`(c("SNP","r")) %>% 
    dplyr::mutate(r2 = r^2)
  return(LD_sub)
}






snp_plot <- function(finemap_DT,
                     LD_matrix,
                     gene,
                     method="original",
                     show_plot=T,
                     subtitle=NA,
                     multi = T,
                     LD_SNP = NA){
  { 
  # finemap_DT <- data.table::fread("Data/GWAS/Nalls23andMe_2019/LRRK2/Multi-finemap/Multi-finemap_results.txt", sep="\t")
  score_dict <- setNames(c("-log10(P-value)", "PIP","PP", "PIP", "Conditional Probability","-log10(P-value)"),
                         c("original","SUSIE", "ABF", "FINEMAP", "COJO","COLOC"))
  # X-tick spacing
  spacing <- if(length(finemap_DT$SNP)>1000){250000}else{50000}
  roundBreaks <- seq(plyr::round_any(min(finemap_DT$POS),10000), max(finemap_DT$POS), spacing) 
  # yLimits2 <- c(0,1.1)
  
  # Merge with LD info
  # If not specified, identify a lead SNP by summing the PPs from each fine-mapping method
  # if(is.na(LD_SNP)){
    LD_SNP <- subset(finemap_DT, leadSNP==T)$SNP
    LD_sub <- LD_with_leadSNP(LD_matrix, LD_SNP)
    DT <- data.table:::merge.data.table(finemap_DT, LD_sub, by = "SNP")
  # }  
  
  if(method=="original"){   
    is.na(DT$Probability) <- 0
    p <- ggplot(data = DT, aes(x=POS, y= -log10(P), label=SNP, color= r2))
# =======
#   if(method=="original"){
#     DT <- finemap_DT 
#     
#     is.na(DT$Probability) <- 0 
#     
#     p <- ggplot(data = DT, aes(x=POS, y= -log10(P), label=SNP, color= -log10(P) ))
# >>>>>>> 1e2aecb9b38f6c049a9c6f1d9baed0f0d268e0b4:echolocatoR/R/plot.R
    title <- paste0(gene," : Before fine-mapping")
    tag_SNPs <- labelSNPs <- construct_SNPs_labels(DT, lead=T, method = F, consensus = T)
    subtitle <- if(is.na(subtitle)){paste0(length(DT$SNP)," SNPs")}else{subtitle}
  } else {
    
    # COLOC
    if(method=="COLOC"){
      DT <- DT %>% dplyr::rename(Credible_Set = "Colocalized" )
      title <- paste0(gene," : Colocalization (",method,")") 
      p <- ggplot(data = DT, aes(x=POS, y= -log10(P), label=SNP, color= r2 ))  
    # Fine-mapping methods
    } else if (multi){
      title <- paste0(gene," : After fine-mapping (",method,")") 
      DT <- DT %>% dplyr::rename(Probability = paste0(method,".Probability"), 
                                 Credible_Set = paste0(method,".Credible_Set")) 
      is.na(DT$Probability) <- 0 
      p <- ggplot(data = DT, aes(x=POS, y=Probability, color= r2 )) + 
# =======
#       DT <- finemap_DT %>% dplyr::rename(Probability = paste0(method,".Probability"), 
#                                          Credible_Set = paste0(method,".Credible_Set")) 
#       is.na(DT$Probability) <- 0 
#       p <- ggplot(data = DT, aes(x=POS, y=Probability, label=SNP, color= -log10(P) )) + 
# >>>>>>> 1e2aecb9b38f6c049a9c6f1d9baed0f0d268e0b4:echolocatoR/R/plot.R
        ylim(c(0,1.1)) 
      subtitle <- if(is.na(subtitle)){
        paste0(length(subset(DT, Credible_Set>0)$SNP), " Candidate SNP(s)")
        }else{subtitle}
    } else {
      title <- paste0(gene," : After fine-mapping (",method,")")  
      p <- ggplot(data = DT, aes(x=POS, y=Probability,  color= r2 )) + 
        ylim(c(0,1.1)) 
    } 
    labelSNPs <- construct_SNPs_labels(DT, lead=T, method = T, consensus = F)
    tag_SNPs <- subset(labelSNPs, Credible_Set>0)
  }
  
    
  p <- p + geom_hline(yintercept=0,alpha=.5, linetype=1, size=.5) +
    stat_smooth(data=DT, aes(x=POS, y=r2, fill=r2),color="darkgray", 
                se = F, formula = y ~ x, 
                method = 'loess', span=.05) + 
    geom_point() + # alpha=.5 
    # geom_segment(aes(xend=POS, yend=0, color= -log10(P)), alpha=.5) +
    geom_point(data=labelSNPs, pch=21, fill=NA, size=4, color=labelSNPs$color, stroke=1) +
    geom_point(data=subset(DT, SNP==LD_SNP), pch=18, fill=NA, size=4, color="red") + 
    # Labels (one for background, one for text)
    geom_label_repel(data=tag_SNPs, aes(label=SNP), 
                     color=NA, 
                     nudge_x = .5, 
                     box.padding = .5,
                     label.size=NA, 
                     alpha=.5, 
                     seed = 1,) +
    geom_label_repel(data=tag_SNPs, aes(label=SNP), 
                     color=tag_SNPs$color,
                     segment.alpha = .5, 
                     nudge_x = .5, 
                     box.padding = .5, 
                     fill = NA, 
                     alpha=1, 
                     seed = 1) +
    labs(title = title,
         subtitle = subtitle,
         y = score_dict[[method]], 
         x = "Position",
         color = bquote(paste(R^2," with ",.(LD_SNP) ) )) +
    theme(plot.title = element_text(hjust = 0.5), 
          plot.subtitle = element_text(hjust = 0.5)) +
    scale_x_continuous(breaks = roundBreaks) +
    scale_color_gradient(low="blue", high="red", limits = c(0,1))
  if(show_plot){print(p)} else{ return(p) }
  }

}
 


multi_finemap_plot <- function(finemap_DT,
                               LD_matrix,
                               results_path,
                               finemap_method_list, 
                               gene, 
                               conditioned_snps,
                               original=T, 
                               save_plot=T,
                               ncols=1,
                               width=500, #500,
                               height=1000 #1000
                               ){ 
  method_list <- if(original){c("original", finemap_method_list)}else{finemap_method_list} 
  
  # Assemble plots in list 
   plot_list <- lapply(method_list, function(method, 
                                            finemap_DT.=finemap_DT,
                                            LD_matrix.=LD_matrix){   
    printer("\n Plotting...",method)
    if(method=="COJO"){
      p <- COJO_plot(cojo_DT = finemap_DT., 
                results_path = results_path, 
                show_plot = F, 
                save_plot = T,
                conditioned_snps = conditioned_snps)
    } else{
      p <- snp_plot(finemap_DT = finemap_DT., 
                    LD_matrix = LD_matrix.,
                    gene = gene, 
                    method = method, 
                    show_plot = F, 
                    multi = T)
    } 
    return(p) 
  }) 
  # Multi-plot
  ## Adjust plotting parameters
  if(length(finemap_method_list)>3){ncols <- 2; width <- width*2}
  
  # grDevices::graphics.off() 
  cp <- cowplot::plot_grid(plotlist = plot_list, ncol = ncols) 
  # print(cp)
  # plot_list <- lapply(1:3,function(e){ggplot(iris, aes(x=Sepal.Width,y=Sepal.Length))+geom_point()})

  # library(gridExtra)
  # cp <- do.call("grid.arrange", c(plot_list, ncol=ncols))
  # plot(cp)
   
  # Rmisc::multiplot(plotlist = plot_list, cols = ncols)
   
  # Save plot
  if(save_plot){
    png(filename=file.path(results_path,"Multi-finemap/multi_finemap_plot.png"),width = width, height)
    print(cp)
    dev.off()
    # ggplot2::ggsave(file=file.path(results_path,"Multi-finemap/multi_finemap_plot.png"),
    #        cp, width = width, height = height)
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



# remotes::install_github("tylermorganwall/rayshader")
# library(rayshader)
# library(ggplot2)
# gg = ggplot(diamonds, aes(x, depth)) +
#   stat_density_2d(aes(fill = stat(nlevel)), 
#                   geom = "polygon",
#                   n = 100,bins = 10,contour = TRUE) +
#   facet_wrap(clarity~.) +
#   scale_fill_viridis_c(option = "A")
# plot_gg(gg,multicore=TRUE,width=5,height=5,scale=250)
# 



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
