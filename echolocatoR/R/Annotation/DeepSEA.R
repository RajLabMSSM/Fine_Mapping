## ~~~~~~~~~~~~~~~~~~        ~~~~~~~~~~~~~~~~~~ 
# ~~~~~~~~~~~~~~~~~~ DeepSEA ~~~~~~~~~~~~~~~~~~ 
## ~~~~~~~~~~~~~~~~~~        ~~~~~~~~~~~~~~~~~~ 

# Results generated using DeepSEA online server(http://deepsea.princeton.edu/help/)
## or the selene python module (https://github.com/FunctionLab/selene/blob/master/tutorials/analyzing_mutations_with_trained_models/analyzing_mutations_with_trained_models.ipynb).

DeepSEA.prepare_data <- function(results_path){
  
  DeepSEA.merge_data <- function(dat){
    dat.merged <- dat %>% dplyr::mutate(CHR=gsub("chr","",chr), POS=pos) %>% 
      merge.data.frame(finemap_DT,#[,c("CHR","POS","SNP")], 
                       by=c("CHR","POS")) %>% 
      data.table::data.table() %>% arrange(SNP)
    return(dat.merged)
  }
  # Fine-mapping results
  finemap_DT <- data.table::fread(file.path(results_path,"Multi-finemap/Multi-finemap_results.txt"))
  # SNP-wise 'Functional Significant Scores'
  funsig <- data.table::fread(file.path(results_path, "DeepSEA/infile.vcf.out.funsig")) %>% DeepSEA.merge_data() 
  # SNP-wise probabilites of being eQTL, GWAS, or HGMD (The Human Gene Mutation Database) hit
  snpclass <- data.table::fread(file.path(results_path, "DeepSEA/infile.vcf.out.snpclass")) %>% DeepSEA.merge_data()
  deepsea.dat <- cbind(snpclass, `Functional significance score`=funsig$`Functional significance score`)
  return(deepsea.dat)
}
# deepsea.dat <- DeepSEA.prepare_data(results_path = "./Data/GWAS/Nalls23andMe_2019/LRRK2")


DeepSEA.corrplot <- function(deepsea.dat, save_path=F){ 
  PP.cols <- grep(".Probability",colnames(deepsea.dat), value = T)
  cor.dat <- deepsea.dat[,c(PP.cols,"mean.PP","eQTL-probability","GWAS-probability","HGMD-probability","Functional significance score")]
  cor.df <- cor(cor.dat)
  cor.df[is.na(cor(cor.dat))] <- 0 
  if(save_path!=F){png(file = file.path(results_path,"DeepSEA","DeepSEA.corrplot.png") )}
  corrplot::corrplot(cor.df,
                     # title = "Correlations Between Fine-mapping PPs and DeepSEA Predictions",
                     # bg = "black",  
                     order = "hclust",
                     # addrect = 2,
                     type = "upper",
                     tl.col = "black", 
                     method = "ellipse",
                     addCoef.col = T,
                     number.cex = 1,
                     tl.srt = 45,
                     tl.cex = .5
                     # col = colorRampPalette(c("red","white","blue"))(200)
                     # col = RColorBrewer::brewer.pal(5,"RdBu")
                     ) 
  dev.off()
}

DeepSEA.plot_predictions <- function(deepsea.dat){ 
  deepsea.melt <- data.table::melt.data.table(data.table::data.table(deepsea.dat.all), 
                                             id.vars = "SNP",
                                             measure.vars = c("eQTL-probability",
                                                              "GWAS-probability",
                                                              "HGMD-probability",
                                                              "Functional significance score"),
                                             variable.name = "Prediction",
                                             value.name = "Probability")
  deepsea.melt$Prediction <- gsub("-| ","\n",deepsea.melt$Prediction)
  deepsea.melt$Prediction <- factor(deepsea.melt$Prediction, levels = unique(deepsea.melt$Prediction) )
  gp <- ggplot(data=deepsea.melt, aes(x=SNP, y=Probability, fill=Probability)) +
    geom_col() + 
    facet_grid(Prediction~.) + 
    theme_dark() + 
    ylim(c(0,1)) +   
    scale_fill_viridis_c() +
    geom_hline(yintercept = 0) + 
    geom_text(stat = 'identity',aes(label=SNP, hjust=-.2, vjust= -.5, srt=45),size=3 ) + 
    theme(axis.title.x=element_blank(), 
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          rect = element_rect(fill = "transparent"),
          panel.background = element_rect(fill = "transparent"),
          strip.text.y = element_text(angle = 0),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +  
    labs(title="DeepSEA Predictions")
  print(gp)
  if(save_path!=F){
    ggsave(filename = file.path(results_path,"DeepSEA","DeepSEA.predictions.png"),
           plot = gp, dpi = 600, height=10, width=10)
  } 
  return(gp)
}
