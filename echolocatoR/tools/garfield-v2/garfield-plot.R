# GARFIELD - GWAS analysis of regulatory or functional information enrichment with LD correction.
# Copyright (C) 2016 Wellcome Trust Sanger Institute / EMBL - European Bioinformatics Institute
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#### Rscript software_plot.R in.file 100000 out.file.prefix link.file Phenotype 10

#### Read input parameters ####
	#options(echo=TRUE) # if you want see commands in output file
	args <- commandArgs(trailingOnly = TRUE)
	id_input <- match("-i",args)+1
	if (is.na(id_input)) {print("Error: No input file specified. Use option -i");quit(save="no",status=1)}
	input_file <- args[id_input]
	if (!file.exists(input_file)){
    	print(paste("Error: Input file ", input_file, " does not exist!", sep=""))
    	quit(save="no",status=1)
	}
	id_output <- match("-o",args)+1
	if (is.na(id_output)) {print("Error: No output file specified. Use option -o");quit(save="no",status=1)}
	output_prefix <- args[id_output]
	#id_link <- match("-l",args)+1
	#if (is.na(id_link)) {print("Error: No link file specified. Use option -l");quit(save="no",status=1)}
	#link_file <- args[id_link]
	#if (!file.exists(link_file)){
	#    print(paste("Error: Link file ", link_file, " does not exist!", sep=""))
	#    quit(save="no",status=1)
	#}
	id_filter <- match("-f",args)+1
	if (is.na(id_filter)) {print("No threshold for the mininum number of variants at P-value threshold was specified.");filter <-0 }
	filter <- as.numeric(args[id_filter])
	id_title <- match("-t",args)+1
	plot_title <- args[id_title]
	if (is.na(id_title)) {plot_title=""}

	id_padjusted <- match("-padj",args)+1
	padjusted <- as.numeric(args[id_padjusted])
    if (is.na(id_padjusted)) {print("Error: No P-value theshold for significant enrichment option specified. Please use -padj option"); quit(save="no",status=1)}
    if (padjusted>1 | padjusted<0) {print("Error: P-value theshold for significant enrichment should be between 0 and 1."); quit(save="no",status=1)}
	id_subset <- match("-s",args)+1
	SUBSET <- args[id_subset]
	if (is.na(id_subset)) {SUBSET <- 0}

	fixed_colours <- 0
	lthresh <- -1
	id_colour <- match("-col",args)+1
	if (!is.na(id_colour)) {
		col.thresh <- unlist(strsplit(args[id_colour], ",", fixed=T))
		fixed_colours <- 1
		lthresh <- length(col.thresh)
	}
	
##### Print arguments #####
print("########################################")
print("########  GARFIELD PLOT PARAMETERS  #########")
print(paste("Input file: ",input_file, sep=""))
print(paste("Output prefix file: ",output_prefix, sep=""))
print(paste("Filter thresholds for minimun number of variants at Pvalue threshold: ",filter, sep=""))
print(paste("Plot title: ",plot_title, sep=""))
print(paste("Threshold for significant enrichment: ",padjusted, sep=""))
print(paste("Annotations to plot (default 0 denoting all): ",SUBSET, sep=""))
if (!is.na(id_colour)) {
print(paste("Using custom colours if of the right length, otherwise using default sets. Custom colours supplied: ", paste(col.thresh,collapse=","), sep=""))
}
print("########################################")

#### Read files
	source("garfield-plot-function.R") # load R function to create plot
	input <- read.table(input_file, header=TRUE, comment.char="") 
	get_indices <- function(x){seq(as.numeric(x[1]), as.numeric(x[2]), 1)}
	if (SUBSET==0){
	    SUBSET <- paste(1, "-", length(unique(input$ID)), sep="")
	} else {
		subsetIndices <- unlist(lapply(strsplit(unlist(strsplit(SUBSET, split=",", fixed=TRUE)), split="-", fixed=TRUE), get_indices))
		input <- input[which(input$ID %in% subsetIndices),]
	}
	input <- unique(input[which(as.numeric(as.character(input$NThresh))>=filter),])
	#link <- read.table(link_file, header=TRUE)

#### Create plots 	
	input$Category <- as.character(input$Category)
	if (length(which(is.na(input$Category)))>1){ 
		input$Category[which(is.na(input$Category))] <- "Custom" 
	}
	input$Category <- as.factor(input$Category)

	for (category in as.character(unique(input$Category))) {
		ids <- which(input$Category==category)
		thresholds <- sort(unique(input$PThresh[ids]))
		annotations <- unique(as.character(input$ID[ids]))
		if (category %in% c("Genic","Histone_Modifications","Chromatin_States")){
			tissues <- as.character(input$Type[ids][match(annotations,input$ID[ids])])
			nms <- as.character(input$Celltype[ids][match(annotations,input$ID[ids])])
			compact <- FALSE
			tissue_label <- "Feature"
			if (category %in% c("Genic")){
				nms <- tissues
			}
		} else if (category %in% c("TFBS","FAIRE","Hotspots","Peaks","Footprints")){
			nms <- as.character(input$Celltype[ids][match(annotations,input$ID[ids])])
			tissues <- as.character(input$Tissue[ids][match(annotations,input$ID[ids])])
			compact <- TRUE
			tissue_label <- "Tissue"
			if (category %in% c("Hotspots","Peaks","Footprints")){
				nms <- tissues
			}		
		} else {
			nms <- as.character(input$Celltype[ids][match(annotations,input$ID[ids])])
			tissues <- as.character(input$Tissue[ids][match(annotations,input$ID[ids])])
			compact <- FALSE
			tissue_label <- "Annotation"
			nms <- tissues = as.character(input$Annotation[ids][match(annotations,input$ID[ids])])
			nc0 <- strsplit(nms,"/",fixed=TRUE)
			nms <- tissues <- matrix(unlist(nc0),nr=length(nc0[[1]]))[length(nc0[[1]]),]
		}

		DATA <- DATA_p <- matrix(NA, nrow=length(thresholds)+1, ncol=length(annotations))
		for (j in 1:length(annotations)){
			for (i in 1:length(thresholds)){
				ii <- which(input$ID==annotations[j] & input$PThresh==thresholds[i])
				DATA[i,j] <- input$OR[ii]
				DATA_p[i,j] <- input$Pvalue[ii]
				if (input$Pvalue[ii]==-1){ DATA_p[i,j] <- 2 }
			}
		}
		DATA[length(thresholds)+1,] <- DATA_p[length(thresholds)+1,] <- 1
		DATA_p <- -log10(DATA_p)

		ann.col <- colorRampPalette(c("tomato","skyblue3","yellow","brown2","lightgreen","lightgoldenrod3","purple","pink","darkblue","gray","darkgreen"))( length(unique(tissues)) )[as.numeric(as.factor(tissues))]

		tr <- -log10(padjusted)
		ann.col.new <- matrix(ann.col,nrow=length(thresholds),ncol=length(ann.col),byrow=TRUE)
		for (i in 1:length(thresholds)){		
			ann.col.new[i,which(DATA_p[i,]<tr)] <- 0
		}

		ord <- order(tissues)

		if (fixed_colours==0 | lthresh!=(length(unique(thresholds))+1)){
			col.thresh <- colorRampPalette(c("black","firebrick3","tomato","RosyBrown2","dodgerblue3","skyblue2","lightskyblue1","gray70","blanchedalmond"))( length(unique(thresholds))+1 )
		}
		if (length(thresholds)<4){
			rws <- length(thresholds):1
		} else {
			rws <- 4:1
		}
		tissues <- gsub("_", " ", as.character(tissues))
		nms <- gsub("_", " ", as.character(nms))

		pdf(paste(output_prefix,".",category,".pdf",sep=""),12,10)
			layout(matrix(c(1,2,3,3),nr=2), widths = c(9,2),heights = c(9, 1), respect = FALSE)
			par(oma = c(0, 0, 0, 0))
			garfield.plot(DATA[,ord],ann.cols=matrix(ann.col.new[rws,ord],nrow=length(rws),ncol=length(ord)), ann.col=ann.col[ord], ann.pch=15, rp.type="p",line.col=col.thresh,show.grid=TRUE, show.radial.grid=TRUE,labels=nms[ord],breaks=tissues[ord], radlab=TRUE,cex.axis=0.1, cex.lab=0.1, mar = c(6, 2, 5, 1), label.prop=1.1,poly.col=col.thresh, compact=compact)
			par(mar = c(0, 0, 0, 0))
			title(main=paste(plot_title," ",category,sep=""),line=7,cex=2) 
			plot(1:2,1:2, type="n", axes=FALSE, xlab="", ylab="")
			legend("bottom", c(thresholds,"1"),col=col.thresh,lty=1,lwd=6,title="GWAS P-value Threshold",horiz=TRUE,cex=1,bty="n")
			plot(1:2,1:2, type="n", axes=FALSE, xlab="", ylab="")
			legend("right",unique(tissues[ord]),col=unique(ann.col[ord]),lty=1,lwd=5,title=tissue_label,cex=1,bty="n")
		dev.off()

	}

print("Figures and table of results created!")
