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

#!/usr/bin/env Rscript

##### Read arguments ######
args <- as.character(commandArgs(trailingOnly=TRUE))
id_input <- match("-i",args)+1
if (is.na(id_input)) {print("Error: No input file specified. Use option -i");quit(save="no",status=1)}
INPUT <- args[id_input]
if (!file.exists(INPUT)){
    print(paste("Error: Input file ", INPUT, " does not exist!", sep=""))
    quit(save="no",status=1)
}
id_output <- match("-o",args)+1
if (is.na(id_output)) {print("Error: No output file specified. Use option -o");quit(save="no",status=1)}
OUTPUT <- args[id_output]
id_link <- match("-l",args)+1
if (is.na(id_link)) {print("Error: No link file specified. Use option -l");quit(save="no",status=1)}
LINK <- args[id_link]
if (!file.exists(LINK)){
    print(paste("Error: Link file ", LINK, " does not exist!", sep=""))
    quit(save="no",status=1)
}
id_pthresh <- match("-pt",args)+1
if (is.na(id_pthresh)) {print("Error: No p-value threshold was specified. Use option -pt, e.g. -pt 1e-5,1e-6,1e-7,1e-8 ");quit(save="no",status=1)}
PTHRESH <- args[id_pthresh]
PTS <- as.numeric(unlist(strsplit(PTHRESH,",",fixed=TRUE)))
id_binning <- match("-b",args)+1
if (is.na(id_binning)) {print("Error: No binning was specified. Use option -b, e.g. -b m5,n5,t5 ");quit(save="no",status=1)}
BINNING <- args[id_binning]
bins <- cbind(substring(unlist(strsplit(BINNING, ",",fixed=T)),1,1),substring(unlist(strsplit(BINNING, ",",fixed=T)),2))
M <- as.numeric(bins[match("m",bins[,1]),2])
N <- as.numeric(bins[match("n",bins[,1]),2])
T <- as.numeric(bins[match("t",bins[,1]),2])
if (is.na(T)){T <-5; print("Error: No number of bins for the TSS distance was specified. Please check option -b, e.g. -b m5,n5,t5 ");quit(save="no",status=1) }
if (is.na(M)){M <-5; print("Error: No number of bins for the MAF distribution was specified. Please check option -b, e.g. -b m5,n5,t5 ");quit(save="no",status=1) }
if (is.na(N)){N <-5; print("Error: No number of bins for the number of LD proxies was specified. Please check option -b, e.g. -b m5,n5,t5 ");quit(save="no",status=1) }
id_subset <- match("-s",args)+1
SUBSET <- args[id_subset]
if (is.na(id_subset)) {print("Note: Testing all annotations, if you are only interested in a subset of them you could use option -s, e.g. -s 1-10 "); SUBSET <- 0}
id_condition <- match("-c",args)+1
CONDITION <- as.numeric(args[id_condition])
if (is.na(id_condition)) {print("Note: No condition option specified. Using no conditioning by default (-c 0). If you are interested in running conditional analysis use option -c 1"); CONDITION <- 0}
if (CONDITION!=1 & CONDITION!=0) {print("Error: Condition option expects a value of 0 (single annotation) or 1 (conditioning)."); quit(save="no",status=1)}
id_condition_thresh <- match("-ct",args)+1
CONDITION_THRESH <- as.numeric(args[id_condition_thresh])
id_padjusted <- match("-padj",args)+1
PADJUSTED <- as.numeric(args[id_padjusted])
if (CONDITION==1){
    if (is.na(id_condition_thresh)) {print("Note: No condition threshold option specified. Using -ct 0.05 by default"); CONDITION_THRESH <- 0.05}
    if (CONDITION_THRESH>1 | CONDITION_THRESH<0) {print("Error: Condition threshold should be between 0 and 1."); quit(save="no",status=1)}
    if (is.na(id_padjusted)) {print("Error: No P-value theshold for significant enrichment option specified. Please use -padj option"); quit(save="no",status=1)}
    if (PADJUSTED>1 | PADJUSTED<0) {print("Error: P-value theshold for significant enrichment should be between 0 and 1."); quit(save="no",status=1)}
    if (!file.exists(OUTPUT)){
        print(paste("Error: File ", OUTPUT, " does not exist. To run conditional analysis please run GARFIELD analysis of one annotation at a time (-c 0)!", sep=""))
        quit(save="no",status=1)
    }
}

##### Print arguments #####
print("########################################")
print("########  GARFIELD PARAMETERS  #########")
print(paste("Input file: ",INPUT, sep=""))
print(paste("Output file: ",OUTPUT, sep=""))
print(paste("Link file: ",LINK, sep=""))
print(paste("Pvalue thresholds: ",PTHRESH, sep=""))
print(paste("Binning parameters: ",BINNING, ". T=",T,", M=",M,", N=",N, sep=""))
print(paste("Annotations to use (default 0 denoting all): ",SUBSET, sep=""))
print(paste("Type of analysis - single annotation (0), conditional (1) :",CONDITION, sep=""))
if (CONDITION==1){
    print(paste("Threshold for significant enrichment: ",PADJUSTED, sep=""))
    print(paste("Threshold for conditional analysis: ",CONDITION_THRESH, sep=""))
    print(paste("Performing analysis at min Pvalue threshold specified: ",min(PTS), sep=""))
}
print("########################################")

##### Read input file ######
print(paste("Reading input file ", INPUT, sep=""))
a <- read.table(INPUT, colClasses="character",header=FALSE)
ncol <- ncol(a)

##### Define bins ######
print("Defining MAF, TSS distance and number of LD proxy bins")
ln <- length(unlist(strsplit(as.character(a[1,7]),"",fixed=TRUE)))
mat <- matrix(as.numeric(unlist(strsplit(as.character(a[,7]),"",fixed=TRUE))),nc=ln,byrow=TRUE)
a[,2] <- as.numeric(a[,2]) # position
a[,3] <- as.numeric(a[,3]) # pval
a[,4] <- as.numeric(a[,4]) # maf
a[,5] <- as.numeric(a[,5]) # tssd
a[,6] <- as.numeric(a[,6]) # ntags r2>0.8

if (T>1) {
    Tq <- cut(a[,5],as.numeric(quantile(a[,5],seq(0,1,length.out=T+1))), include.lowest = TRUE)
} else {
    Tq = rep(1,nrow(a))
}

if (M>1){
    tmp1 <- quantile(a[,4],seq(0,1,length.out=M+1))
    qm <- NULL
    if (tmp1[2]!=tmp1[1]){
        qm <- c(-1,tmp1[2])
        for (i in M:2){
            tmp2 <- quantile(a[which(a[,4] > tmp1[2]),4],seq(0,1,length.out=i))
            qm <- c(qm,tmp2[2])
            tmp1 <- tmp2
        }
    } else {
        qm <- c(-1,tmp1[2])
        for (i in (M):2){
            tmp2 <- quantile(a[which(a[,4] > tmp1[2]),4],seq(0,1,length.out=i))
            qm <- c(qm,tmp2[2])
            tmp1 <- tmp2
        }
    }
    Mq <- cut(a[,4],as.numeric(qm), include.lowest = TRUE)
} else {
    Mq = rep(1,nrow(a))
}

if (N>1) {
    tmp1 <- quantile(a[,6],seq(0,1,length.out=N+1))
    qn <- NULL
    if (tmp1[2]!=tmp1[1]){
        qn <- c(-1,tmp1[2])
        for (i in N:2){
	       tmp2 <- quantile(a[which(a[,6] > tmp1[2]),6],seq(0,1,length.out=i))
	       qn <- c(qn,tmp2[2])
	       tmp1 <- tmp2
        }
    } else {
        qn <- c(-1,tmp1[2])
        for (i in (N):2){
	       tmp2 <- quantile(a[which(a[,6] > tmp1[2]),6],seq(0,1,length.out=i))
	       qn <- c(qn,tmp2[2])
	       tmp1 <- tmp2
        }
    }
    Nq <- cut(a[,6],as.numeric(qn), include.lowest = TRUE)
} else {
    Nq = rep(1,nrow(a))
}


features <- ""
if (M>1) {features <- paste(features,"Mq+", sep="")} 
if (T>1) {features <- paste(features,"Tq+", sep="")} 
if (N>1) {features <- paste(features,"Nq+", sep="")} 
print(paste("Accounting for features : ",features,sep=""))

#### Define response variables ######
print("Defining response variables according to p-value thresholds")
Pmat <- NULL
nt <- NULL
for (t in PTS){
	Pmat <- cbind(Pmat,as.numeric(a[,3]<=t))
	nt <- c(nt, sum(as.numeric(a[,3]<=t)))
}
n <- nrow(Pmat)

##### Process the subset of annotations to run enrichment for ######
get_indices <- function(x){seq(as.numeric(x[1]), as.numeric(x[2]), 1)}
if (SUBSET==0){
    SUBSET <- paste(1, "-", ncol(mat), sep="")
}
subsetIndices <- unlist(lapply(strsplit(unlist(strsplit(SUBSET, split=",", fixed=TRUE)), split="-", fixed=TRUE), get_indices))

##### Read link file ######
print(paste("Reading link file ", LINK, sep=""))
lnk <- read.table(LINK, header =TRUE)

if (CONDITION==0){
    print("Performing enrichment analysis one annotation at a time")
    write.table("ID PThresh OR Pvalue Beta SE CI95_lower CI95_upper NAnnotThesh NAnnot NThresh N linkID Annotation Celltype Tissue Type Category",file=OUTPUT, row.names=FALSE, col.names=FALSE, append=FALSE, quote=FALSE)
    for (i in subsetIndices){
	   na <- sum(mat[,i])
	   dat <- NULL
	   for (j in 1:length(PTS)){
		  t <- PTS[j]
		  formula <- paste("Pmat[,j]~",features,"mat[,i]",sep="")
		  model <- glm(formula, family='binomial')
		  smr <- summary(model)$coefficients
		  smrC <- smr[nrow(smr),]
		  nat <- sum(mat[,i]*as.numeric(a[,3]<=t))
		  dat <- rbind(dat,c(i,t,exp(smrC[1]), smrC[4],smrC[1],smrC[2],smrC[1]-1.96*smrC[2],smrC[1]+1.96*smrC[2], nat, na, nt[j], n))
	   }
	   print(i)
       dat = as.data.frame(dat)
       dat$linkID = lnk[i,1]
       dat$Annotation = lnk[i,2]
       dat$Celltype = lnk[i,3]
       dat$Tissue = lnk[i,4]
       dat$Type = lnk[i,5]
       dat$Category = lnk[i,6]
       write.table(dat,file=OUTPUT, row.names=FALSE, col.names=FALSE, append=TRUE, quote=FALSE)
    }
    print("Analysis complete")
} else {
    print("Performing conditional analysis - model selection")
    dat <- read.table(OUTPUT, header=TRUE, comment.char="")
    dat <- dat[which(dat[,1] %in% (subsetIndices)),] ## Filter based on SUBSET
    dat <- dat[which(dat[,2] %in% min(PTS)),] ## Filter based on min PTHRESH
    indices <- order(dat[,4], decreasing = FALSE)[which(dat[,4][order(dat[,4], decreasing = FALSE)]<PADJUSTED)] # order p-values in increasing order
    if (length(indices)==0) {print("No annotations reached significance. Nothing to be done.");quit(save="no",status=0)}
    mat <- mat[,subsetIndices]
    m0 <- mat[,indices]
    t <- min(PTS)
    dataS <- data.frame(P=Pmat[,1],MAF=Mq,TSS=Tq,NTAGS=Nq, m0)
    dataS0 <- data.frame(P=Pmat[,1],MAF=Mq,TSS=Tq,NTAGS=Nq, m0[,1])
    cols <- 1
    if (M>1) cols <- c(cols,2)
    if (T>1) cols <- c(cols,3)
    if (N>1) cols <- c(cols,4)
    cols <- c(cols,5)
    model0 <- glm(P ~ . , family='binomial', data = dataS0)
    smr <- summary(model0)$coefficients
    smrC <- smr[nrow(smr),]
    dat2 <- NULL
    write.table("ID_c ID PThresh OR_c Pvalue_c Beta_c SE_c CI95_lower_c CI95_upper_c OR Pvalue Beta SE CI95_lower CI95_upper NAnnotThesh NAnnot NThresh N linkID Annotation Celltype Tissue Type Category",file=paste(OUTPUT,".",CONDITION_THRESH,".cond.indep",sep=""), row.names=FALSE, col.names=FALSE, append=FALSE, quote=FALSE)
    write.table("ID_c ID PThresh OR_c Pvalue_c Beta_c SE_c CI95_lower_c CI95_upper_c OR Pvalue Beta SE CI95_lower CI95_upper NAnnotThesh NAnnot NThresh N linkID Annotation Celltype Tissue Type Category",file=paste(OUTPUT,".",CONDITION_THRESH,".cond.dep",sep=""), row.names=FALSE, col.names=FALSE, append=FALSE, quote=FALSE)
    dat2 <- rbind(dat2,c(1,dat[indices,1][1],t,exp(smrC[1]), smrC[4],smrC[1],smrC[2],smrC[1]-1.96*smrC[2],smrC[1]+1.96*smrC[2]))
    dat2 <- as.data.frame(dat2)
    dat2 <- cbind(dat2,dat[indices[1],3:18])
    write.table(dat2,file=paste(OUTPUT,".",CONDITION_THRESH,".cond.indep",sep=""), row.names=FALSE, col.names=FALSE, append=TRUE, quote=FALSE)
    if (length(indices)==1) {print("Only one annotations reached significance. Nothing more to be done.");quit(save="no",status=0)}
    modeltemp <- model0
    for (i in 2:(length(indices))){
        cols1 <- c(cols,i+4)
        dataS1 <- dataS[,cols1]
        dat2 <- NULL
        dat3 <- NULL
        model1 <- glm(P ~ . , family='binomial', data = dataS1)
        smr <- summary(model1)$coefficients
        smrC <- smr[nrow(smr),]
        pvan <- anova(model1,modeltemp,test="Chisq")[[5]][2]
        if (pvan > CONDITION_THRESH){
            cols1 <- cols
            dat3 <- rbind(dat3,c(i,dat[indices,1][i],t,exp(smrC[1]), smrC[4],smrC[1],smrC[2],smrC[1]-1.96*smrC[2],smrC[1]+1.96*smrC[2]))
            dat3 <- as.data.frame(dat3)
            dat3 <- cbind(dat3,dat[indices[i],3:18])  
            write.table(dat3,file=paste(OUTPUT,".",CONDITION_THRESH,".cond.dep",sep=""), row.names=FALSE, col.names=FALSE, append=TRUE, quote=FALSE)
        } else {
        	modeltemp <- model1
            cols <- cols1
            dat2 <- rbind(dat2,c(i,dat[indices,1][i],t,exp(smrC[1]), smrC[4],smrC[1],smrC[2],smrC[1]-1.96*smrC[2],smrC[1]+1.96*smrC[2]))
            dat2 <- as.data.frame(dat2)
            dat2 <- cbind(dat2,dat[indices[i],3:18]) 
            write.table(dat2,file=paste(OUTPUT,".",CONDITION_THRESH,".cond.indep",sep=""), row.names=FALSE, col.names=FALSE, append=TRUE, quote=FALSE)
        }
        print(i)
    }
    dataS1 <- dataS[,cols]
    nms = dat[indices,][cols[-c(grep('X',names(dataS1),invert=TRUE))]-4,]
    nms = paste(nms[,1],nms[,15],nms[,16],nms[,17],nms[,18],sep=",")
    model1 <- glm(P ~ . , family='binomial', data = dataS1)
    finalsummary <- summary(model1)$coefficients
    finalsummary <- as.data.frame(finalsummary)
    finalsummary$Terms = rownames(finalsummary)
    finalsummary$Terms[grep('X', finalsummary$Terms)] = nms
    write.table(finalsummary[,c(5,1:4)], file=paste(OUTPUT,".",CONDITION_THRESH,".model.summary",sep=""), row.names=FALSE, col.names=TRUE, append=FALSE, quote=FALSE, sep="\t")
    print("Analysis complete")
}
