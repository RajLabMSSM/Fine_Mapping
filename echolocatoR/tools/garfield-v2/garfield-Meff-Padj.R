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

## script to find the effective number of annotations and the adjusted for multiple testing P-value
args <- commandArgs(trailingOnly = TRUE)
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
id_subset <- match("-s",args)+1
SUBSET <- args[id_subset]
if (is.na(id_subset)) {print("Note: Using all annotations, if you are only interested in a subset of them you could use option -s, e.g. -s 1-10 "); SUBSET <- 0}

a <- read.table(INPUT, header = FALSE, colClasses = "character")
mx <- strsplit(a[,ncol(a)], "", fixed=TRUE)
mx2 <- matrix(as.numeric(unlist(mx)), byrow = TRUE, nr = nrow(a))
get_indices <- function(x){seq(as.numeric(x[1]), as.numeric(x[2]), 1)}
if (SUBSET == 0){
    SUBSET <- paste(1, "-", ncol(mx2), sep="")
}
subsetIndices <- unlist(lapply(strsplit(unlist(strsplit(SUBSET, split=",", fixed=TRUE)), split="-", fixed=TRUE), get_indices))
cr <- cor(mx2[,subsetIndices])
ids <- is.na(cr)
cr[ids] <- 0
if(sum(ids) > 0){
	print("Presence of some NA values in correlation matrix calculation!! Check your annotation files!! Setting those correlations to 0!!")
}
e <- eigen(cr)
e$values[which(e$values<0)] <- 0
Meff <- (sum(sqrt(e$values)))^2/sum(e$values)
Padj <- 0.05/Meff

write.table(rbind(Meff,Padj), file = OUTPUT, sep="\t",row.names=TRUE,col.names=FALSE, quote=FALSE, append=FALSE)
