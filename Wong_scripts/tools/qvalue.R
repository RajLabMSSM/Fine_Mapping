args = commandArgs(trailingOnly=TRUE)

library(qvalue)
myData <- read.table(args[1], header=as.logical(args[2]))
myPvals <- as.numeric(myData[,1])
qobj <- qvalue(p = myPvals)
write.table(qobj$qvalues, quote=FALSE, row.names=FALSE, col.names=FALSE)

