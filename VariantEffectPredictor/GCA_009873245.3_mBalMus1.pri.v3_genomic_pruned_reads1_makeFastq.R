inReads <- readLines("GCA_009873245.3_mBalMus1.pri.v3_genomic_pruned_reads1")
inReads <- unlist(inReads)
outVec <- rep(NA,length(inReads)*2)
readInserts <- sort(c(seq(1,length(outVec),4),seq(2,length(outVec),4)))
outVec[readInserts] <- inReads
outVec[seq(3,length(outVec),4)] <- "+"
outVec[seq(4,length(outVec),4)] <- paste(rep("K",70),collapse="")
outMat <- NULL
outMat <- cbind(outMat,outVec)
write.table(outMat,file="GCA_009873245.3_mBalMus1.pri.v3_genomic_pruned_reads1.fastq",quote=FALSE,col.names=FALSE,row.names=FALSE)



