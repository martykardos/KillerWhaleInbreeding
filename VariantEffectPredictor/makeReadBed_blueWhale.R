# make bed files to make fastq from reference fastas
setwd("/scratch/mkardos/orca/cetaceanGenomes/cactus")
options(scipen=999999)
chroms <- read.table("chromSizes.GCA_009873245.3_mBalMus1.pri.v3_genomic_pruned.fasta",header=FALSE)
readCount <- sum(chroms[,2])/10
readBed <- matrix(NA,nrow=round(readCount + 0.1 *readCount),ncol=3)

rowIter <- 1
for(i in 1:nrow(chroms)){
  chrLeng <- chroms[i,2]
  chrName <- chroms[i,1]
  
  starts <- seq(0,chrLeng-60,10)
  ends <- starts + 70
  
  if(ends[length(ends)] > chrLeng) {
    ends[length(ends)] <- chrLeng
    starts[length(starts)] <- chrLeng - 69
  }
  
  chrBed <- cbind(rep(chrName,length(starts)),starts,ends)
  
  readBed[rowIter:(rowIter + nrow(chrBed) - 1),] <- chrBed
  rowIter <- rowIter + nrow(chrBed)
  print(i)
}
if(sum(readBed[,2] < 0,na.rm=TRUE) > 0)readBed <- readBed[-which(readBed[,2] < 0),]

readBed <- readBed[-which(is.na(readBed[,1])),]
write.table(readBed,file="GCA_009873245.3_mBalMus1.pri.v3_genomic_pruned_read.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")

