# convert Phred-scaled genotype likelihoods to probabilities for ROH analysis
setwd("~/Documents/orca/analyses_31May2022/ROH3")
library(reshape2)
library(data.table)
library(matrixStats)
PLs <- read.table("kw_151.snp.final.passOnly.biallelicOnly.noIndels.mac1_indivFiltered.hetDepthFiltered.PL",header=FALSE,colClasses="character")
depths <- read.table("kw_151.snp.final.passOnly.biallelicOnly.noIndels.mac1_indivFiltered.hetDepthFiltered.idepth",header=TRUE)
PLs <- as.matrix(PLs)

infoMat <- PLs[,1:4]
colnames(infoMat) <- c("chrom","pos","all1","all2")
marker <- paste(infoMat[,1],"_",infoMat[,2],sep="")
alleles <- infoMat[,3:4]

PLs <- cbind(marker,alleles,PLs[,5:ncol(PLs)])


gc()       # cleanup space

genoProbs <- matrix(NA,nrow=nrow(PLs),ncol=147*3) # store genotype probabilities
fstCols <- seq(1,ncol(genoProbs),3)
for(i in 1:147){             # loop through individuals
  theseVals1 <- as.matrix(colsplit(PLs[,i+3],pattern=",",names=as.character(1:3)))
  if(sum(theseVals1[,1] == ".",na.rm=TRUE) > 0){
    theseVals1[which=(theseVals1[,1] == "."),1] <- NA
  }
  theseVals <- NULL
  for(j in 1:3){
    theseVals <- cbind(theseVals,as.numeric(theseVals1[,j]))
  }
  # denormalize
  normMat <- theseVals
  normMat[which(normMat == 0)] <- NA
  normVals <- rowMins(normMat,na.rm=TRUE)
  normOutMat <- matrix(rep(normVals,3),nrow=nrow(normMat),ncol=3,byrow=FALSE)
  theseVals <- theseVals + normOutMat
  
  #calculate raw probs
  theseProbs <- 10^(theseVals/-10)
  
  # rescale the probs so that they add up to one
  scaler <- matrix(rep(1/(rowSums(theseProbs)),3),nrow=nrow(theseProbs),ncol=3) 
  theseProbs <- scaler*theseProbs
  genoProbs[,(fstCols[i]:(fstCols[i]+2))] <- theseProbs
  print(i)
}

namesVec <- rep(NA,ncol(genoProbs))
starts <- seq(1,length(namesVec),3)
for(i in 1:nrow(depths)){
  namesVec[starts[i]:(starts[i]+2)]<- paste(depths[i,1],".",1:3,sep="")
}
colnames(genoProbs) <- namesVec

outGenoProbs <- cbind(PLs[,1:3],genoProbs)

# separate genotype probs by chromosome and save each of them separately
chroms <- c(1:5,7:22) # exclude chromosome 6 (the sex chromosome)
chromMat1 <- colsplit(outGenoProbs[,1],pattern="_",names=c("col1","col2"))
chromVec <- colsplit(chromMat1[,1],pattern="r",names=c("col1","col2"))[,2]
for(i in 1:length(chroms)){
  thisChromOut <- outGenoProbs[which(chromVec == chroms[i]),]
  write.table(thisChromOut,file=paste("orca_genoLikes_chr",chroms[i],sep=""),quote=FALSE,row.names=FALSE)
  print(i)
}
  

for(i in c(1:5,7:22)){
  print(sum(infoMat[,1] == paste("chr",c(1:5,7:22)[i],sep="")))
}