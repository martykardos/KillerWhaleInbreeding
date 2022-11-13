setwd("/Users/martin.kardos/Documents/orca/analyses_31May2022/RXY")
kwIDs <- read.table("kw_151.snp.final.passOnly.biallelicOnly.noIndels.mac1_indivFiltered.hetDepthFiltered_maxMiss0.25.tfam")
library(data.table)
popIDs <- read.table("orca_pops_ids.txt",header=TRUE)      # link ids to populations

setwd("~/Documents/orca/analyses_31May2022/Ne/GONE/lowMissingIndivs")
####################
# srkws
####################
ped <- fread("kw_151.snp.final.passOnly.biallelicOnly.noIndels.mac1_indivFiltered.hetDepthFiltered_minGQ20_maxMiss0.25_autosomesOnly_srkwOnly_lowMissing_mac2WithinPop.tped")
ped <- as.matrix(ped)

ped2 <- ped[,5:ncol(ped)]
for(i in 1:ncol(ped2)){
  ped2[which(ped2[,i] == "0"),i] <- "NO"
  print(i)
}

refAlls <- ped2[,1]
i <- 1
while(sum(refAlls == "NO") > 0){
  refAlls[which(refAlls == "NO")] <- ped2[which(refAlls == "NO"),i]
  i <- i+ 1
  print(i)
}

outGenos <- matrix("NO",nrow=nrow(ped2),ncol=ncol(ped2))
for(i in 1:ncol(outGenos)){
  outGenos[which(ped2[,i] == refAlls & ped2[,i] != "NO"),i] <- "0"
  outGenos[which(ped2[,i] != refAlls & ped2[,i] != "NO"),i] <- "1"
  print(i)
}



info <- ped[,1:4]

freqs <- rowSums((outGenos == "1" & outGenos != "NO"),na.rm=TRUE)/(rowSums((outGenos == "NO") == FALSE))
missingness <- rowSums(outGenos == "NO")/ncol(outGenos)
singles <- rowSums(outGenos == "1") == 1
keepLoci <- which(freqs > 0 & freqs < 1 & missingness < 0.1 & info[,1] != "0" & singles == FALSE)
outDat <- outGenos[keepLoci,]
outInfo <- info[keepLoci,]
outNumLoci <- 10000
for(i in 1:reps){
  theseLoci <- sort(sample(1:nrow(outDat),outNumLoci,replace=FALSE))
  thisInInfo <- outInfo[theseLoci,]
  thisInDat <- outDat[theseLoci,]
  thisOutInfo <- thisInInfo[,1:2]
  thisOutInfo[,1] <- paste(rep("c",nrow(thisOutInfo)),thisOutInfo[,1],sep="")
  
  # convert to Genepop
  
  # genePopGenos
  fstAlls <- thisInDat[,seq(1,ncol(thisInDat),2)]
  secAlls <- thisInDat[,seq(2,ncol(thisInDat),2)]
  pasteGenos <- matrix(NA,nrow=ncol(fstAlls),ncol=nrow(fstAlls))
  for(j in 1:nrow(pasteGenos)){
    pasteGenos[j,] <- paste(as.character(fstAlls[,j]),as.character(secAlls[,j]),sep="")
  }
  genePopGenos <- matrix(NA,nrow=nrow(pasteGenos),ncol=ncol(pasteGenos) )
  for(j in 1:nrow(genePopGenos)){
    if(sum(pasteGenos[j,] == "00") > 0) genePopGenos[j,which(pasteGenos[j,] == "00")] <- "0101" 
    if(sum(pasteGenos[j,] == "01") > 0) genePopGenos[j,which(pasteGenos[j,] == "01")] <- "0102" 
    if(sum(pasteGenos[j,] == "10") > 0) genePopGenos[j,which(pasteGenos[j,] == "10")] <- "0201"
    if(sum(pasteGenos[j,] == "11") > 0) genePopGenos[j,which(pasteGenos[j,] == "11")] <- "0202"
    if(sum(pasteGenos[j,] == "NONO") > 0) genePopGenos[j,which(pasteGenos[j,] == "NONO")] <- "0000"
  }
  
  outMat <- rep(NA,nrow(pasteGenos)+3)
  outMat[1] <- "snpData"
  outMat[2] <- paste(paste(thisOutInfo[,2],",",sep=""),collapse=" ")
  outMat[3] <- "pop"
  for(j in 4:length(outMat)){
    outMat[j] <- paste(c(paste("ind_",j,",",sep=""),genePopGenos[j-3,]),collapse=" ")
  }
  writeLines(outMat,con=paste("srkwGenePop_",i,".GEN",sep=""))
  write.table(thisOutInfo,file=paste("srkwGenePop_linkage_",i,sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
}


####################
# arkws
####################
ped <- fread("kw_151.snp.final.passOnly.biallelicOnly.noIndels.mac1_indivFiltered.hetDepthFiltered_minGQ20_maxMiss0.25_autosomesOnly_arkwOnly_lowMissing_mac2WithinPop.tped")
ped <- as.matrix(ped)

ped2 <- ped[,5:ncol(ped)]
for(i in 1:ncol(ped2)){
  ped2[which(ped2[,i] == "0"),i] <- "NO"
  print(i)
}

refAlls <- ped2[,1]
i <- 1
while(sum(refAlls == "NO") > 0){
  refAlls[which(refAlls == "NO")] <- ped2[which(refAlls == "NO"),i]
  i <- i+ 1
  print(i)
}

outGenos <- matrix("NO",nrow=nrow(ped2),ncol=ncol(ped2))
for(i in 1:ncol(outGenos)){
  outGenos[which(ped2[,i] == refAlls & ped2[,i] != "NO"),i] <- "0"
  outGenos[which(ped2[,i] != refAlls & ped2[,i] != "NO"),i] <- "1"
  print(i)
}



info <- ped[,1:4]

freqs <- rowSums((outGenos == "1" & outGenos != "NO"),na.rm=TRUE)/(rowSums((outGenos == "NO") == FALSE))
missingness <- rowSums(outGenos == "NO")/ncol(outGenos)
singles <- rowSums(outGenos == "1") == 1
keepLoci <- which(freqs > 0 & freqs < 1 & missingness < 0.1 & info[,1] != "0" & singles == FALSE)
outDat <- outGenos[keepLoci,]
outInfo <- info[keepLoci,]
outNumLoci <- 10000
for(i in 1:reps){
  theseLoci <- sort(sample(1:nrow(outDat),outNumLoci,replace=FALSE))
  thisInInfo <- outInfo[theseLoci,]
  thisInDat <- outDat[theseLoci,]
  thisOutInfo <- thisInInfo[,1:2]
  thisOutInfo[,1] <- paste(rep("c",nrow(thisOutInfo)),thisOutInfo[,1],sep="")
  
  # convert to Genepop
  
  # genePopGenos
  fstAlls <- thisInDat[,seq(1,ncol(thisInDat),2)]
  secAlls <- thisInDat[,seq(2,ncol(thisInDat),2)]
  pasteGenos <- matrix(NA,nrow=ncol(fstAlls),ncol=nrow(fstAlls))
  for(j in 1:nrow(pasteGenos)){
    pasteGenos[j,] <- paste(as.character(fstAlls[,j]),as.character(secAlls[,j]),sep="")
  }
  genePopGenos <- matrix(NA,nrow=nrow(pasteGenos),ncol=ncol(pasteGenos) )
  for(j in 1:nrow(genePopGenos)){
    if(sum(pasteGenos[j,] == "00") > 0) genePopGenos[j,which(pasteGenos[j,] == "00")] <- "0101" 
    if(sum(pasteGenos[j,] == "01") > 0) genePopGenos[j,which(pasteGenos[j,] == "01")] <- "0102" 
    if(sum(pasteGenos[j,] == "10") > 0) genePopGenos[j,which(pasteGenos[j,] == "10")] <- "0201"
    if(sum(pasteGenos[j,] == "11") > 0) genePopGenos[j,which(pasteGenos[j,] == "11")] <- "0202"
    if(sum(pasteGenos[j,] == "NONO") > 0) genePopGenos[j,which(pasteGenos[j,] == "NONO")] <- "0000"
  }
  
  outMat <- rep(NA,nrow(pasteGenos)+3)
  outMat[1] <- "snpData"
  outMat[2] <- paste(paste(thisOutInfo[,2],",",sep=""),collapse=" ")
  outMat[3] <- "pop"
  for(j in 4:length(outMat)){
    outMat[j] <- paste(c(paste("ind_",j,",",sep=""),genePopGenos[j-3,]),collapse=" ")
  }
  writeLines(outMat,con=paste("arkwGenePop_",i,".GEN",sep=""))
  write.table(thisOutInfo,file=paste("arkwGenePop_linkage_",i,sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
}






####################
# tkws
####################
ped <- fread("kw_151.snp.final.passOnly.biallelicOnly.noIndels.mac1_indivFiltered.hetDepthFiltered_minGQ20_maxMiss0.25_autosomesOnly_tkwOnly_lowMissing_mac2WithinPop.tped")
ped <- as.matrix(ped)

ped2 <- ped[,5:ncol(ped)]
for(i in 1:ncol(ped2)){
  ped2[which(ped2[,i] == "0"),i] <- "NO"
  print(i)
}

remLoci <- which(rowSums(ped2 == "NO") == ncol(ped2))

if(length(remLoci) > 0){
  ped <- ped[-remLoci,]
  ped2 <- ped2[-remLoci,]
}


refAlls <- ped2[,1]
i <- 1
while(sum(refAlls == "NO") > 0 & (i <= ncol(ped2)) ){
  refAlls[which(refAlls == "NO")] <- ped2[which(refAlls == "NO"),i]
  i <- i+ 1
  print(i)
}

outGenos <- matrix("NO",nrow=nrow(ped2),ncol=ncol(ped2))
for(i in 1:ncol(outGenos)){
  outGenos[which(ped2[,i] == refAlls & ped2[,i] != "NO"),i] <- "0"
  outGenos[which(ped2[,i] != refAlls & ped2[,i] != "NO"),i] <- "1"
  print(i)
}



info <- ped[,1:4]

freqs <- rowSums((outGenos == "1" & outGenos != "NO"),na.rm=TRUE)/(rowSums((outGenos == "NO") == FALSE))
missingness <- rowSums(outGenos == "NO")/ncol(outGenos)
singles <- rowSums(outGenos == "1") == 1
keepLoci <- which(freqs > 0 & freqs < 1 & missingness < 0.1 & info[,1] != "0" & singles == FALSE)
outDat <- outGenos[keepLoci,]
outInfo <- info[keepLoci,]
outNumLoci <- 10000
for(i in 1:reps){
  theseLoci <- sort(sample(1:nrow(outDat),outNumLoci,replace=FALSE))
  thisInInfo <- outInfo[theseLoci,]
  thisInDat <- outDat[theseLoci,]
  thisOutInfo <- thisInInfo[,1:2]
  thisOutInfo[,1] <- paste(rep("c",nrow(thisOutInfo)),thisOutInfo[,1],sep="")
  
  # convert to Genepop
  
  # genePopGenos
  fstAlls <- thisInDat[,seq(1,ncol(thisInDat),2)]
  secAlls <- thisInDat[,seq(2,ncol(thisInDat),2)]
  pasteGenos <- matrix(NA,nrow=ncol(fstAlls),ncol=nrow(fstAlls))
  for(j in 1:nrow(pasteGenos)){
    pasteGenos[j,] <- paste(as.character(fstAlls[,j]),as.character(secAlls[,j]),sep="")
  }
  genePopGenos <- matrix(NA,nrow=nrow(pasteGenos),ncol=ncol(pasteGenos) )
  for(j in 1:nrow(genePopGenos)){
    if(sum(pasteGenos[j,] == "00") > 0) genePopGenos[j,which(pasteGenos[j,] == "00")] <- "0101" 
    if(sum(pasteGenos[j,] == "01") > 0) genePopGenos[j,which(pasteGenos[j,] == "01")] <- "0102" 
    if(sum(pasteGenos[j,] == "10") > 0) genePopGenos[j,which(pasteGenos[j,] == "10")] <- "0201"
    if(sum(pasteGenos[j,] == "11") > 0) genePopGenos[j,which(pasteGenos[j,] == "11")] <- "0202"
    if(sum(pasteGenos[j,] == "NONO") > 0) genePopGenos[j,which(pasteGenos[j,] == "NONO")] <- "0000"
  }
  
  outMat <- rep(NA,nrow(pasteGenos)+3)
  outMat[1] <- "snpData"
  outMat[2] <- paste(paste(thisOutInfo[,2],",",sep=""),collapse=" ")
  outMat[3] <- "pop"
  for(j in 4:length(outMat)){
    outMat[j] <- paste(c(paste("ind_",j,",",sep=""),genePopGenos[j-3,]),collapse=" ")
  }
  writeLines(outMat,con=paste("tkwGenePop_",i,".GEN",sep=""))
  write.table(thisOutInfo,file=paste("tkwGenePop_linkage_",i,sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
}





