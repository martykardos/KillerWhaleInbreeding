# calculate Froh
setwd("~/Documents/orca/analyses_31May2022/ROH3")

minLengs <- c(1000000,10000000)
chrLengs <- read.table("orcaChrLengs",header=TRUE)
rohFiles <- paste("IBD_100window_50minSNPs_10stepSize_31December2020_ibdTracts_100_chrom_",c(1:5,7:22),sep="")
popKey <- read.table("popKey.txt",header=TRUE) # population key
IDs <- read.table("orcaIDForROHAnalysis",header=TRUE)

tracts <- NULL
for(i in 1:length(rohFiles)){
  theseTracts <- read.table(rohFiles[i],header=FALSE)
  tracts <- rbind(tracts,theseTracts)
  print(i)
}

lengs <- tracts[,4] - tracts[,3]
tracts <- cbind(tracts,lengs)

rohMat <- NULL

frohMat <- NULL
for(i in 1:nrow(IDs)){
  outFroh  <- c(NA,NA,rep(0,2))
  theseTracts <- tracts[which(tracts[,1] == IDs[i,1]),]
  for(j in 1:2){
    if(sum(theseTracts[,5] >= minLengs[j]) > 0){
      lengTracts <- theseTracts [which(theseTracts[,5] >= minLengs[j]),]
      outFroh[j+2] <- sum(lengTracts[,5])/sum(as.numeric(chrLengs[,2]))
    }
  }
  id <- IDs[i,2]
  pop <- popKey[which(popKey[,2] == IDs[i,2]),1]
  outFroh[1] <- id
  outFroh[2] <- pop
  frohMat <- rbind(frohMat,outFroh)
  
  idVec <- rep(id,nrow(theseTracts))
  chrVec <- paste("chr",theseTracts[,2],sep="")
  
  thisOutROH <- cbind(idVec,chrVec,theseTracts[,3:5])
  colnames(thisOutROH) <- c("id","chrom","start","stop","length")
  rohMat <- rbind(rohMat,thisOutROH)
  # save ROH
  
  print(i)
}


write.table(frohMat,file="Froh_7August2022",quote=FALSE,row.names=FALSE)
write.table(rohMat,file="allROH_18August2022",col.names=TRUE,row.names=FALSE)


