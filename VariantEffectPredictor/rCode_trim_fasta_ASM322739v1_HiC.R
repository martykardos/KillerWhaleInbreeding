inFastaNamePrefix <- "ASM322739v1_HiC"
thisFasta <- readLines(paste(inFastaNamePrefix,".fasta",sep=""))
nameIndex <- grep(">",thisFasta)
seqStarts <- nameIndex + 1
seqEnds <- nameIndex[2] - 1
seqEnds <- c(seqEnds,nameIndex[2:length(nameIndex)] - 1,length(thisFasta))
seqLeng <- nchar(thisFasta[2])

##############################################
# identify locations in thisFasta with Ns
##############################################

NLocs <- grep("N",thisFasta)
if(sum(NLocs %in% nameIndex) > 0) NLocs <- NLocs[-which(NLocs %in% nameIndex == TRUE)]

###################################################
# remove all sequences that are completely masked
###################################################
library(stringr)
NCounts <- str_count(thisFasta[NLocs],"N")
if(sum(NCounts == seqLeng) > 0){
  remSeqs <- NLocs[which(NCounts == seqLeng)]
  thisFasta <- thisFasta[-remSeqs]
  
  #################################################
  # revise sequence starts and ends
  #################################################
  nameIndex <- grep(">",thisFasta)
  seqStarts <- nameIndex + 1
  seqEnds <- nameIndex[2] - 1
  seqEnds <- c(seqEnds,nameIndex[3:length(nameIndex)] - 1,length(thisFasta))
  
  # remove singleton read scaffolds if necessary
  #if(sum(nameIndex[2:length(nameIndex)] - nameIndex[1:(length(nameIndex) - 1)] <= 2) > 0){
  #  thisFasta <- thisFasta[ -nameIndex [which(nameIndex[1:(length(nameIndex) - 1)] == nameIndex[2:length(nameIndex)] -1)]]
  #  nameIndex <- grep(">",thisFasta)
  #  seqStarts <- nameIndex + 1
  #  seqEnds <- nameIndex[2] - 1
  #  seqEnds <- c(seqEnds,nameIndex[2:length(nameIndex)] - 1,length(thisFasta))
  #}
  NLocs <- grep("N",thisFasta)
  if(sum(NLocs %in% nameIndex) > 0) NLocs <- NLocs[-which(NLocs %in% nameIndex == TRUE)]
}

# identify scaffolds that need trimming of Ns
outFasta <- rep(NA,length(thisFasta)*2)
ticker <- 1

for(i in 1:length(nameIndex)){
  if(sum(NLocs >= seqStarts[i] & NLocs <= seqEnds[i]) == 0) outFasta[ticker:(ticker + length(thisFasta[nameIndex[i]:seqEnds[i]])-1)] <- thisFasta[nameIndex[i]:seqEnds[i]]
  if(sum(NLocs >= seqStarts[i] & NLocs <= seqEnds[i]) > 0){
    thisInScaf <- paste(thisFasta[seqStarts[i]:seqEnds[i]],collapse="")
    thisOutScaf <- NULL
    # break the scaffold at any occurrence of Ns
    if(gregexpr("N",thisInScaf)[[1]][1] != -1){
      theseNLocs <- gregexpr("N",thisInScaf)[[1]]
      theseBaseLocs <- (1:nchar(thisInScaf))[-theseNLocs]
      
      firstNLocs <- theseNLocs[1]
      if(sum(theseNLocs[1:(length(theseNLocs)-1)] != theseNLocs[2:(length(theseNLocs))] -1) > 0 & length(theseNLocs) > 1){
        firstNLocs <- c(firstNLocs,theseNLocs[which(theseNLocs[1:(length(theseNLocs)-1)] != theseNLocs[2:(length(theseNLocs))] -1) +1])
        lastNLocs <- theseNLocs[which(theseNLocs[2:length(theseNLocs)] != theseNLocs[1:(length(theseNLocs)-1)] +1)]
        if(sum(theseNLocs[2:length(theseNLocs)] != theseNLocs[1:(length(theseNLocs)-1)] +1) > 0)lastNLocs<- c(lastNLocs,theseNLocs[length(theseNLocs)])
      }
      if( length(theseNLocs) == 1){
        lastNLocs <- theseNLocs[1]
      }
      firstBaseLocs <- NULL # vectors of the start and stop of base positions of non masked DNA
      lastBaseLocs <- NULL
      
      if(firstNLocs[1] != 1){  # if the first N block is not at the beginning of the scaffold
        for(j in 1:length(firstNLocs)){
          if(j == 1){
            firstBaseLocs <- 1
            lastBaseLocs <- firstNLocs[j] - 1
          }
          if(j > 1) {
            firstBaseLocs <- c(firstBaseLocs,lastNLocs[j - 1] + 1)
            lastBaseLocs <- c(lastBaseLocs,firstNLocs[j] - 1)
          }
        }
        if(lastNLocs[length(lastNLocs)] != nchar(thisInScaf)){
          firstBaseLocs <- c(firstBaseLocs,lastNLocs[length(lastNLocs)] + 1)
          lastBaseLocs <- c(lastBaseLocs,nchar(thisInScaf))
        }
      }
      
      if(firstNLocs[1] == 1){  # if the first N block IS at the beginning of the scaffold
        for(j in 1:length(firstNLocs)){
          if(j < length(firstNLocs)){
            firstBaseLocs <- c(firstBaseLocs,lastNLocs[j - 1] + 1)
            lastBaseLocs <- c(lastBaseLocs,firstNLocs[j] - 1)
          }
          
          if(j == length(firstNLocs)){
            if(lastNLocs[j] < nchar(thisInScaf)){
              firstBaseLocs <- c(firstBaseLocs,lastNLocs[j - 1] + 1)
              lastBaseLocs <- c(lastBaseLocs,firstNLocs[j] - 1)
            }
          }
        }
        if(lastNLocs[length(lastNLocs)] != nchar(thisInScaf)){
          firstBaseLocs <- c(firstBaseLocs,lastNLocs[length(lastNLocs)] + 1)
          lastBaseLocs <- c(lastBaseLocs,nchar(thisInScaf))
        }
      }
    }
    outScafNames <- paste(thisFasta[nameIndex[i]],"_",1:length(firstBaseLocs),sep="")
    theseOutScafs <- NULL
    for(j in 1:length(outScafNames)){
      bases <- toupper(substring(thisInScaf,firstBaseLocs[j],lastBaseLocs[j]))
      basesLeng <- nchar(bases)
      if(basesLeng > seqLeng){
        readStarts <- seq(1,basesLeng,seqLeng)
        readEnds <- seq(seqLeng,basesLeng,seqLeng)
        if(length(readEnds) < length(readStarts))readEnds <- c(readEnds,basesLeng)
        reads <- substring(bases,first=readStarts,last=readEnds)
      }
      theseOutScafs <- c(theseOutScafs,outScafNames[j],reads)
    }
    outFasta[ticker:(ticker + length(theseOutScafs) - 1)] <- theseOutScafs
  }
  ticker <- min(which(is.na(outFasta)))
  print(i)
}
outFasta <- outFasta[-which(is.na(outFasta))]
finalFasta <- NULL
finalFasta <- cbind(finalFasta,outFasta)
write.table(finalFasta,file=paste(inFastaNamePrefix,"_pruned.fasta",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)

