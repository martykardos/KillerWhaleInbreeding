####
#### R script implementing the LOD method of Pemberton et al., (2012, Am J Hum Genet), 
#### Wang et al., (2009, Genet Epidemiol), Kardos et al. (2017, Genetics), Kardos et al. (2018, Nature Ecology & Evolution) to identify runs of homozygosity (for ConGen 2018)
####
####
setwd("~/Documents/orca/analyses_31May2022/ROH3")
library(reshape2)
library(data.table)
##########################################################
# preliminaries: fill these out before running the script
##########################################################
allChroms <- as.character(c(1:5,7:22))                       
winSize <- 100                                                # number of SNPS to use in a single window for IBD LOD score
minSNP <- 50                                                  # minimum umber of non-missing SNPs to include a window LOD score 
stepSize <- 10                                                # the step size (in SNPs) you want to use for sliding the window of size winSize across the genome
e <- 0.00001                                                  # probability of a heterozygous genotype having the highest likelihood under IBD: this can be due to mutation or a sequencing error
minMaf <- 0.05                                                # minimum minor allele frequency

# get allele frequencies
allFreqs <- read.table("kw_151.snp.final.passOnly.biallelicOnly.noIndels.mac1_indivFiltered.hetDepthFiltered.frq",header=TRUE,row.names=NULL)
colnames(allFreqs) <- c("chrom","pos","nAlleles","genomes","allAFreq","allBFreq")
allAFreqs <- colsplit(allFreqs[,5],pattern=":",names=c("allele","freq"))
allFreqMat <- cbind(allFreqs[,1:2],allAFreqs[,2])

# get minor allele frequency
maf <- allFreqMat[,3]
maf[which(maf > 0.5)] <- 1 - maf[which(maf > 0.5)]
removeLoci <- paste(allFreqMat[which(maf <= minMaf),1],"_",allFreqMat[which(maf <= minMaf),2],sep="")

#output name extension
outName <- paste("IBD","_",winSize,"window_",minSNP,"minSNPs_",stepSize,"stepSize_31December2020",sep="")   # name for the output files
chrLengs <- rep(NA,length(allChroms))
for(zb in 1:length(allChroms)){
  LODMat <- NULL            # matrix to store the LOD scores for each window for each individual
  HetMat <- NULL            # matrix to store window heterozygosity
  nSNPMat <- NULL           # matrix to store the number of typed SNPs for each individual in each window
  
  #######################################################
  # get the genotype likelihoods and allele frequencies
  #######################################################
  thisChrom <- allChroms[zb]
  library(data.table)
  genoLike <- read.table(paste("orca_genoLikes_chr",thisChrom,sep=""), header=TRUE)    # read in the genotypes
  genoLike <- as.matrix(genoLike)
  genoLikePos <- colsplit(as.character(genoLike[,1]),pattern="_",names=c("chrom","pos"))[,2]
  markerSplit <- colsplit(string=as.character(genoLike[,1]),pattern="_",names=c("chrom","pos"))
  chromVec <- markerSplit[,1]                        # vectors of the chromosome names and locus positions
  posVec <- as.numeric(markerSplit[,2])
  freqs <- allFreqMat[allFreqMat[,1] == paste("chr",allChroms[zb],sep=""),]
  freq1 <- as.numeric(freqs[,3])

  # remove loci with maf < minMaf
  theseRemLoci <- which(genoLike[,1] %in% removeLoci)
  genoLike <- genoLike[-theseRemLoci,]
  genoLikePos <- genoLikePos[-theseRemLoci]
  posVec <- posVec[-theseRemLoci]
  freqs <- freqs[-theseRemLoci,]
  freq1 <- freq1[-theseRemLoci]
  NumLoci <- nrow(genoLike)
  chrLengs[zb] <- max(genoLikePos)
  print(NumLoci)    # print the number of loci in the data set
  
  numInds <- (length(genoLike[1,]) - 3)/3         # count individuals in the data set
  ##########################################################################################################
  # build vectors of the probability of the data at each locus
  # assuming non-IBD and assuming IBDStatus... equations are from Wang et al. (2009, Genetic Epidemiology), and 
  # 
  ##########################################################################################################
  # note that "a" stands for first of the two observed alleles (arbitrarily), and r stands for the second of the two  
  aaAuto <- NULL # genotype probs assuming autozygosity
  arAuto <- NULL
  rrAuto <- NULL
  missAuto <- NULL
  
  aaNAuto <- NULL # genotype probs assuming non-autozygosity
  arNAuto <- NULL
  rrNAuto <- NULL
  missNAuto <- 1
  
  aaAuto <- freq1^2+freq1*(1-freq1)    # calculate genotype probabilities under IBD
  arAuto <- rep(e,length(aaAuto))
  rrAuto <- (1-freq1)^2+freq1*(1-freq1)
  
  aaNAuto <- freq1^2  ### non-IBD
  arNAuto <- 2*(freq1)*(1-freq1)
  rrNAuto <- (1-freq1)^2
  missNAuto <- 1
  
  #### for each individual, build two vectors of data probalities: one assuming autozygosity, and one assuming non-autozyosity
  ids <- 1:numInds      # arbitrary individual IDs
  autoProbs <- NULL     # object to store the autozygous genotype probs
  nAutoProbs <- NULL    # object to store the non-autozygous probs
  indColStarts <- seq(4,ncol(genoLike),3)
  for(i in 1:length(ids)){
    thisInd <- NULL
    thisInd <- matrix(as.numeric(genoLike[,indColStarts[i]:(indColStarts[i] + 2)]),ncol = 3)
    thisAutoProbs <- rowSums(cbind(aaAuto,arAuto,rrAuto)*thisInd)        
    thisNAutoProbs <- rowSums(cbind(aaNAuto,arNAuto,rrNAuto)*thisInd)
    autoProbs <- cbind(autoProbs,thisAutoProbs)
    nAutoProbs <- cbind(nAutoProbs,thisNAutoProbs)
    print(i)
  }
  
  ############################################################################################################
  # split the genome into windows of SNPs and calculate the LOD score for each window for each individual
  ############################################################################################################
  
  #### get the number of SNPs on each chromosome
  nSNPs <- sum(chromVec == thisChrom)
  winVec <- NULL            # vector to store the ID of the IBD window for each snp
  saveChromVec <- NULL      # vector to store the chromosome ID for each window
  startVec <- NULL          # vectors to store the ID of the start and end SNPs of the IBD windows
  endVec <- NULL            # vector to store the number of missing genotypes for each 60 SNP window
  startPosVec <- NULL       # vector to store the bp position of the start of each window
  endPosVec <- NULL         # vector to store the bp position of the end of each window
  winIterator <- 0          # iterator for windows
  
  for(i in 1:length(ids)){
    thisGenoLikes <- genoLike[,indColStarts[i]:(indColStarts[i] + 2)]
    thisAutoProbs<- NULL
    thisAutoProbs <- autoProbs[,i]
    thisNAutoProbs<- NULL
    thisNAutoProbs <- nAutoProbs[,i]
    LODVec <- NULL
    HetVec <- NULL
    nSNPVec <- NULL
    thisMissVec <- NULL
    chrGenos <- NULL
    chrGenos <- thisGenoLikes
    chrAutoProbs <- thisAutoProbs
    chrNAutoProbs <- thisNAutoProbs
    starts <- NULL
    starts <- seq(1,length(chrAutoProbs),stepSize)     # start SNP of each IBD window in the genome
    starts <- starts[-which(starts > (length(chrAutoProbs)-winSize+1))]   # trim the end of starts so that there are no windows with fewer then winSize SNPs
    ends <- starts + winSize - 1
    if(sum(ends > length(chrAutoProbs)) >0){ends[which(ends > length(chrAutoProbs))] <- length(chrAutoProbs)}
    chrSNPPos <- NULL
    startPos <- NULL 
    endPos <- NULL
    chrSNPPos <- posVec
    startPos <- chrSNPPos[starts]
    endPos <- chrSNPPos[ends]
    
    #-------------------------------
    # save the window information
    #-------------------------------
    if(i == 1){
      saveChromVec <- append(saveChromVec,rep(thisChrom,length(starts)))
      startVec <- append(startVec,starts)
      endVec <- append(endVec,ends)
      startPosVec <- append(startPosVec,startPos)
      endPosVec <- append(endPosVec,endPos)
      theseWins <- NULL
      winStartMat <- as.data.frame(matrix(rep(1:length(starts),winSize),nrow=winSize,ncol=length(starts),byrow=TRUE))
      theseWins <- unlist(winStartMat)
      if(length(theseWins) > length(chrAutoProbs)){theseWins <- theseWins[1:length(chrAutoProbs)]}
      theseWins <- theseWins + winIterator  
    }
    
    #--------------------------------------------------------
    # get the LOD scores and heterozygosity with each window
    #--------------------------------------------------------
    chrLODVec <- rep(NA,length(starts))   # vector to store the LOD score for each window
    chrHetVec <- rep(NA,length(starts))   # vector to store proportion het snps within each window
    chrNSNPVec <- rep(NA,length(starts))
    for (j in 1:length(starts)){
      thisMissVec <- append(thisMissVec,sum(is.na(chrGenos[starts[j]:ends[j],1])))
      chrNSNPVec[j] <- sum(is.na(chrGenos[starts[j]:ends[j],1]) == FALSE)
      
      if(chrNSNPVec[j] >= minSNP){	    
        #########################
        ####### get LOD score
        #########################
        aProbs <- NULL
        aProbs <- chrAutoProbs[starts[j]:ends[j]]
        naProbs <- NULL
        naProbs <- chrNAutoProbs[starts[j]:ends[j]]
        chrLODVec [j] <- sum( log10( aProbs/naProbs),na.rm=TRUE)
        
        ### get proportion of SNPs that are het in this window
        chrHetVec[j] <- sum(as.numeric(chrGenos[starts[j]:ends[j],2]) > as.numeric(chrGenos[starts[j]:ends[j],1]) & as.numeric(chrGenos[starts[j]:ends[j],2]) > as.numeric(chrGenos[starts[j]:ends[j],3]),na.rm=TRUE)/sum(is.na(chrGenos[starts[j]:ends[j],2]) == FALSE)
      }    
    }
    LODVec <- append(LODVec,chrLODVec)
    HetVec <- append(HetVec,chrHetVec)
    nSNPVec <- append(nSNPVec,chrNSNPVec)
    LODMat <- cbind(LODMat,LODVec)
    HetMat <- cbind(HetMat,HetVec)
    nSNPMat <- cbind(nSNPMat,nSNPVec)
    print(paste("done with individual ",i,sep=""))
  }
  
  ###################################################################################################################################
  # Identify the LOD score threshold above which you will call a segment IBD using gaussian kernel density model...
  # This is the local mimimum density between the two modes in the distribution of LOD scores... This follows Pemberton et al., 2012
  ###################################################################################################################################
  
  lodVec <- NULL
  for(i in 1:ncol(LODMat)){
    lodVec <- append(lodVec,LODMat[,i])
  }
  lodVec2 <- lodVec
  if(sum(is.na(lodVec)) > 0){lodVec2 <- lodVec[-which(is.na(lodVec) == TRUE)]}
  dens <- density(lodVec2,kernel="gaussian",na.rm=TRUE)
  plot(dens)
  thresh <- 0
  IBDMat <- LODMat >= thresh
  IBDMat <- as.data.frame(IBDMat)
  winMat <- cbind(saveChromVec,startVec,endVec,startPosVec,endPosVec)    # matrix storing information on the order and location of each IBD window in the genome
  winMat <- as.data.frame(winMat)
  
  ####################################################################################################################
  # concatenate contiguous IBD chromosome segments and save the locations of IBD tracts for each individual 
  ####################################################################################################################
  
  IBDTractMat <- NULL    # store the start, stop positions, and the length of each RoH for each individudal
  for (i in 1:ncol(IBDMat)){
   thisDat <- NULL       # the ith individual's IBD statuses in each window across the genome
   thisDat <- IBDMat[,i]
   outDat <- NULL  # matrix to store the IBD tract information for this individual
   if(sum(is.na(thisDat)) > 0){thisDat[is.na(thisDat)] <- FALSE}
    begins <- NULL          # the first windows in each RoH
    ends <- NULL            # the last window in each RoH
    begins <- which(thisDat[2:length(thisDat)] == TRUE & thisDat[1:(length(thisDat)-1)] == FALSE )+1  # which windows are TRUE (for IBD) and the previous one is FALSE (for IBD)
    if(is.na(thisDat[1]) == FALSE) {
      if(thisDat[1] == TRUE){
        begins <- c(1,begins)
      }
    }
    ends <- which(thisDat[1:(length(thisDat)-1)] == TRUE & thisDat[2:length(thisDat)] == FALSE ) # which windows are true for IBD and the next one is FALSE (for IBD)
    if(is.na(thisDat[length(thisDat)]) == FALSE) {
      if(thisDat[length(thisDat)] == TRUE){
        ends <- c(ends,length(thisDat))
      }
    }
    tractStartPos <- NULL
    tractEndPos <- NULL
    tractStartPos <- as.character(winMat$startPosVec [begins])
    tractEndPos <- as.character(winMat$endPosVec[ends])
    
    #-------------------------------------------------------------------------------------
    # if there are multiple ROH, see if any of them overlap and concatenate them
    #-------------------------------------------------------------------------------------
    newStartPos <- NULL
    newEndPos <- NULL
    if(length(tractStartPos) > 1){
      overs <- NULL   # test if each RoH window overlaps with the previous one
      overs <- as.numeric(tractStartPos[2:length(tractStartPos)]) <= as.numeric(tractEndPos[1:(length(tractStartPos)-1)])
      overs <- c(NA,overs)
      groups <- c(1,rep(NA,length(overs)-1))       # vector of the identity of unique RoH
      for(k in 2:length(groups)){
        if(is.na(overs[k]) == FALSE){
          if(overs[k] == FALSE){groups[k] <- 1 + groups[k-1]}else
            if(overs[k] == TRUE) {groups[k] <- groups[k-1]}
        }
      }
      uniqGroups <- NULL  # unique RoH identifiers
      uniqGroups <- unique(groups)
      for(k in 1:length(uniqGroups)){
        newStartPos <- append(newStartPos, tractStartPos[which(groups == uniqGroups[k])[1]])
        newEndPos   <- append(newEndPos,tractEndPos[which(groups == uniqGroups[k])[length(which(groups == uniqGroups[k]))]])
      }
    }
      
    if(length(tractStartPos) == 1){
      newStartPos <- NULL
      newEndPos <- NULL
      newStartPos <- tractStartPos
      newEndPos <- tractEndPos
    }
    outDat <- rbind(outDat,cbind(rep(ids[i],length(newStartPos)),rep(thisChrom,length(newStartPos)),newStartPos,newEndPos))
    print(i)
    IBDTractMat <- rbind(IBDTractMat,outDat)
  }
  
  colnames(HetMat) <- paste("ind",ids,sep="")
  colnames(IBDTractMat)[1:2] <- c("id","chromosome")
  lengths <- as.numeric(IBDTractMat[,4]) - as.numeric(IBDTractMat[,3])
  
  #########################################################
  # save files
  #########################################################
  write.table(IBDMat,file=paste(outName,"_IBDStatusMatrix_",winSize,"_chrom_",allChroms[zb],sep=""),quote=FALSE,row.names=FALSE)
  write.table(winMat,file=paste(outName,"_windowInformationMatrix_",winSize,"_chrom_",allChroms[zb],sep=""),quote=FALSE,row.names=FALSE)
  write.table(IBDTractMat,file=paste(outName,"_ibdTracts_",winSize,"_chrom_",allChroms[zb],sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)
  write.table(LODMat,file=paste(outName,"_windowLODScoreMatrix_",winSize,"_chrom_",allChroms[zb],sep=""),quote=FALSE,row.names=FALSE)
  write.table(nSNPMat,file=paste(outName,"_nSNPMatrix_",winSize,"_chrom_",allChroms[zb],sep=""),quote=FALSE,row.names=FALSE)
  write.table(HetMat,file=paste(outName,"_hetMat_",winSize,"_chrom_",allChroms[zb],sep=""),quote=FALSE,row.names=FALSE)
}

###############
# make id key
###############
outIDs <- colnames(genoLike[,seq(4,ncol(genoLike),3)])
outIDs2 <- rep(NA,length(outIDs))
for(i in 1:length(outIDs2)){
  thisLength <- nchar(outIDs[i])
  outIDs2[i] <- substr(outIDs[i],start=1,stop=thisLength-2)
}
idNums <- 1:length(outIDs)
outKey <- cbind(as.character(idNums),outIDs2)
colnames(outKey) <- c("idNum","ID")
write.table(outKey,file="orcaIDForROHAnalysis",quote=FALSE,row.names=FALSE)

# save chromosome lengths
outLengs <- cbind(allChroms,chrLengs)
write.table(outLengs,file="orcaChrLengs",quote=FALSE,row.names=FALSE)


