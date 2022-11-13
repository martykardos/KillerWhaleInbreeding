# identify ancestral and derived alleles in killer whales
setwd("~/Documents/orca/analyses_31May2022/RXY")

cetGenos <- read.table("cetaceans.var.alleles.14April2021",header=TRUE)

library(data.table)
kwGenos <- fread("kw_151.snp.final.passOnly.biallelicOnly.noIndels.mac1_indivFiltered.hetDepthFiltered_maxMiss0.25.tped")
kwGenos <- as.matrix(kwGenos)
kwIDs <- read.table("kw_151.snp.final.passOnly.biallelicOnly.noIndels.mac1_indivFiltered.hetDepthFiltered_maxMiss0.25.tfam")



#############################################################################
# identify the ancestral allele as the majority among the non-killer whale
# dolphins
##############################################################################



# count occurrence of first allele
whtSideAlls <- as.character(cetGenos[,7])
rightWhaleAlls <- as.character(cetGenos[,11])
indoPacAlls <- as.character(cetGenos[,17])
dolphAlls <- cbind(whtSideAlls,rightWhaleAlls,indoPacAlls)

for(i in 1:3){
  dolphAlls[dolphAlls[,i] == ".",i] <- NA
}

whtSideCounts <- rowSums(dolphAlls[,2:3] == cbind(dolphAlls[,1],dolphAlls[,1]),na.rm=TRUE)
rightWhaleCounts <- rowSums(dolphAlls[,c(1,3)] == cbind(dolphAlls[,2],dolphAlls[,2]),na.rm=TRUE)
indoPacCounts <- rowSums(dolphAlls[,c(1,2)] == cbind(dolphAlls[,3],dolphAlls[,3]),na.rm=TRUE)

ancAll <- rep(NA,nrow(dolphAlls))
ancAll[which(whtSideCounts == 1 | whtSideCounts == 2)] <- dolphAlls[which(whtSideCounts == 1 | whtSideCounts == 2),1]
ancAll[which(rightWhaleCounts == 1 | rightWhaleCounts == 2)] <- dolphAlls[which(rightWhaleCounts == 1 | rightWhaleCounts == 2),1]
ancAll[which(indoPacCounts == 1 | indoPacCounts == 2)] <- dolphAlls[which(indoPacCounts == 1 | indoPacCounts == 2),1]


outAncDat <- cbind(cetGenos[,1:3],ancAll)

# remove indels
remRows <- which(nchar(cetGenos[,3]) > 1) # remove indels
outAncDat <- outAncDat[-remRows,]

# remove non-polarized loci

remRows <- which(is.na(outAncDat[,4]))
outAncDat <- outAncDat[-remRows,]


####################################################################
# link derived allele information to the empirical genotypes
####################################################################
ancOutNames <- paste(outAncDat[,1],":",outAncDat[,2],sep="")

outAncDat <- cbind(ancOutNames,outAncDat)
# remove duplicate ancOutDat records
dups <- which(duplicated(outAncDat[,1]) == TRUE)
outAncDat <- outAncDat[-dups,]

rm(dups)

# make missing genotypes NA

for(i in 5:ncol(kwGenos)){
  if(sum(kwGenos[,i] == "0") > 0){
    kwGenos[which(kwGenos[,i] == "0"),i] <- NA
  }
  print(i)
}

#-----------------------------------------------------
# get rid of loci that are not polarized from kwGenos
#-----------------------------------------------------
remLoci <- which(kwGenos[,2] %in% ancOutNames == FALSE)
kwGenos <- kwGenos[-remLoci,]


# remove ancestral allele loci not present in kwGenos
remAncDat <- which(outAncDat[,1] %in% kwGenos[,2] == FALSE)
outAncDat <- outAncDat[-remAncDat,]

write.table(outAncDat,file="ancestralAlleles_Killers",quote=FALSE,row.names=FALSE)


# collect the garbage
gc()


################################################################
# identify the derived allele(s) and remove any loci where
# all present alleles are derived
#################################################################
ancAllVec <- as.character(outAncDat[,5])
derAllsVec <- rep(NA,nrow(kwGenos))
colIter <- 5
while((sum(is.na(derAllsVec) == FALSE) < length(derAllsVec)) & (colIter<=ncol(kwGenos)) ) {
  if(sum(kwGenos[,colIter] != ancAllVec,na.rm=TRUE) > 0){
    derAllsVec[which((kwGenos[,colIter] != ancAllVec) & (is.na(kwGenos[,colIter]) == FALSE) )] <- kwGenos[which((kwGenos[,colIter] != ancAllVec) & (is.na(kwGenos[,colIter]) == FALSE) ),colIter]
  }
  colIter <- colIter + 1
  print(colIter)
}

#--------------------------------------------
# remove loci fixed for the ancestral allele
#--------------------------------------------
remLoci <- which(is.na(derAllsVec))
derAllsVec <- derAllsVec[-remLoci]
kwGenos <- kwGenos[-remLoci,]
rm(outAncDat)
rm(ancAllVec)

################################################################
# make a matrix of derived allele indicator genotypes
################################################################
derAllMat <- matrix(rep(derAllsVec,294),nrow=nrow(kwGenos),ncol=294)
rm(derAllsVec)
gc()

derGenos <- kwGenos[,5:ncol(kwGenos)] == derAllMat
# make missing genotypes missing
for(i in 1:ncol(derGenos)){
  if(sum(is.na(kwGenos[,i+4])) > 0){
    derGenos[is.na(kwGenos[,i+4]),i] <- NA
  }
  print(i)
}
derGenos <- derGenos * 1                # convert genotypes to numeric format

# write the derived allele matrix
write.table(kwGenos[,1:4],file="kw_DerivedAlleleGenotypeInformationMatrix_1August2022",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(derGenos,file="kw_DerivedAlleleGenotypes_1August2022",quote=FALSE,row.names=FALSE,col.names=FALSE)



