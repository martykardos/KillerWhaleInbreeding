


setwd("/Users/martin.kardos/Documents/orca/analyses_31May2022/RXY")
kwIDs <- read.table("kw_151.snp.final.passOnly.biallelicOnly.noIndels.mac1_indivFiltered.hetDepthFiltered_maxMiss0.25.tfam")
library(data.table)
derGenos <- fread("kw_DerivedAlleleGenotypes_1August2022",colClasses=rep("numeric",294))         # derived/ancestral genotypes (0 = ancestral, 1 = derived)
derGenoInfo <- fread("kw_DerivedAlleleGenotypeInformationMatrix_1August2022",colClasses=rep("character",4))     # chromosome and position of each locus in the genome
popIDs <- read.table("orca_pops_ids.txt",header=TRUE)      # link ids to populations
varEffs <- read.table("variant_effect_output_skinny")
# get IDs to keep in the analysis

#IDs_keep <- read.table("orcaIDForROHAnalysis",header=TRUE)



#######################################################
# organize the genotypes by population
#######################################################
derGenos <- as.matrix(derGenos)
derGenoInfo <- as.matrix(derGenoInfo)
srkwIDs <- popIDs[which(popIDs[,1] == "SRKW"),2]
arkwIDs <- popIDs[which(popIDs[,1] == "ARKW"),2]
tkwIDs <- popIDs[which(popIDs[,1] == "TKW"),2]

offshoreIDs <- popIDs[which(popIDs[,1] == "offshore"),2]
nrkwIDs <-  popIDs[which(popIDs[,1] == "NRKW"),2]

idVec <- NULL
for(i in 1:nrow(kwIDs)){
  idVec <- c(idVec,rep(as.character(kwIDs[i,1]),2))
}

colnames(derGenoInfo) <- c("chrom","name","blank","position")
colnames(derGenos) <- idVec


srkwGenos <- derGenos[,which(colnames(derGenos) %in% srkwIDs)]
arkwGenos <- derGenos[,which(colnames(derGenos) %in% arkwIDs)]
tkwGenos <- derGenos[,which(colnames(derGenos) %in% tkwIDs)]
nrkwGenos <- derGenos[,which(colnames(derGenos) %in% nrkwIDs)]
offshoreGenos <- derGenos[,which(colnames(derGenos) %in% offshoreIDs)]

############################################################################
# calculate MLH for each individual in each popualtion
# using only loci that are polymorphic within that particular population
############################################################################
chrLengs <- read.table("orcaChrLengs",header=TRUE)
genSize <- sum(chrLengs[,2]) # genome size is 

#arkw
theseGenos <- arkwGenos
allFreqs <- rowSums(theseGenos,na.rm=TRUE)/rowSums(is.na(theseGenos) == FALSE)
polys <- which(allFreqs > 0 & allFreqs < 1)
polyGenos <- theseGenos[polys,]
fstCols <- seq(1,ncol(theseGenos),2)
secCols <- fstCols + 1
hetMat <- theseGenos[,fstCols] != theseGenos[,secCols]
nonMiss <- colSums(is.na(hetMat) == FALSE)

hetMatPoly <- polyGenos[,fstCols] != polyGenos[,secCols]       # only include polymorphic loci as usual 
indHet <- NULL
lociHet <- c(5000,10000,15000,20000,50000,nrow(hetMatPoly))      # numbers of loci to use to calculate individual heterozygosity

for(i in 1:length(lociHet)){
  theseLoci <- sample(1:nrow(hetMatPoly),size=lociHet[i],replace=FALSE)             # randomly select polymorphic loci
  nonMissTheseLoci <- colSums(is.na(hetMatPoly[theseLoci,]) == FALSE)             # count missing genotypes for each individual
  theseHets <- colSums(hetMatPoly[theseLoci,],na.rm=TRUE)/nonMissTheseLoci # save heterozygosity
  indHet<- cbind(indHet,theseHets)
}

arkwIndHet <- indHet
arkwIndHet <- cbind(colnames(theseGenos)[seq(1,ncol(theseGenos),2)],arkwIndHet)   # attache names to heterozygosity
arkwHet <- mean( ((colSums(hetMat,na.rm=TRUE)/nonMiss)*nrow(hetMat))/genSize )


chrLengs <- read.table("orcaChrLengs",header=TRUE)
genSize <- sum(chrLengs[,2]) # genome size is 

#srkw
theseGenos <- srkwGenos
allFreqs <- rowSums(theseGenos,na.rm=TRUE)/rowSums(is.na(theseGenos) == FALSE)
polys <- which(allFreqs > 0 & allFreqs < 1)
polyGenos <- theseGenos[polys,]
fstCols <- seq(1,ncol(theseGenos),2)
secCols <- fstCols + 1
hetMat <- theseGenos[,fstCols] != theseGenos[,secCols]
nonMiss <- colSums(is.na(hetMat) == FALSE)

hetMatPoly <- polyGenos[,fstCols] != polyGenos[,secCols]       # only include polymorphic loci as usual 
indHet <- NULL
lociHet <- c(5000,10000,15000,20000,50000,nrow(hetMatPoly))      # numbers of loci to use to calculate individual heterozygosity

for(i in 1:length(lociHet)){
  theseLoci <- sample(1:nrow(hetMatPoly),size=lociHet[i],replace=FALSE)             # randomly select polymorphic loci
  nonMissTheseLoci <- colSums(is.na(hetMatPoly[theseLoci,]) == FALSE)             # count missing genotypes for each individual
  theseHets <- colSums(hetMatPoly[theseLoci,],na.rm=TRUE)/nonMissTheseLoci # save heterozygosity
  indHet<- cbind(indHet,theseHets)
}

srkwIndHet <- indHet
srkwIndHet <- cbind(colnames(theseGenos)[seq(1,ncol(theseGenos),2)],srkwIndHet)   # attache names to heterozygosity
srkwHet <- mean( ((colSums(hetMat,na.rm=TRUE)/nonMiss)*nrow(hetMat))/genSize )

chrLengs <- read.table("orcaChrLengs",header=TRUE)
genSize <- sum(chrLengs[,2]) # genome size is 

#tkw
theseGenos <- tkwGenos
allFreqs <- rowSums(theseGenos,na.rm=TRUE)/rowSums(is.na(theseGenos) == FALSE)
polys <- which(allFreqs > 0 & allFreqs < 1)
polyGenos <- theseGenos[polys,]
fstCols <- seq(1,ncol(theseGenos),2)
secCols <- fstCols + 1
hetMat <- theseGenos[,fstCols] != theseGenos[,secCols]
nonMiss <- colSums(is.na(hetMat) == FALSE)

hetMatPoly <- polyGenos[,fstCols] != polyGenos[,secCols]       # only include polymorphic loci as usual 
indHet <- NULL
lociHet <- c(5000,10000,15000,20000,50000,nrow(hetMatPoly))      # numbers of loci to use to calculate individual heterozygosity

for(i in 1:length(lociHet)){
  theseLoci <- sample(1:nrow(hetMatPoly),size=lociHet[i],replace=FALSE)             # randomly select polymorphic loci
  nonMissTheseLoci <- colSums(is.na(hetMatPoly[theseLoci,]) == FALSE)             # count missing genotypes for each individual
  theseHets <- colSums(hetMatPoly[theseLoci,],na.rm=TRUE)/nonMissTheseLoci # save heterozygosity
  indHet<- cbind(indHet,theseHets)
}

tkwIndHet <- indHet
tkwIndHet <- cbind(colnames(theseGenos)[seq(1,ncol(theseGenos),2)],tkwIndHet)   # attache names to heterozygosity
tkwHet <- mean( ((colSums(hetMat,na.rm=TRUE)/nonMiss)*nrow(hetMat))/genSize )


chrLengs <- read.table("orcaChrLengs",header=TRUE)
genSize <- sum(chrLengs[,2]) # genome size is 

#offshore
theseGenos <- offshoreGenos
allFreqs <- rowSums(theseGenos,na.rm=TRUE)/rowSums(is.na(theseGenos) == FALSE)
polys <- which(allFreqs > 0 & allFreqs < 1)
polyGenos <- theseGenos[polys,]
fstCols <- seq(1,ncol(theseGenos),2)
secCols <- fstCols + 1
hetMat <- theseGenos[,fstCols] != theseGenos[,secCols]
nonMiss <- colSums(is.na(hetMat) == FALSE)

hetMatPoly <- polyGenos[,fstCols] != polyGenos[,secCols]       # only include polymorphic loci as usual 
indHet <- NULL
lociHet <- c(5000,10000,15000,20000,50000,nrow(hetMatPoly))      # numbers of loci to use to calculate individual heterozygosity

for(i in 1:length(lociHet)){
  theseLoci <- sample(1:nrow(hetMatPoly),size=lociHet[i],replace=FALSE)             # randomly select polymorphic loci
  nonMissTheseLoci <- colSums(is.na(hetMatPoly[theseLoci,]) == FALSE)             # count missing genotypes for each individual
  theseHets <- colSums(hetMatPoly[theseLoci,],na.rm=TRUE)/nonMissTheseLoci # save heterozygosity
  indHet<- cbind(indHet,theseHets)
}

offshoreIndHet <- indHet
offshoreIndHet <- cbind(colnames(theseGenos)[seq(1,ncol(theseGenos),2)],offshoreIndHet)   # attache names to heterozygosity
offshoreHet <- mean( ((colSums(hetMat,na.rm=TRUE)/nonMiss)*nrow(hetMat))/genSize )

#tkw
theseGenos <- nrkwGenos
allFreqs <- rowSums(theseGenos,na.rm=TRUE)/rowSums(is.na(theseGenos) == FALSE)
polys <- which(allFreqs > 0 & allFreqs < 1)
polyGenos <- theseGenos[polys,]
fstCols <- seq(1,ncol(theseGenos),2)
secCols <- fstCols + 1
hetMat <- theseGenos[,fstCols] != theseGenos[,secCols]
nonMiss <- colSums(is.na(hetMat) == FALSE)

hetMatPoly <- polyGenos[,fstCols] != polyGenos[,secCols]       # only include polymorphic loci as usual 
indHet <- NULL
lociHet <- c(5000,10000,15000,20000,50000,nrow(hetMatPoly))      # numbers of loci to use to calculate individual heterozygosity

for(i in 1:length(lociHet)){
  theseLoci <- sample(1:nrow(hetMatPoly),size=lociHet[i],replace=FALSE)             # randomly select polymorphic loci
  nonMissTheseLoci <- colSums(is.na(hetMatPoly[theseLoci,]) == FALSE)             # count missing genotypes for each individual
  theseHets <- colSums(hetMatPoly[theseLoci,],na.rm=TRUE)/nonMissTheseLoci # save heterozygosity
  indHet<- cbind(indHet,theseHets)
}

nrkwIndHet <- indHet
nrkwIndHet <- cbind(colnames(theseGenos)[seq(1,ncol(theseGenos),2)],nrkwIndHet)   # attache names to heterozygosity
nrkwHet <- mean( ((colSums(hetMat,na.rm=TRUE)/nonMiss)*nrow(hetMat))/genSize )


chrLengs <- read.table("orcaChrLengs",header=TRUE)
genSize <- sum(chrLengs[,2]) # genome size is 

colnames(srkwIndHet) <- c("ID",paste(as.character(lociHet[1:(length(lociHet)-1)]),"_loci",sep=""),"allLoci")
colnames(arkwIndHet) <- c("ID",paste(as.character(lociHet[1:(length(lociHet)-1)]),"_loci",sep=""),"allLoci")
colnames(nrkwIndHet) <- c("ID",paste(as.character(lociHet[1:(length(lociHet)-1)]),"_loci",sep=""),"allLoci")
colnames(tkwIndHet) <- c("ID",paste(as.character(lociHet[1:(length(lociHet)-1)]),"_loci",sep=""),"allLoci")
colnames(offshoreIndHet) <- c("ID",paste(as.character(lociHet[1:(length(lociHet)-1)]),"_loci",sep=""),"allLoci")

write.table(srkwIndHet,file="srkw_indivHet",row.names=FALSE,quote=FALSE)
write.table(arkwIndHet,file="arkw_indivHet",row.names=FALSE,quote=FALSE)
write.table(nrkwIndHet,file="nrkw_indivHet",row.names=FALSE,quote=FALSE)
write.table(tkwIndHet,file="tkw_indivHet",row.names=FALSE,quote=FALSE)
write.table(offshoreIndHet,file="offshore_indivHet",row.names=FALSE,quote=FALSE)


srkwOutHet <- cbind(srkwIDs,srkwIndHet)
arkwOutHet <- cbind(arkwIDs,arkwIndHet)
nrkwOutHet <- cbind(nrkwIDs,nrkwIndHet)
tkwOutHet <- cbind(tkwIDs,tkwIndHet)
offshoreOutHet <- cbind(offshoreIDs,offshoreIndHet)

par(mfrow=c(2,3))


#######################################################
# organize the variant effects
#######################################################



# remove annotation information for loci not found in the genotypes
#varEffs <- cbind(varEffs,delVec)
varEffs <- varEffs[-which(varEffs[,2] %in% derGenoInfo[,2] == FALSE),]

# remove genotypes not found in the annotation information
remGenos <- which(derGenoInfo[,2] %in% varEffs[,2] == FALSE)

derGenoInfo <- derGenoInfo[-remGenos,]
srkwGenos <- srkwGenos[-remGenos,]
arkwGenos <- arkwGenos[-remGenos,]
tkwGenos <- tkwGenos[-remGenos,]
nrkwGenos <- nrkwGenos[-remGenos,]
offshoreGenos <- offshoreGenos[-remGenos,]


varEffVec <- rep("0",nrow(varEffs))
varEffVec[grep("HIGH",varEffs[,5])] <- "HIGH"
varEffVec[grep("MODERATE",varEffs[,5])] <- "MOD"
varEffVec[grep("intergenic",varEffs[,4])] <- "NEUT"


outEffVec <- varEffVec[match(derGenoInfo[,2],varEffs[,2])]

derGenoInfo <- cbind(derGenoInfo,outEffVec)

######################################################
# clean up the workspace
######################################################
rm(derGenos)
gc()

###############################################################
# get the SFS for ARKW
###############################################################

# arkw first
delLoci <- arkwGenos[derGenoInfo[,5] == "HIGH",]   # extract deleterious loci
modLoci <- arkwGenos[derGenoInfo[,5] == "MOD",]
neutLoci <- arkwGenos[derGenoInfo[,5] == "NEUT",]  # extract the neutral loci

# get missingness for each type of locus
delMiss <- rowSums(is.na(delLoci))/2
hist(delMiss[which(delMiss < 10)])

neutMiss <- rowSums(is.na(neutLoci))/2
hist(neutMiss[which(neutMiss < 10)])

modMiss <- rowSums(is.na(modLoci))/2
hist(modMiss[which(modMiss < 10)])

sampInds <- 12         # number of individuals to downsample to

delLoci <- delLoci[-which(delMiss >= 2),]        # remove loci with too many missing genotypes
neutLoci <- neutLoci[-which(neutMiss >= 2),]
modLoci <- modLoci[-which(modMiss >= 2),]

#-------------------------------------
# downsample genotypes
#-------------------------------------

# highly deleterious variants
newDelLoci <- matrix(NA,nrow=nrow(delLoci),ncol=2*sampInds)
fstCols <- seq(1,ncol(delLoci),2)
for(i in 1:nrow(newDelLoci)){
  grabCols <- sample(which(is.na(delLoci[i,]) == FALSE),sampInds*2)
  newDelLoci[i,] <- delLoci[i,grabCols]
  print(i)
}

delDerCounts <- rowSums(newDelLoci)
delDerCounts[which(delDerCounts == 0)] <- NA
delDerCounts[which(delDerCounts == ncol(newDelLoci))] <- NA


# netural variants
newNeutLoci <- matrix(NA,nrow=nrow(neutLoci),ncol=2*sampInds)
fstCols <- seq(1,ncol(neutLoci),2)
for(i in 1:nrow(newNeutLoci)){
  grabCols <- sample(which(is.na(neutLoci[i,]) == FALSE),sampInds*2)
  newNeutLoci[i,] <- neutLoci[i,grabCols]
  print(i)
}
neutDerCounts <- rowSums(newNeutLoci)
neutDerCounts[which(neutDerCounts == 0)] <- NA
neutDerCounts[which(neutDerCounts == ncol(newNeutLoci))] <- NA


# moderate effect variants
newModLoci <- matrix(NA,nrow=nrow(modLoci),ncol=2*sampInds)
fstCols <- seq(1,ncol(modLoci),2)
for(i in 1:nrow(newModLoci)){
  grabCols <- sample(which(is.na(modLoci[i,]) == FALSE),sampInds*2)
  newModLoci[i,] <- modLoci[i,grabCols]
  print(i)
}

modDerCounts <- rowSums(newModLoci)
modDerCounts[which(modDerCounts == 0)] <- NA
modDerCounts[which(modDerCounts == ncol(newModLoci))] <- NA



neutSFS <- rep(NA,length(min(neutDerCounts,na.rm=TRUE):max(neutDerCounts,na.rm=TRUE)))
delSFS <- rep(NA,length(min(delDerCounts,na.rm=TRUE):max(delDerCounts,na.rm=TRUE)))
modSFS <- rep(NA,length(min(modDerCounts,na.rm=TRUE):max(modDerCounts,na.rm=TRUE)))
for(i in 1:length(neutSFS)){
  neutSFS[i] <- sum(neutDerCounts == i,na.rm=TRUE)
  delSFS[i] <- sum(delDerCounts == i,na.rm=TRUE)
  modSFS[i] <- sum(modDerCounts == i,na.rm=TRUE)

}


arkw_neutSFS <- neutSFS
arkw_delSFS <- delSFS
arkw_modSFS <- modSFS
allCount <- 1:length(arkw_neutSFS)
sfsARKW <- cbind(allCount,arkw_neutSFS,arkw_delSFS,arkw_modSFS)


###############################################################
# get the SFS for SRKW
###############################################################

delLoci <- srkwGenos[derGenoInfo[,5] == "HIGH",]   # extract delterious loci
modLoci <- srkwGenos[derGenoInfo[,5] == "MOD",]
neutLoci <- srkwGenos[derGenoInfo[,5] == "NEUT",]  # extract the neutral loci

# get missingness for each type of locus
delMiss <- rowSums(is.na(delLoci))/2
hist(delMiss[which(delMiss < 10)])

neutMiss <- rowSums(is.na(neutLoci))/2
hist(neutMiss[which(neutMiss < 10)])

modMiss <- rowSums(is.na(modLoci))/2
hist(modMiss[which(modMiss < 10)])

delLoci <- delLoci[-which(delMiss >= 2),]        # remove loci with too many missing genotypes
neutLoci <- neutLoci[-which(neutMiss >= 2),]
modLoci <- modLoci[-which(modMiss >= 2),]

#-------------------------------------
# downsample genotypes
#-------------------------------------

# highly deleterious variants
newDelLoci <- matrix(NA,nrow=nrow(delLoci),ncol=2*sampInds)
fstCols <- seq(1,ncol(delLoci),2)
for(i in 1:nrow(newDelLoci)){
  grabCols <- sample(which(is.na(delLoci[i,]) == FALSE),sampInds*2)
  newDelLoci[i,] <- delLoci[i,grabCols]
  print(i)
}

delDerCounts <- rowSums(newDelLoci)
delDerCounts[which(delDerCounts == 0)] <- NA
delDerCounts[which(delDerCounts == ncol(newDelLoci))] <- NA


# netural variants
newNeutLoci <- matrix(NA,nrow=nrow(neutLoci),ncol=2*sampInds)
fstCols <- seq(1,ncol(neutLoci),2)
for(i in 1:nrow(newNeutLoci)){
  grabCols <- sample(which(is.na(neutLoci[i,]) == FALSE),sampInds*2)
  newNeutLoci[i,] <- neutLoci[i,grabCols]
  print(i)
}
neutDerCounts <- rowSums(newNeutLoci)
neutDerCounts[which(neutDerCounts == 0)] <- NA
neutDerCounts[which(neutDerCounts == ncol(newNeutLoci))] <- NA


# moderate effect variants
newModLoci <- matrix(NA,nrow=nrow(modLoci),ncol=2*sampInds)
fstCols <- seq(1,ncol(modLoci),2)
for(i in 1:nrow(newModLoci)){
  grabCols <- sample(which(is.na(modLoci[i,]) == FALSE),sampInds*2)
  newModLoci[i,] <- modLoci[i,grabCols]
  print(i)
}

modDerCounts <- rowSums(newModLoci)
modDerCounts[which(modDerCounts == 0)] <- NA
modDerCounts[which(modDerCounts == ncol(newModLoci))] <- NA



neutSFS <- rep(NA,length(min(neutDerCounts,na.rm=TRUE):max(neutDerCounts,na.rm=TRUE)))
delSFS <- rep(NA,length(min(delDerCounts,na.rm=TRUE):max(delDerCounts,na.rm=TRUE)))
modSFS <- rep(NA,length(min(modDerCounts,na.rm=TRUE):max(modDerCounts,na.rm=TRUE)))
for(i in 1:length(neutSFS)){
  neutSFS[i] <- sum(neutDerCounts == i,na.rm=TRUE)
  delSFS[i] <- sum(delDerCounts == i,na.rm=TRUE)
  modSFS[i] <- sum(modDerCounts == i,na.rm=TRUE)

}


srkw_neutSFS <- neutSFS
srkw_delSFS <- delSFS
srkw_modSFS <- modSFS
allCount <- 1:length(srkw_neutSFS)
sfsSRKW <- cbind(allCount,srkw_neutSFS,srkw_delSFS,srkw_modSFS)


###############################################################
# get the SFS for TKW
###############################################################

delLoci <- tkwGenos[derGenoInfo[,5] == "HIGH",]   # extract delterious loci
modLoci <- tkwGenos[derGenoInfo[,5] == "MOD",]
neutLoci <- tkwGenos[derGenoInfo[,5] == "NEUT",]  # extract the neutral loci

# get missingness for each type of locus
delMiss <- rowSums(is.na(delLoci))/2
hist(delMiss[which(delMiss < 10)])

neutMiss <- rowSums(is.na(neutLoci))/2
hist(neutMiss[which(neutMiss < 10)])

modMiss <- rowSums(is.na(modLoci))/2
hist(modMiss[which(modMiss < 10)])


delLoci <- delLoci[-which(delMiss >= 1),]        # remove loci with too many missing genotypes
neutLoci <- neutLoci[-which(neutMiss >= 1),]
modLoci <- modLoci[-which(modMiss >= 1),]

#-------------------------------------
# downsample genotypes
#-------------------------------------

# highly deleterious variants
newDelLoci <- matrix(NA,nrow=nrow(delLoci),ncol=2*sampInds)
fstCols <- seq(1,ncol(delLoci),2)
for(i in 1:nrow(newDelLoci)){
  grabCols <- sample(which(is.na(delLoci[i,]) == FALSE),sampInds*2)
  newDelLoci[i,] <- delLoci[i,grabCols]
  print(i)
}

delDerCounts <- rowSums(newDelLoci)
delDerCounts[which(delDerCounts == 0)] <- NA
delDerCounts[which(delDerCounts == ncol(newDelLoci))] <- NA


# netural variants
newNeutLoci <- matrix(NA,nrow=nrow(neutLoci),ncol=2*sampInds)
fstCols <- seq(1,ncol(neutLoci),2)
for(i in 1:nrow(newNeutLoci)){
  grabCols <- sample(which(is.na(neutLoci[i,]) == FALSE),sampInds*2)
  newNeutLoci[i,] <- neutLoci[i,grabCols]
  print(i)
}
neutDerCounts <- rowSums(newNeutLoci)
neutDerCounts[which(neutDerCounts == 0)] <- NA
neutDerCounts[which(neutDerCounts == ncol(newNeutLoci))] <- NA


# moderate effect variants
newModLoci <- matrix(NA,nrow=nrow(modLoci),ncol=2*sampInds)
fstCols <- seq(1,ncol(modLoci),2)
for(i in 1:nrow(newModLoci)){
  grabCols <- sample(which(is.na(modLoci[i,]) == FALSE),sampInds*2)
  newModLoci[i,] <- modLoci[i,grabCols]
  print(i)
}

modDerCounts <- rowSums(newModLoci)
modDerCounts[which(modDerCounts == 0)] <- NA
modDerCounts[which(modDerCounts == ncol(newModLoci))] <- NA



neutSFS <- rep(NA,length(min(neutDerCounts,na.rm=TRUE):max(neutDerCounts,na.rm=TRUE)))
delSFS <- rep(NA,length(min(delDerCounts,na.rm=TRUE):max(delDerCounts,na.rm=TRUE)))
modSFS <- rep(NA,length(min(modDerCounts,na.rm=TRUE):max(modDerCounts,na.rm=TRUE)))
for(i in 1:length(neutSFS)){
  neutSFS[i] <- sum(neutDerCounts == i,na.rm=TRUE)
  delSFS[i] <- sum(delDerCounts == i,na.rm=TRUE)
  modSFS[i] <- sum(modDerCounts == i,na.rm=TRUE)

}


tkw_neutSFS <- neutSFS
tkw_delSFS <- delSFS
tkw_modSFS <- modSFS
allCount <- 1:length(tkw_neutSFS)
sfsTKW <- cbind(allCount,tkw_neutSFS,tkw_delSFS,tkw_modSFS)

###################################################
# write the SFS
###################################################
write.table(sfsSRKW,file="sfs_SRKW_downSamp_12",quote=FALSE,row.names=FALSE)
write.table(sfsTKW,file="sfs_TKW_downSamp_12",quote=FALSE,row.names=FALSE)
write.table(sfsARKW,file="sfs_ARKW_downSamp_12",quote=FALSE,row.names=FALSE)




################################################
# plot the sfs
################################################
sfsSRKW <- read.table("sfs_SRKW_downSamp_12",header=TRUE)
sfsARKW <- read.table("sfs_ARKW_downSamp_12",header=TRUE)
sfsTKW <- read.table("sfs_TKW_downSamp_12",header=TRUE)


par(mfrow=c(1,3))
library(scales)
# SRKW
plot(c(0,1),c(0,0.22),type="n",ylab="Proportion of Loci",xlab="Derived Allele Frequency",main="Southern Residents",cex.lab=1.3)
lines(sfsSRKW[,1]/(12*2),sfsSRKW[,2]/sum(sfsSRKW[,2]),col=alpha("orange",alpha=0.5),lwd=3)
lines(sfsSRKW[,1]/(12*2),(sfsSRKW[,3]+sfsSRKW[,4])/sum(sfsSRKW[,3]+sfsSRKW[,4]),col=alpha("blue",alpha=0.5),lwd=3)

legend(x=0.4,y=0.17,xjust=FALSE,yjust=FALSE,legend=c("Deleterious","Neutral"),lty="solid",
       lwd=2,col=c(alpha("blue",alpha=0.5),alpha("orange",alpha=0.5)),bty="n")


# ARKW
plot(c(0,1),c(0,0.22),type="n",ylab="",xlab="Derived Allele Frequency",main="Alaska Residents",cex.lab=1.3)
lines(sfsARKW[,1]/(12*2),sfsARKW[,2]/sum(sfsARKW[,2]),col=alpha("orange",alpha=0.5),lwd=3)
lines(sfsARKW[,1]/(12*2),(sfsARKW[,3]+sfsARKW[,4])/sum(sfsARKW[,3]+sfsARKW[,4]),col=alpha("blue",alpha=0.5),lwd=3)

# TKW
plot(c(0,1),c(0,0.22),type="n",ylab="",xlab="Derived Allele Frequency",main="Transients",cex.lab=1.3)
lines(sfsTKW[,1]/(12*2),sfsTKW[,2]/sum(sfsTKW[,2]),col=alpha("orange",alpha=0.5),lwd=3)
lines(sfsTKW[,1]/(12*2),(sfsTKW[,3]+sfsTKW[,4])/sum(sfsTKW[,3]+sfsTKW[,4]),col=alpha("blue",alpha=0.5),lwd=3)




###########################################################
# compare population to population for each type of locus
###########################################################



par(mfrow=c(1,2))
library(scales)

#---------------------
# Neutral SFS
#---------------------
# SRKW
plot(c(0,1),c(0,0.22),type="n",ylab="Proportion of Loci",xlab="Derived Allele Frequency",main="Neutral Variants",cex.lab=1.3)
lines(sfsSRKW[,1]/(12*2),sfsSRKW[,2]/sum(sfsSRKW[,2]),col="#a6cee3",lwd=3)

# ARKW
lines(sfsARKW[,1]/(12*2),sfsARKW[,2]/sum(sfsARKW[,2]),col="#1f78b4",lwd=3)

# TKW
lines(sfsTKW[,1]/(12*2),sfsTKW[,2]/sum(sfsTKW[,2]),col="#b2df8a",lwd=3)


legend(x=0.3,y=0.17,xjust=FALSE,yjust=FALSE,legend=c("Transients","Alaska Residents","Southern Residents"),lty="solid",
       lwd=2,col=c("#b2df8a","#1f78b4","#a6cee3"),bty="n",cex=0.8)



#-----------------------------
#  deleterious SFS
#-----------------------------
# SRKW
plot(c(0,1),c(0,0.22),type="n",ylab="",xlab="Derived Allele Frequency",main="Deleterious Variants",cex.lab=1.3)
lines(sfsSRKW[,1]/(12*2),(sfsSRKW[,3]+sfsSRKW[,4])/sum(sfsSRKW[,3]+sfsSRKW[,4]),col="#a6cee3",lwd=3)

# ARKW
lines(sfsARKW[,1]/(12*2),(sfsARKW[,3]+sfsARKW[,4])/sum(sfsARKW[,3]+sfsARKW[,4]),col="#1f78b4",lwd=3)



# TKW
lines(sfsTKW[,1]/(12*2),(sfsTKW[,3]+sfsTKW[,4])/sum(sfsTKW[,3]+sfsTKW[,4]),col="#b2df8a",lwd=3)


######################################################################
# calculate the number of homozygous deleterious variants for each
# individual
######################################################################

# srkw
delLoci <- srkwGenos[derGenoInfo[,5] == "HIGH" | derGenoInfo[,5] == "MOD",]   # extract delterious loci
ids <- unique(colnames(delLoci))
homs <- rep(NA,length(ids))
for(i in 1:length(homs)){
  thisDat <- delLoci[,colnames(delLoci) == ids[i]]
  homs[i] <- sum(thisDat[,1] == 1 & thisDat[,2] == 1,na.rm=TRUE)
  print(i)
}
outHoms <- cbind(ids,homs)
write.table(outHoms,file="deleteriousHomozygousCount_srkw",quote=FALSE,row.names=FALSE)


# arkw
delLoci <- arkwGenos[derGenoInfo[,5] == "HIGH" | derGenoInfo[,5] == "MOD",]   # extract delterious loci
ids <- unique(colnames(delLoci))
homs <- rep(NA,length(ids))
for(i in 1:length(homs)){
  thisDat <- delLoci[,colnames(delLoci) == ids[i]]
  homs[i] <- sum(thisDat[,1] == 1 & thisDat[,2] == 1,na.rm=TRUE)
  print(i)
}
outHoms <- cbind(ids,homs)
write.table(outHoms,file="deleteriousHomozygousCount_arkw",quote=FALSE,row.names=FALSE)


# tkw
delLoci <- tkwGenos[derGenoInfo[,5] == "HIGH" | derGenoInfo[,5] == "MOD",]   # extract delterious loci
ids <- unique(colnames(delLoci))
homs <- rep(NA,length(ids))
for(i in 1:length(homs)){
  thisDat <- delLoci[,colnames(delLoci) == ids[i]]
  homs[i] <- sum(thisDat[,1] == 1 & thisDat[,2] == 1,na.rm=TRUE)
  print(i)
}
outHoms <- cbind(ids,homs)
write.table(outHoms,file="deleteriousHomozygousCount_tkw",quote=FALSE,row.names=FALSE)


# nrkw
delLoci <- nrkwGenos[derGenoInfo[,5] == "HIGH" | derGenoInfo[,5] == "MOD",]   # extract delterious loci
ids <- unique(colnames(delLoci))
homs <- rep(NA,length(ids))
for(i in 1:length(homs)){
  thisDat <- delLoci[,colnames(delLoci) == ids[i]]
  homs[i] <- sum(thisDat[,1] == 1 & thisDat[,2] == 1,na.rm=TRUE)
  print(i)
}
outHoms <- cbind(ids,homs)
write.table(outHoms,file="deleteriousHomozygousCount_nrkw",quote=FALSE,row.names=FALSE)


# offshore
delLoci <- offshoreGenos[derGenoInfo[,5] == "HIGH" | derGenoInfo[,5] == "MOD",]   # extract delterious loci
ids <- unique(colnames(delLoci))
homs <- rep(NA,length(ids))
for(i in 1:length(homs)){
  thisDat <- delLoci[,colnames(delLoci) == ids[i]]
  homs[i] <- sum(thisDat[,1] == 1 & thisDat[,2] == 1,na.rm=TRUE)
  print(i)
}
outHoms <- cbind(ids,homs)
write.table(outHoms,file="deleteriousHomozygousCount_offshore",quote=FALSE,row.names=FALSE)





######################################################################
# calculate the number of homozygous LOF variants for each
# individual
######################################################################

# srkw
delLoci <- srkwGenos[derGenoInfo[,5] == "HIGH",]   # extract delterious loci
ids <- unique(colnames(delLoci))
homs <- rep(NA,length(ids))
for(i in 1:length(homs)){
  thisDat <- delLoci[,colnames(delLoci) == ids[i]]
  homs[i] <- sum(thisDat[,1] == 1 & thisDat[,2] == 1,na.rm=TRUE)
  print(i)
}
outHoms <- cbind(ids,homs)
write.table(outHoms,file="lofHomozygousCount_srkw",quote=FALSE,row.names=FALSE)


# arkw
delLoci <- arkwGenos[derGenoInfo[,5] == "HIGH",]   # extract delterious loci
ids <- unique(colnames(delLoci))
homs <- rep(NA,length(ids))
for(i in 1:length(homs)){
  thisDat <- delLoci[,colnames(delLoci) == ids[i]]
  homs[i] <- sum(thisDat[,1] == 1 & thisDat[,2] == 1,na.rm=TRUE)
  print(i)
}
outHoms <- cbind(ids,homs)
write.table(outHoms,file="lofHomozygousCount_arkw",quote=FALSE,row.names=FALSE)


# tkw
delLoci <- tkwGenos[derGenoInfo[,5] == "HIGH",]   # extract delterious loci
ids <- unique(colnames(delLoci))
homs <- rep(NA,length(ids))
for(i in 1:length(homs)){
  thisDat <- delLoci[,colnames(delLoci) == ids[i]]
  homs[i] <- sum(thisDat[,1] == 1 & thisDat[,2] == 1,na.rm=TRUE)
  print(i)
}
outHoms <- cbind(ids,homs)
write.table(outHoms,file="lofHomozygousCount_tkw",quote=FALSE,row.names=FALSE)


# nrkw
delLoci <- nrkwGenos[derGenoInfo[,5] == "HIGH",]   # extract delterious loci
ids <- unique(colnames(delLoci))
homs <- rep(NA,length(ids))
for(i in 1:length(homs)){
  thisDat <- delLoci[,colnames(delLoci) == ids[i]]
  homs[i] <- sum(thisDat[,1] == 1 & thisDat[,2] == 1,na.rm=TRUE)
  print(i)
}
outHoms <- cbind(ids,homs)
write.table(outHoms,file="lofHomozygousCount_nrkw",quote=FALSE,row.names=FALSE)


# offshore
delLoci <- offshoreGenos[derGenoInfo[,5] == "HIGH",]   # extract delterious loci
ids <- unique(colnames(delLoci))
homs <- rep(NA,length(ids))
for(i in 1:length(homs)){
  thisDat <- delLoci[,colnames(delLoci) == ids[i]]
  homs[i] <- sum(thisDat[,1] == 1 & thisDat[,2] == 1,na.rm=TRUE)
  print(i)
}
outHoms <- cbind(ids,homs)
write.table(outHoms,file="lofHomozygousCount_offshore",quote=FALSE,row.names=FALSE)




################################################################################
# calculate the number of homozygous moderately deleterious variants for each
# individual
################################################################################

# srkw
delLoci <- srkwGenos[derGenoInfo[,5] == "MOD",]   # extract delterious loci
ids <- unique(colnames(delLoci))
homs <- rep(NA,length(ids))
for(i in 1:length(homs)){
  thisDat <- delLoci[,colnames(delLoci) == ids[i]]
  homs[i] <- sum(thisDat[,1] == 1 & thisDat[,2] == 1,na.rm=TRUE)
  print(i)
}
outHoms <- cbind(ids,homs)
write.table(outHoms,file="moderateHomozygousCount_srkw",quote=FALSE,row.names=FALSE)


# arkw
delLoci <- arkwGenos[derGenoInfo[,5] == "MOD",]   # extract delterious loci
ids <- unique(colnames(delLoci))
homs <- rep(NA,length(ids))
for(i in 1:length(homs)){
  thisDat <- delLoci[,colnames(delLoci) == ids[i]]
  homs[i] <- sum(thisDat[,1] == 1 & thisDat[,2] == 1,na.rm=TRUE)
  print(i)
}
outHoms <- cbind(ids,homs)
write.table(outHoms,file="moderateHomozygousCount_arkw",quote=FALSE,row.names=FALSE)


# tkw
delLoci <- tkwGenos[derGenoInfo[,5] == "MOD",]   # extract delterious loci
ids <- unique(colnames(delLoci))
homs <- rep(NA,length(ids))
for(i in 1:length(homs)){
  thisDat <- delLoci[,colnames(delLoci) == ids[i]]
  homs[i] <- sum(thisDat[,1] == 1 & thisDat[,2] == 1,na.rm=TRUE)
  print(i)
}
outHoms <- cbind(ids,homs)
write.table(outHoms,file="moderateHomozygousCount_tkw",quote=FALSE,row.names=FALSE)


# nrkw
delLoci <- nrkwGenos[derGenoInfo[,5] == "MOD",]   # extract delterious loci
ids <- unique(colnames(delLoci))
homs <- rep(NA,length(ids))
for(i in 1:length(homs)){
  thisDat <- delLoci[,colnames(delLoci) == ids[i]]
  homs[i] <- sum(thisDat[,1] == 1 & thisDat[,2] == 1,na.rm=TRUE)
  print(i)
}
outHoms <- cbind(ids,homs)
write.table(outHoms,file="moderateHomozygousCount_nrkw",quote=FALSE,row.names=FALSE)


# offshore
delLoci <- offshoreGenos[derGenoInfo[,5] == "MOD",]   # extract delterious loci
ids <- unique(colnames(delLoci))
homs <- rep(NA,length(ids))
for(i in 1:length(homs)){
  thisDat <- delLoci[,colnames(delLoci) == ids[i]]
  homs[i] <- sum(thisDat[,1] == 1 & thisDat[,2] == 1,na.rm=TRUE)
  print(i)
}
outHoms <- cbind(ids,homs)
write.table(outHoms,file="moderateHomozygousCount_offshore",quote=FALSE,row.names=FALSE)




############################################################################
#---------------------------------------------------------------------------
# repeat the counting of homozygous deleterious variants for each individual
# after removing fixed loci from the data in each population
#---------------------------------------------------------------------------
############################################################################


######################################################################
# calculate the number of homozygous deleterious variants for each
# individual
######################################################################

fixedVec <- NULL
segVec <- NULL


# srkw
delLoci <- srkwGenos[derGenoInfo[,5] == "HIGH" | derGenoInfo[,5] == "MOD",]   # extract delterious loci

seg <- which(rowSums(delLoci,na.rm=TRUE) < rowSums(is.na(delLoci) == FALSE) & rowSums(delLoci,na.rm=TRUE) > 0)
segVec <- c(segVec,length(seg))
fixed <- which(rowSums(delLoci,na.rm=TRUE) == rowSums(is.na(delLoci) == FALSE))
fixedVec <- c(fixedVec,length(fixed))
delLoci <- delLoci[-fixed,]

ids <- unique(colnames(delLoci))
homs <- rep(NA,length(ids))
for(i in 1:length(homs)){
  thisDat <- delLoci[,colnames(delLoci) == ids[i]]
  homs[i] <- sum(thisDat[,1] == 1 & thisDat[,2] == 1,na.rm=TRUE)
  print(i)
}
outHoms <- cbind(ids,homs)
write.table(outHoms,file="deleteriousHomozygousCount_noFixedLoci_srkw",quote=FALSE,row.names=FALSE)


# arkw
delLoci <- arkwGenos[derGenoInfo[,5] == "HIGH" | derGenoInfo[,5] == "MOD",]   # extract delterious loci
seg <- which(rowSums(delLoci,na.rm=TRUE) < rowSums(is.na(delLoci) == FALSE) & rowSums(delLoci,na.rm=TRUE) > 0)
segVec <- c(segVec,length(seg))
fixed <- which(rowSums(delLoci,na.rm=TRUE) == rowSums(is.na(delLoci) == FALSE))
fixedVec <- c(fixedVec,length(fixed))
delLoci <- delLoci[-fixed,]

ids <- unique(colnames(delLoci))
homs <- rep(NA,length(ids))
for(i in 1:length(homs)){
  thisDat <- delLoci[,colnames(delLoci) == ids[i]]
  homs[i] <- sum(thisDat[,1] == 1 & thisDat[,2] == 1,na.rm=TRUE)
  print(i)
}
outHoms <- cbind(ids,homs)
write.table(outHoms,file="deleteriousHomozygousCount_noFixedLoci_arkw",quote=FALSE,row.names=FALSE)


# tkw
delLoci <- tkwGenos[derGenoInfo[,5] == "HIGH" | derGenoInfo[,5] == "MOD",]   # extract delterious loci
seg <- which(rowSums(delLoci,na.rm=TRUE) < rowSums(is.na(delLoci) == FALSE) & rowSums(delLoci,na.rm=TRUE) > 0)
segVec <- c(segVec,length(seg))
fixed <- which(rowSums(delLoci,na.rm=TRUE) == rowSums(is.na(delLoci) == FALSE))
fixedVec <- c(fixedVec,length(fixed))
delLoci <- delLoci[-fixed,]

ids <- unique(colnames(delLoci))
homs <- rep(NA,length(ids))
for(i in 1:length(homs)){
  thisDat <- delLoci[,colnames(delLoci) == ids[i]]
  homs[i] <- sum(thisDat[,1] == 1 & thisDat[,2] == 1,na.rm=TRUE)
  print(i)
}
outHoms <- cbind(ids,homs)
write.table(outHoms,file="deleteriousHomozygousCount_noFixedLoci_tkw",quote=FALSE,row.names=FALSE)


# nrkw
delLoci <- nrkwGenos[derGenoInfo[,5] == "HIGH" | derGenoInfo[,5] == "MOD",]   # extract delterious loci
seg <- which(rowSums(delLoci,na.rm=TRUE) < rowSums(is.na(delLoci) == FALSE) & rowSums(delLoci,na.rm=TRUE) > 0)
segVec <- c(segVec,length(seg))
fixed <- which(rowSums(delLoci,na.rm=TRUE) == rowSums(is.na(delLoci) == FALSE))
fixedVec <- c(fixedVec,length(fixed))
delLoci <- delLoci[-fixed,]

ids <- unique(colnames(delLoci))
homs <- rep(NA,length(ids))
for(i in 1:length(homs)){
  thisDat <- delLoci[,colnames(delLoci) == ids[i]]
  homs[i] <- sum(thisDat[,1] == 1 & thisDat[,2] == 1,na.rm=TRUE)
  print(i)
}
outHoms <- cbind(ids,homs)
write.table(outHoms,file="deleteriousHomozygousCount_noFixedLoci_nrkw",quote=FALSE,row.names=FALSE)


# offshore
delLoci <- offshoreGenos[derGenoInfo[,5] == "HIGH" | derGenoInfo[,5] == "MOD",]   # extract delterious loci
seg <- which(rowSums(delLoci,na.rm=TRUE) < rowSums(is.na(delLoci) == FALSE) & rowSums(delLoci,na.rm=TRUE) > 0)
segVec <- c(segVec,length(seg))
fixed <- which(rowSums(delLoci,na.rm=TRUE) == rowSums(is.na(delLoci) == FALSE))
fixedVec <- c(fixedVec,length(fixed))
delLoci <- delLoci[-fixed,]

ids <- unique(colnames(delLoci))
homs <- rep(NA,length(ids))
for(i in 1:length(homs)){
  thisDat <- delLoci[,colnames(delLoci) == ids[i]]
  homs[i] <- sum(thisDat[,1] == 1 & thisDat[,2] == 1,na.rm=TRUE)
  print(i)
}
outHoms <- cbind(ids,homs)
write.table(outHoms,file="deleteriousHomozygousCount_noFixedLoci_offshore",quote=FALSE,row.names=FALSE)
popID <- c("SRKW","ARKW","TKW","NRKW","Offshore")

segFixed <- cbind(popID,segVec,fixedVec)
write.table(segFixed,file="numDelFixedAndSegregating_allPops",quote=FALSE,row.names=FALSE)


