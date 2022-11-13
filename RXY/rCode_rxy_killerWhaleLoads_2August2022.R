setwd("/Users/martin.kardos/Documents/orca/analyses_31May2022/RXY")
kwIDs <- read.table("kw_151.snp.final.passOnly.biallelicOnly.noIndels.mac1_indivFiltered.hetDepthFiltered_maxMiss0.25.tfam")
library(data.table)
derGenos <- fread("kw_DerivedAlleleGenotypes_1August2022",colClasses=rep("numeric",294))         # derived/ancestral genotypes (0 = ancestral, 1 = derived)
derGenoInfo <- fread("kw_DerivedAlleleGenotypeInformationMatrix_1August2022",colClasses=rep("character",4))     # chromosome and position of each locus in the genome
popIDs <- read.table("orca_pops_ids.txt",header=TRUE)      # link ids to populations
varEffs <- read.table("variant_effect_output_skinny")
# get IDs to keep in the analysis
IDs_keep <- read.table("orcaIDForROHAnalysis",header=TRUE)


#######################################################
# organize the genotypes by population
#######################################################
derGenos <- as.matrix(derGenos)
derGenoInfo <- as.matrix(derGenoInfo)
srkwIDs <- popIDs[which(popIDs[,1] == "SRKW"),2]
arkwIDs <- popIDs[which(popIDs[,1] == "ARKW"),2]
tkwIDs <- popIDs[which(popIDs[,1] == "TKW"),2]


idVec <- NULL
for(i in 1:nrow(kwIDs)){
  idVec <- c(idVec,rep(as.character(kwIDs[i,1]),2))
}

colnames(derGenoInfo) <- c("chrom","name","blank","position")
colnames(derGenos) <- idVec

srkwGenos <- derGenos[,which(colnames(derGenos) %in% srkwIDs)]
arkwGenos <- derGenos[,which(colnames(derGenos) %in% arkwIDs)]
tkwGenos <- derGenos[,which(colnames(derGenos) %in% tkwIDs)]

#######################################################
# organize the variant effects
#######################################################


# remove genotypes not found in the annotation information
remGenos <- which(derGenoInfo[,2] %in% varEffs[,2] == FALSE)

derGenoInfo <- derGenoInfo[-remGenos,]
srkwGenos <- srkwGenos[-remGenos,]
arkwGenos <- arkwGenos[-remGenos,]
tkwGenos <- tkwGenos[-remGenos,]
# remove annotation information for loci not found in the genotypes
varEffs <- varEffs[-which(varEffs[,2] %in% derGenoInfo[,2] == FALSE),]
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



#########################################
# subsample individual genomes
#########################################
srkwNonMiss <- rowSums(is.na(srkwGenos)== FALSE)
arkwNonMiss <- rowSums(is.na(arkwGenos) == FALSE)
tkwNonMiss <- rowSums(is.na(tkwGenos)== FALSE)
minToKeep <- 24

remLoci <- which(srkwNonMiss < 24 | arkwNonMiss < 24 | tkwNonMiss < 24)

derGenoInfo <- derGenoInfo[-remLoci,]
srkwGenos <- srkwGenos[-remLoci,]
tkwGenos <- tkwGenos[-remLoci,]
arkwGenos <- arkwGenos[-remLoci,]

gc()

###########################
# sub sample genomes
###########################

# SRKW first
newSRKWGenos <- matrix(NA,nrow=nrow(srkwGenos),ncol=minToKeep)
fstCols <- seq(1,ncol(newSRKWGenos),2)
for(i in 1:nrow(newSRKWGenos)){
  grabCols <- sample(which(is.na(srkwGenos[i,]) == FALSE),minToKeep)
  newSRKWGenos[i,] <- srkwGenos[i,grabCols]
  print(i)
}


# ARKW
newARKWGenos <- matrix(NA,nrow=nrow(arkwGenos),ncol=minToKeep)
fstCols <- seq(1,ncol(newARKWGenos),2)
for(i in 1:nrow(newARKWGenos)){
  grabCols <- sample(which(is.na(arkwGenos[i,]) == FALSE),minToKeep)
  newARKWGenos[i,] <- arkwGenos[i,grabCols]
  print(i)
}


# TKW
newTKWGenos <- matrix(NA,nrow=nrow(tkwGenos),ncol=minToKeep)
fstCols <- seq(1,ncol(newTKWGenos),2)
for(i in 1:nrow(newTKWGenos)){
  grabCols <- sample(which(is.na(tkwGenos[i,]) == FALSE),minToKeep)
  newTKWGenos[i,] <- tkwGenos[i,grabCols]
  print(i)
}

arkwGenos <- newARKWGenos
srkwGenos <- newSRKWGenos
tkwGenos <- newTKWGenos


###############################################################
#  Organize ARKW data
###############################################################

# arkw first
arkw_delLoci <- arkwGenos[derGenoInfo[,5] == "HIGH",]   # extract delterious loci
arkw_modLoci <- arkwGenos[derGenoInfo[,5] == "MOD",]
arkw_neutLoci <- arkwGenos[derGenoInfo[,5] == "NEUT",]  # extract the neutral loci

# get missingness for each type of locus
arkw_delMiss <- rowSums(is.na(arkw_delLoci))
arkw_neutMiss <- rowSums(is.na(arkw_neutLoci))
arkw_modMiss <- rowSums(is.na(arkw_modLoci))


# get allele frequencies for each type of locus
arkw_delFreqs  <- rowSums(arkw_delLoci,na.rm=TRUE)/(ncol(arkw_delLoci)-arkw_delMiss)
arkw_neutFreqs  <- rowSums(arkw_neutLoci,na.rm=TRUE)/(ncol(arkw_neutLoci)-arkw_neutMiss)
arkw_modFreqs  <- rowSums(arkw_modLoci,na.rm=TRUE)/(ncol(arkw_modLoci)-arkw_modMiss)



###############################################################
#  Organize SRKW data
###############################################################

# arkw first
srkw_delLoci  <- srkwGenos[derGenoInfo[,5] == "HIGH",]   # extract delterious loci
srkw_modLoci  <- srkwGenos[derGenoInfo[,5] == "MOD",]
srkw_neutLoci <- srkwGenos[derGenoInfo[,5] == "NEUT",]  # extract the neutral loci

# get missingness for each type of locus
srkw_delMiss  <- rowSums(is.na(srkw_delLoci))
srkw_neutMiss <- rowSums(is.na(srkw_neutLoci))
srkw_modMiss  <- rowSums(is.na(srkw_modLoci))

# get allele frequencies for each type of locus
srkw_delFreqs  <- rowSums(srkw_delLoci,na.rm=TRUE)/(ncol(srkw_delLoci)-srkw_delMiss)
srkw_neutFreqs  <- rowSums(srkw_neutLoci,na.rm=TRUE)/(ncol(srkw_neutLoci)-srkw_neutMiss)
srkw_modFreqs  <- rowSums(srkw_modLoci,na.rm=TRUE)/(ncol(srkw_modLoci)-srkw_modMiss)

###############################################################
#  Organize TKW data
###############################################################

# arkw first
tkw_delLoci  <- tkwGenos[derGenoInfo[,5] == "HIGH",]   # extract deleterious loci
tkw_modLoci  <- tkwGenos[derGenoInfo[,5] == "MOD",]
tkw_neutLoci <- tkwGenos[derGenoInfo[,5] == "NEUT",]  # extract the neutral loci

# get missingness for each type of locus
tkw_delMiss  <- rowSums(is.na(tkw_delLoci))
tkw_neutMiss <- rowSums(is.na(tkw_neutLoci))
tkw_modMiss  <- rowSums(is.na(tkw_modLoci))

# get allele frequencies for each type of locus
tkw_delFreqs  <- rowSums(tkw_delLoci,na.rm=TRUE)/(ncol(tkw_delLoci)-tkw_delMiss)
tkw_neutFreqs  <- rowSums(tkw_neutLoci,na.rm=TRUE)/(ncol(tkw_neutLoci)-tkw_neutMiss)
tkw_modFreqs  <- rowSums(tkw_modLoci,na.rm=TRUE)/(ncol(tkw_modLoci)-tkw_modMiss)



#############################################################
# block jackknifing over loci
#############################################################
boots <- 100
modStarts <- seq(1,nrow(srkw_modLoci),round(nrow(srkw_modLoci)/boots))
modEnds <- modStarts + round(nrow(srkw_modLoci)/boots) - 1
modEnds [length(modEnds)]<- nrow(srkw_modLoci)

neutStarts <- seq(1,nrow(srkw_neutLoci),round(nrow(srkw_neutLoci)/boots))
neutEnds <- neutStarts + round(nrow(srkw_neutLoci)/boots) - 1
neutEnds [length(neutEnds)]<- nrow(srkw_neutLoci)


delStarts <- seq(1,nrow(srkw_delLoci),floor(nrow(srkw_delLoci)/boots))
delEnds <- delStarts + floor(nrow(srkw_delLoci)/boots) - 1
delEnds [length(delEnds)] <- nrow(srkw_delLoci)



# calculate L(X not Y) in the notation of Do et al.
rxyMat_raw <- matrix(NA,ncol=3,nrow=3)
# LOF mutations Rxy
rxyMat_raw[1,1] <- sum(srkw_delFreqs * (1-arkw_delFreqs),na.rm=TRUE)
rxyMat_raw[1,2] <- sum(srkw_delFreqs * (1-tkw_delFreqs),na.rm=TRUE)
rxyMat_raw[1,3] <- sum(arkw_delFreqs * (1-tkw_delFreqs),na.rm=TRUE)

# missense mutations Rxy
rxyMat_raw[2,1] <- sum(srkw_modFreqs * (1-arkw_modFreqs),na.rm=TRUE)
rxyMat_raw[2,2] <- sum(srkw_modFreqs * (1-tkw_modFreqs),na.rm=TRUE)
rxyMat_raw[2,3] <- sum(arkw_modFreqs * (1-tkw_modFreqs),na.rm=TRUE)

# neutral mutations Rxy
rxyMat_raw[3,1] <- sum(srkw_neutFreqs * (1-arkw_neutFreqs),na.rm=TRUE)
rxyMat_raw[3,2] <- sum(srkw_neutFreqs * (1-tkw_neutFreqs),na.rm=TRUE)
rxyMat_raw[3,3] <- sum(arkw_neutFreqs * (1-tkw_neutFreqs),na.rm=TRUE)



# calculate L(Y not X) in the notation of Do et al.
rxyMat_raw2 <- matrix(NA,ncol=3,nrow=3)
# LOF mutations Rxy
rxyMat_raw2[1,1] <- sum(arkw_delFreqs * (1-srkw_delFreqs),na.rm=TRUE)
rxyMat_raw2[1,2] <- sum(tkw_delFreqs * (1-srkw_delFreqs),na.rm=TRUE)
rxyMat_raw2[1,3] <- sum(tkw_delFreqs * (1-arkw_delFreqs),na.rm=TRUE)

# missense mutations Rxy
rxyMat_raw2[2,1] <- sum(arkw_modFreqs * (1-srkw_modFreqs),na.rm=TRUE)
rxyMat_raw2[2,2] <- sum(tkw_modFreqs * (1-srkw_modFreqs),na.rm=TRUE)
rxyMat_raw2[2,3] <- sum(tkw_modFreqs * (1-arkw_modFreqs),na.rm=TRUE)

# neutral mutations Rxy
rxyMat_raw2[3,1] <- sum(arkw_neutFreqs * (1-srkw_neutFreqs),na.rm=TRUE)
rxyMat_raw2[3,2] <- sum(tkw_neutFreqs * (1-srkw_neutFreqs),na.rm=TRUE)
rxyMat_raw2[3,3] <- sum(tkw_neutFreqs * (1-arkw_neutFreqs),na.rm=TRUE)

rxyFinalMat <- rxyMat_raw/rxyMat_raw2   # this is 



# Take the jackknife sample outcomes
boots <- 100   # number of replicates
# L(X not Y) in the notation of Do et al.
RxyMat_del1 <- matrix(NA,nrow=boots,ncol=3)
RxyMat_mod1 <- matrix(NA,nrow=boots,ncol=3)
RxyMat_neut1 <- matrix(NA,nrow=boots,ncol=3)

# L(Y not X) in the notation of Do et al.
RxyMat_del2 <- matrix(NA,nrow=boots,ncol=3)
RxyMat_mod2 <- matrix(NA,nrow=boots,ncol=3)
RxyMat_neut2 <- matrix(NA,nrow=boots,ncol=3)

for(i in 1:nrow(RxyMat_del1)){
  delRows <- (1:nrow(srkw_delLoci))[-(delStarts[i]:delEnds[i])]
  srkwDelBoot <- srkw_delFreqs[delRows]
  arkwDelBoot <- arkw_delFreqs[delRows]
  tkwDelBoot <- tkw_delFreqs[delRows]
  
  modRows <- (1:nrow(srkw_modLoci))[-(modStarts[i]:modEnds[i])]
  srkwModBoot <- srkw_modFreqs[modRows]
  arkwModBoot <- arkw_modFreqs[modRows]
  tkwModBoot <- tkw_modFreqs[modRows]
  
  neutRows <- (1:nrow(srkw_neutLoci))[-(neutStarts[i]:neutEnds[i])]
  srkwNeutBoot <- srkw_neutFreqs[neutRows]
  arkwNeutBoot <- arkw_neutFreqs[neutRows]
  tkwNeutBoot <- tkw_neutFreqs[neutRows]
  
  RxyMat_neut1[i,] <-c(                                    # this is equivalent to L_XnotY and L_YnotX in Do et al. (2015, Nat Gen)
    sum(srkwNeutBoot * (1-arkwNeutBoot),na.rm=TRUE),
    sum(srkwNeutBoot * (1-tkwNeutBoot),na.rm=TRUE),
    sum(arkwNeutBoot * (1-tkwNeutBoot),na.rm=TRUE))
  RxyMat_del1[i,] <-c(
    sum(srkwDelBoot * (1-arkwDelBoot),na.rm=TRUE),
    sum(srkwDelBoot * (1-tkwDelBoot),na.rm=TRUE),
    sum(arkwDelBoot * (1-tkwDelBoot),na.rm=TRUE))
  RxyMat_mod1[i,] <-c(
    sum(srkwModBoot * (1-arkwModBoot),na.rm=TRUE),
    sum(srkwModBoot * (1-tkwModBoot),na.rm=TRUE),
    sum(arkwModBoot * (1-tkwModBoot),na.rm=TRUE))
  
  
  RxyMat_neut2[i,] <-c(
    sum(arkwNeutBoot * (1-srkwNeutBoot),na.rm=TRUE),
    sum(tkwNeutBoot * (1-srkwNeutBoot),na.rm=TRUE),
    sum(tkwNeutBoot * (1-arkwNeutBoot),na.rm=TRUE))
  RxyMat_del2[i,] <-c(
    sum(arkwDelBoot * (1-srkwDelBoot),na.rm=TRUE),
    sum(tkwDelBoot * (1-srkwDelBoot),na.rm=TRUE),
    sum(tkwDelBoot * (1-arkwDelBoot),na.rm=TRUE))
  RxyMat_mod2[i,] <-c(
    sum(arkwModBoot * (1-srkwModBoot),na.rm=TRUE),
    sum(tkwModBoot * (1-srkwModBoot),na.rm=TRUE),
    sum(tkwModBoot * (1-arkwModBoot),na.rm=TRUE))
  
  print(i)
}

RxyMat_mod_final <- RxyMat_mod1/RxyMat_mod2                   # this is the ratio L_XnotY/L_YnotX in Do et al (2015 Nat Gen), i.e., R_X/Y
RxyMat_del_final <- RxyMat_del1/RxyMat_del2
RxyMat_neut_final <- RxyMat_neut1/RxyMat_neut2
###############################################
# plot results
###############################################
lofRxyEsts <- rxyFinalMat[1,]
lofRxyCIs <- matrix(NA,nrow=3,ncol=2)
for(i in 1:3){
  lofRxyCIs[i,] <- quantile(RxyMat_del_final[,i],probs=c(0.025,0.975))
}

missRxyEsts <- rxyFinalMat[2,]
missRxyCIs <- matrix(NA,nrow=3,ncol=2)
for(i in 1:3){
  missRxyCIs[i,] <- quantile(RxyMat_mod_final[,i],probs=c(0.025,0.975))
}

par(mfrow=c(1,1),xpd=TRUE,mar=c(6,6,2,2))

#LOF
plot(c(0.5,3.5),c(0.90,1.2),type="n",axes=FALSE,xlab="",ylab=expression(italic(""*R*"")["X/Y"]),cex.lab=1.3)
lines(c(0.5,3.5),c(1,1),lty="dashed")
points((1:3)-0.1,lofRxyEsts,pch=16,col="darkred",cex=2)
points((1:3)+0.1,missRxyEsts,pch=16,col="darkgray",cex=2)

for(i in 1:3){
  arrows(x0=i-0.1,x1=i-0.1,y0=lofRxyEsts[i],y1=lofRxyCIs[i,1],col="darkred",angle=90,length=0.1)
  arrows(x0=i-0.1,x1=i-0.1,y0=lofRxyEsts[i],y1=lofRxyCIs[i,2],col="darkred",angle=90,length=0.1)
  
  arrows(x0=i+0.1,x1=i+0.1,y0=missRxyEsts[i],y1=missRxyCIs[i,1],col="darkgray",angle=90,length=0.1)
  arrows(x0=i+0.1,x1=i+0.1,y0=missRxyEsts[i],y1=missRxyCIs[i,2],col="darkgray",angle=90,length=0.1)
}

axis(side=2,at=seq(0.9,1.2,0.05))
axis(side=1,at=seq(0.5,3.5,1),labels=rep("",4))
text(x=1,y=0.86,labels="SRKW / ARKW")
text(x=2,y=0.86,labels="SRKW / TKW")
text(x=3,y=0.86,labels="ARKW / TKW")

legend(x=0.5,y=1.16,legend=c("Loss of function","Missense"),col=c("darkred","darkgray"),pch=16,cex=1.3,
       xjust=FALSE,yjust=FALSE,bty="n")

write.table(missRxyEsts,file="missRxyEsts",quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(lofRxyEsts,file="lofRxyEsts",quote=FALSE,col.names=FALSE,row.names=FALSE)

write.table(missRxyCIs,file="missRxyCIs",quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(lofRxyCIs,file="lofRxyCIs",quote=FALSE,col.names=FALSE,row.names=FALSE)


############################################################################
# plot allele frequencies for each population and class of loci
############################################################################
# allele frequency table

lofMeanFreqs <- c(mean(srkw_delFreqs),sd(srkw_delFreqs),mean(arkw_delFreqs),sd(arkw_delFreqs),mean(tkw_delFreqs),sd(tkw_delFreqs))
modMeanFreqs <- c(mean(srkw_modFreqs),sd(srkw_modFreqs),mean(arkw_modFreqs),sd(arkw_modFreqs),mean(tkw_modFreqs),sd(tkw_modFreqs))
neutMeanFreqs <- c(mean(srkw_neutFreqs),sd(srkw_neutFreqs),mean(arkw_neutFreqs),sd(arkw_neutFreqs),mean(tkw_neutFreqs,na.rm=TRUE),sd(tkw_neutFreqs,na.rm=TRUE))

delFreqMat <- cbind(srkw_delFreqs,arkw_delFreqs,tkw_delFreqs)
modFreqMat <- cbind(srkw_modFreqs,arkw_modFreqs,tkw_modFreqs)
neutFreqMat <- cbind(srkw_neutFreqs,arkw_neutFreqs,tkw_neutFreqs)
neutFreqMat <- neutFreqMat[-which(rowSums(is.na(neutFreqMat)) > 0),]

delSDs <- rep(NA,3)
modSDs <- rep(NA,3)
neutSDs <- rep(NA,3)

for(i in 1:3){
  delSDs[i] <- sd(delFreqMat[,i])
  modSDs[i] <- sd(modFreqMat[,i])
  neutSDs[i] <- sd(neutFreqMat[,i])
}

lofFreqVec <- c(srkw_delFreqs,arkw_delFreqs,tkw_delFreqs)
popVec <- c(rep("SRKW",length(srkw_delFreqs)),rep("ARKW",length(arkw_delFreqs)),rep("TKW",length(tkw_delFreqs)))

pops <- c("SRKW_mean","SRKW_sd","ARKW_mean","ARKW_sd","TKW_mean","TKW_sd")
outFreqs <- cbind(pops,lofMeanFreqs,modMeanFreqs,neutMeanFreqs)
write.table(outFreqs,file="meanAlleleFreqs_LOD_missense_neutral",quote=FALSE,col.names=TRUE,row.names=FALSE)







