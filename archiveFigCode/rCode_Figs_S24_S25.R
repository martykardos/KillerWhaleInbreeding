setwd("~/Documents/orca/analyses_31May2022/RXY")
kwIDs <- read.table("kw_151.snp.final.passOnly.biallelicOnly.noIndels.mac1_indivFiltered.hetDepthFiltered_maxMiss0.25.tfam")
varEffs <- read.table("variant_effect_output_skinny")             # variant effects inferred from the Ensembl Variant Effect Predictor
library(data.table)
derGenos <- fread("kw_DerivedAlleleGenotypes_1August2022",colClasses="numeric")         # derived/ancestral genotypes (0 = ancestral, 1 = derived)
derGenoInfo <- fread("kw_DerivedAlleleGenotypeInformationMatrix_1August2022",colClasses=rep("character",4))     # chromosome and position of each locus in the genome
derGenos <- as.matrix(derGenos)
derGenoInfo <- as.matrix(derGenoInfo)
popIDs <- read.table("orca_pops_ids.txt",header=TRUE)      # link ids to populations

chroms <- as.numeric(unique(derGenoInfo[,1]))

chrLengs <- rep(NA,length(chroms))
for(i in 1:length(chroms)){
  chrLengs[i] <- max(as.numeric(derGenoInfo[derGenoInfo[,1] == as.character(chroms[i]),4]))
}
#######################################################
# organize the genotypes by population
#######################################################

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

#######################################################
# organize the variant effects
#######################################################
# remove annotation information for loci not found in the genotypes
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
# split up different classes of genotypes for each population
###############################################################

remThresh <- 0.2  # missingness frequency above which to remove a locus from any population
#-----------------
# arkw first
#-----------------
delLoci <- arkwGenos[derGenoInfo[,5] == "HIGH" | derGenoInfo[,5] == "MOD",]   # extract delterious loci
info <- derGenoInfo[derGenoInfo[,5] == "HIGH" | derGenoInfo[,5] == "MOD",]

# get missingness 
miss <- (rowSums(is.na(delLoci)))/(ncol(delLoci))

# remove loci with too many missing genotypes
delLoci <- delLoci[-which(miss >= remThresh ),]        
info <- info[-which(miss >= remThresh ),]  
arkwDelLoci <- delLoci
arkwInfo <- info
#-----------------
# srkw now
#-----------------
delLoci <- srkwGenos[derGenoInfo[,5] == "HIGH" | derGenoInfo[,5] == "MOD",]   # extract delterious loci
info <- derGenoInfo[derGenoInfo[,5] == "HIGH" | derGenoInfo[,5] == "MOD",]

# get missingness 
miss <- (rowSums(is.na(delLoci)))/(ncol(delLoci))

# remove loci with too many missing genotypes
delLoci <- delLoci[-which(miss >= remThresh ),]        
info <- info[-which(miss >= remThresh ),]  
srkwDelLoci <- delLoci
srkwInfo <- info

#-----------------
# tkw now
#-----------------
delLoci <- tkwGenos[derGenoInfo[,5] == "HIGH" | derGenoInfo[,5] == "MOD",]   # extract delterious loci
info <- derGenoInfo[derGenoInfo[,5] == "HIGH" | derGenoInfo[,5] == "MOD",]

# get missingness 
miss <- (rowSums(is.na(delLoci)))/(ncol(delLoci))

# remove loci with too many missing genotypes
delLoci <- delLoci[-which(miss >= remThresh ),]        
info <- info[-which(miss >= remThresh ),]  
tkwDelLoci <- delLoci
tkwInfo <- info
#-----------------
# offshore now
#-----------------
delLoci <- offshoreGenos[derGenoInfo[,5] == "HIGH" | derGenoInfo[,5] == "MOD",]   # extract delterious loci
info <- derGenoInfo[derGenoInfo[,5] == "HIGH" | derGenoInfo[,5] == "MOD",]

# get missingness 
miss <- (rowSums(is.na(delLoci)))/(ncol(delLoci))

# remove loci with too many missing genotypes
delLoci <- delLoci[-which(miss >= remThresh ),]        
info <- info[-which(miss >= remThresh ),]  
offshoreDelLoci <- delLoci
offshoreInfo <- info


chroms <- chroms[-6]                    # chr six is the sex chromosome
chrLengs <- chrLengs[-6]







############################################################################################
# Make a manhattan plot of the deleterious allele frequencies across the genome
# include a line showing the number of deleterious loci per 100 Kb or something like that
############################################################################################

gap <- 20              # gap between chromosomes in Mb


# southern residents
xIter <- 0 # iterate the starting positions for each chromosome on the x-axis



winSize <- 1 # 1 Mb windows
library(scales)



# srkw
par(mfrow=c(4,1),mar=c(4,6,2,1),xpd=TRUE)
xIter <- 0
plot(c(0,sum(chrLengs/1000000)+gap*(length(chroms)-1)),c(0,1),type="n",cex.lab=1.3,
     axes=FALSE,xlab="",ylab="Allele Frequency",main="Southern Resident",cex.main=1.5)
for(i in 1:length(chroms)){
  
  thisDat <- srkwDelLoci[srkwInfo[,1] == as.character(chroms[i]),]
  thisInfo <- srkwInfo[srkwInfo[,1] == as.character(chroms[i]),]
  freqs <- rowSums(thisDat,na.rm=TRUE)/rowSums(is.na(thisDat) == FALSE)
  pos <- as.numeric(thisInfo[,4])/1000000
  points(xIter + pos,freqs,col=alpha("#a6cef9",alpha=0.2),cex=0.5)
  
  #proportion fixed in 1 Mb windows
  starts <- seq(0,chrLengs[i]/1000000)
  ends <- starts + winSize
  propFixed <- rep(NA,length(starts))
  for(j in 1:length(starts)){
    if(sum(pos > starts[j] & pos <= ends[j]) > 0){
      winDat <- freqs[pos > starts[j] & pos <= ends[j]]
      propFixed[j] <- sum(winDat == 1)/length(winDat)
    }
  }
  if(sum(is.na(propFixed)) > 0){
    starts <- starts[-which(is.na(propFixed))]
    ends <- ends[-which(is.na(propFixed))]
    propFixed <- propFixed[-which(is.na(propFixed))]
  }
  lines(xIter + rowMeans(cbind(starts,ends)),propFixed,col=alpha("black",alpha=0.5),lwd=1)
  
  xIter <- xIter + gap + chrLengs[i]/1000000
}
axis(side=2,at=seq(0,1,0.5))


# arkw
xIter <- 0
plot(c(0,sum(chrLengs/1000000)+gap*(length(chroms)-1)),c(0,1),cex.lab=1.3,
     type="n",axes=FALSE,xlab="",ylab="Allele Frequency",main="Alaska Resident",cex.main=1.5)
for(i in 1:length(chroms)){
  
  thisDat <- arkwDelLoci[arkwInfo[,1] == as.character(chroms[i]),]
  thisInfo <- arkwInfo[arkwInfo[,1] == as.character(chroms[i]),]
  freqs <- rowSums(thisDat,na.rm=TRUE)/rowSums(is.na(thisDat) == FALSE)
  pos <- as.numeric(thisInfo[,4])/1000000
  points(xIter + pos,freqs,col=alpha("#1f78b4",alpha=0.2),cex=0.5)
  
  #proportion fixed in 1 Mb windows
  starts <- seq(0,chrLengs[i]/1000000)
  ends <- starts + winSize
  propFixed <- rep(NA,length(starts))
  for(j in 1:length(starts)){
    if(sum(pos > starts[j] & pos <= ends[j]) > 0){
      winDat <- freqs[pos > starts[j] & pos <= ends[j]]
      propFixed[j] <- sum(winDat == 1)/length(winDat)
    }
  }
  if(sum(is.na(propFixed)) > 0){
    starts <- starts[-which(is.na(propFixed))]
    ends <- ends[-which(is.na(propFixed))]
    propFixed <- propFixed[-which(is.na(propFixed))]
  }
  lines(xIter + rowMeans(cbind(starts,ends)),propFixed,col=alpha("black",alpha=0.5),lwd=1)
  
  xIter <- xIter + gap + chrLengs[i]/1000000
}
axis(side=2,at=seq(0,1,0.5))



# tkw
xIter <- 0
plot(c(0,sum(chrLengs/1000000)+gap*(length(chroms)-1)),c(0,1),cex.lab=1.3,
     type="n",axes=FALSE,xlab="",ylab="Allele Frequency",main="Transient",cex.main=1.5)
for(i in 1:length(chroms)){
  thisDat <- tkwDelLoci[tkwInfo[,1] == as.character(chroms[i]),]
  thisInfo <- tkwInfo[tkwInfo[,1] == as.character(chroms[i]),]
  freqs <- rowSums(thisDat,na.rm=TRUE)/rowSums(is.na(thisDat) == FALSE)
  pos <- as.numeric(thisInfo[,4])/1000000
  points(xIter + pos,freqs,col=alpha("#b2df8a",alpha=0.2),cex=0.5)
  
  #proportion fixed in 1 Mb windows
  starts <- seq(0,chrLengs[i]/1000000)
  ends <- starts + winSize
  propFixed <- rep(NA,length(starts))
  for(j in 1:length(starts)){
    if(sum(pos > starts[j] & pos <= ends[j]) > 0){
      winDat <- freqs[pos > starts[j] & pos <= ends[j]]
      propFixed[j] <- sum(winDat == 1)/length(winDat)
    }
  }
  if(sum(is.na(propFixed)) > 0){
    starts <- starts[-which(is.na(propFixed))]
    ends <- ends[-which(is.na(propFixed))]
    propFixed <- propFixed[-which(is.na(propFixed))]
  }
  lines(xIter + rowMeans(cbind(starts,ends)),propFixed,col=alpha("black",alpha=0.5),lwd=1)
  
  xIter <- xIter + gap + chrLengs[i]/1000000
}
axis(side=2,at=seq(0,1,0.5))



# offshore
xIter <- 0
par(xpd=TRUE)
plot(c(0,sum(chrLengs/1000000)+gap*(length(chroms)-1)),c(0,1),type="n",cex.lab=1.3,
     axes=FALSE,ylab="Allele Frequency",main="Offshore",cex.main=1.5,xlab="")
for(i in 1:length(chroms)){
  thisDat <- offshoreDelLoci[offshoreInfo[,1] == as.character(chroms[i]),]
  thisInfo <- offshoreInfo[offshoreInfo[,1] == as.character(chroms[i]),]
  freqs <- rowSums(thisDat,na.rm=TRUE)/rowSums(is.na(thisDat) == FALSE)
  pos <- as.numeric(thisInfo[,4])/1000000
  points(xIter + pos,freqs,col=alpha("#33a02c",alpha=0.2),cex=0.5)
  
  #proportion fixed in 1 Mb windows
  starts <- seq(0,chrLengs[i]/1000000)
  ends <- starts + winSize
  propFixed <- rep(NA,length(starts))
  for(j in 1:length(starts)){
    if(sum(pos > starts[j] & pos <= ends[j]) > 0){
      winDat <- freqs[pos > starts[j] & pos <= ends[j]]
      propFixed[j] <- sum(winDat == 1)/length(winDat)
    }
  }
  if(sum(is.na(propFixed)) > 0){
    starts <- starts[-which(is.na(propFixed))]
    ends <- ends[-which(is.na(propFixed))]
    propFixed <- propFixed[-which(is.na(propFixed))]
  }
  lines(xIter + rowMeans(cbind(starts,ends)),propFixed,col=alpha("black",alpha=0.5),lwd=1)
  text(x=xIter + mean(ends),y=-0.15,labels=chroms[i],cex=1.2)
  xIter <- xIter + gap + chrLengs[i]/1000000
  
}
axis(side=2,at=seq(0,1,0.5))
text(x=(sum(chrLengs/1000000)+gap*(length(chroms)-1))*0.5,labels="Chromsome",cex=2,y=-0.35)


####################################################################
# plot allele frequency differences for each pair of populations
####################################################################

par(mar=c(5,6,2,4))
srkwFreqs <-rowSums(srkwDelLoci,na.rm=TRUE)/rowSums(is.na(srkwDelLoci) == FALSE)
arkwFreqs <-rowSums(arkwDelLoci,na.rm=TRUE)/rowSums(is.na(arkwDelLoci) == FALSE)
tkwFreqs <-rowSums(tkwDelLoci,na.rm=TRUE)/rowSums(is.na(tkwDelLoci) == FALSE)
offshoreFreqs <-rowSums(offshoreDelLoci,na.rm=TRUE)/rowSums(is.na(offshoreDelLoci) == FALSE)

srkwFreqMat <- cbind(srkwInfo,srkwFreqs)
arkwFreqMat <- cbind(arkwInfo,arkwFreqs)
tkwFreqMat <- cbind(tkwInfo,tkwFreqs)
offshoreFreqMat <- cbind(offshoreInfo,offshoreFreqs)


par(mfrow=c(2,3))
#-------------------
# srkw vs. arkw
#-------------------

thisSrkw <- srkwFreqMat[which(srkwFreqMat[,2] %in% arkwFreqMat[,2]),]
thisArkw <- arkwFreqMat[which(arkwFreqMat[,2] %in% srkwFreqMat[,2]),]
thisSrkw <- thisSrkw[match(thisSrkw[,2],thisArkw[,2]),]
plot(as.numeric(thisArkw[,6]),as.numeric(thisSrkw[,6]),col=alpha("blue",alpha=0.02),xlab="",ylab="",pch=16,cex=1.5)
    text(x=0.5,y=-0.18,expression(paste(italic(""*p*""),"(ARKW)")),col="#1f78b4",font=2,cex=2)
    text(x=-0.21,y=0.5,expression(paste(italic(""*p*""),"(SRKW)")),col="#a6cef9",font=2,cex=2,srt=90)
    
cor(as.numeric(thisArkw[,6]),as.numeric(thisSrkw[,6]))^2
text(x=0.2,y=0.8,labels = expression(paste(italic(""*r*"")^2," = ","0.85",sep="")),cex=2)    


#-------------------
# srkw vs. tkw
#-------------------
    
    thisSrkw <- srkwFreqMat[which(srkwFreqMat[,2] %in% tkwFreqMat[,2]),]
    thisTkw <- tkwFreqMat[which(tkwFreqMat[,2] %in% srkwFreqMat[,2]),]
    thisSrkw <- thisSrkw[match(thisSrkw[,2],thisTkw[,2]),]
    plot(as.numeric(thisTkw[,6]),as.numeric(thisSrkw[,6]),col=alpha("blue",alpha=0.02),xlab="",ylab="",pch=16,cex=1.5)
    text(x=0.5,y=-0.18,expression(paste(italic(""*p*""),"(TKW)")),col="#b2df8a",font=2,cex=2)
    text(x=-0.21,y=0.5,expression(paste(italic(""*p*""),"(SRKW)")),col="#a6cef9",font=2,cex=2,srt=90)
    
    cor(as.numeric(thisSrkw[,6]),as.numeric(thisTkw[,6]))^2
    text(x=0.2,y=0.8,labels = expression(paste(italic(""*r*"")^2," = ","0.34",sep="")),cex=2)    
    
    
    
    
#-------------------
# srkw vs. offshore
#-------------------
    
    thisSrkw <- srkwFreqMat[which(srkwFreqMat[,2] %in% offshoreFreqMat[,2]),]
    thisOff <- offshoreFreqMat[which(offshoreFreqMat[,2] %in% srkwFreqMat[,2]),]
    thisSrkw <- thisSrkw[match(thisSrkw[,2],thisOff[,2]),]
    plot(as.numeric(thisOff[,6]),as.numeric(thisSrkw[,6]),col=alpha("blue",alpha=0.02),xlab="",ylab="",pch=16,cex=1.5)
    text(x=0.5,y=-0.18,expression(paste(italic(""*p*""),"(Offshore)")),col="#33a02c",font=2,cex=2)
    text(x=-0.21,y=0.5,expression(paste(italic(""*p*""),"(SRKW)")),col="#a6cef9",font=2,cex=2,srt=90)
    
    cor(as.numeric(thisSrkw[,6]),as.numeric(thisOff[,6]))^2
    text(x=0.2,y=0.8,labels = expression(paste(italic(""*r*"")^2," = ","0.35",sep="")),cex=2)    
    
    
    
#-------------------
# arkw vs. tkw
#-------------------
    
    thisArkw <- arkwFreqMat[which(arkwFreqMat[,2] %in% tkwFreqMat[,2]),]
    thisTkw <- tkwFreqMat[which(tkwFreqMat[,2] %in% arkwFreqMat[,2]),]
    thisArkw <- thisArkw[match(thisArkw[,2],thisTkw[,2]),]
    plot(as.numeric(thisTkw[,6]),as.numeric(thisArkw[,6]),col=alpha("blue",alpha=0.02),xlab="",ylab="",pch=16,cex=1.5)
    text(x=0.5,y=-0.18,expression(paste(italic(""*p*""),"(TKW)")),col="#b2df8a",font=2,cex=2)
    text(x=-0.21,y=0.5,expression(paste(italic(""*p*""),"(ARKW)")),col="#1f78b4",font=2,cex=2,srt=90)
    
    cor(as.numeric(thisArkw[,6]),as.numeric(thisTkw[,6]))^2
    text(x=0.2,y=0.8,labels = expression(paste(italic(""*r*"")^2," = ","0.38",sep="")),cex=2)    
    
    
#-------------------
# arkw vs. offshore
#-------------------
    
    thisArkw <- arkwFreqMat[which(arkwFreqMat[,2] %in% offshoreFreqMat[,2]),]
    thisOff <- offshoreFreqMat[which(offshoreFreqMat[,2] %in% arkwFreqMat[,2]),]
    thisArkw <- thisArkw[match(thisArkw[,2],thisOff[,2]),]
    plot(as.numeric(thisOff[,6]),as.numeric(thisArkw[,6]),col=alpha("blue",alpha=0.02),xlab="",ylab="",pch=16,cex=1.5)
    text(x=0.5,y=-0.18,expression(paste(italic(""*p*""),"(Offshore)")),col="#33a02c",font=2,cex=2)
    text(x=-0.21,y=0.5,expression(paste(italic(""*p*""),"(ARKW)")),col="#1f78b4",font=2,cex=2,srt=90)
    
    cor(as.numeric(thisArkw[,6]),as.numeric(thisOff[,6]))^2
    text(x=0.2,y=0.8,labels = expression(paste(italic(""*r*"")^2," = ","0.39",sep="")),cex=2)    
    
#-------------------
# tkw vs. offshore
#-------------------
    
    thisTkw <- tkwFreqMat[which(tkwFreqMat[,2] %in% offshoreFreqMat[,2]),]
    thisOff <- offshoreFreqMat[which(offshoreFreqMat[,2] %in% tkwFreqMat[,2]),]
    thisTkw <- thisTkw[match(thisTkw[,2],thisOff[,2]),]
    plot(as.numeric(thisOff[,6]),as.numeric(thisTkw[,6]),col=alpha("blue",alpha=0.02),xlab="",ylab="",pch=16,cex=1.5)
    text(x=0.5,y=-0.18,expression(paste(italic(""*p*""),"(Offshore)")),col="#33a02c",font=2,cex=2)
    text(x=-0.21,y=0.5,expression(paste(italic(""*p*""),"(TKW)")),col="#b2df8a",font=2,cex=2,srt=90)

    cor(as.numeric(thisTkw[,6]),as.numeric(thisOff[,6]))^2
    text(x=0.2,y=0.8,labels = expression(paste(italic(""*r*"")^2," = ","0.44",sep="")),cex=2)    
    