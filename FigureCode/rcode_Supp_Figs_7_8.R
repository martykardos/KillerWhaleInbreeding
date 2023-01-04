# analyze the genomic distribution of ROH in killer whales

setwd("~/Documents/orca/analyses_31May2022/ROH3")
ibdTracts <- read.table("allROH_18August2022",header=TRUE)

numChroms <- c(1:5,7:22)

################################################################################
## get the roh data
################################################################################

Froh <- read.table("Froh_7August2022",header=TRUE)
froh <- Froh

popKey <- read.table("popKey.txt",header=TRUE)
minROHSize <- 1000000
winMat <- NULL
allInds <- unique(ibdTracts[,1])
chroms <- unique(ibdTracts[,2])
winSize <- 1000000
chrLengs <- rep(NA,length(chroms))
indWinInbs <- NULL
for(i in 1:length(chroms)){
  thisDat <- ibdTracts[which(ibdTracts[,2] == chroms[i] & (ibdTracts[,4] - ibdTracts[,3]) >= minROHSize),]
  chrLeng <- max(thisDat[,4])
  chrLengs[i] <- chrLeng
  starts <- seq(0,chrLeng,winSize)
  ends <- starts + winSize
  ends[length(ends)] <- chrLeng
  rohCounts <- rep(NA,length(starts))
  winInb <- rep(NA,length(starts))
  chromIndWinInbs <- matrix(0,nrow=length(allInds),ncol=length(starts))

  for(j in 1:length(starts)){
    rohCounts[j] <- sum(starts[j] < thisDat[,4] & ends[j] > thisDat[,3])

    if(rohCounts[j] > 0){
      theseROH <- thisDat[starts[j] < thisDat[,4] & ends[j] > thisDat[,3],c(1,3,4)]
      if(sum(theseROH[,2] < starts[j]) > 0) theseROH[which(theseROH[,2] < starts[j]),2] <- starts[j]
      if(sum(theseROH[,3] > ends[j]) > 0) theseROH[which(theseROH[,3] > ends[j]),3] <- ends[j]

      thisInb <- sum((theseROH[,3] - theseROH[,2])/(ends[j] - starts[j]))/length(unique(ibdTracts[,1]))
      if(thisInb <= 1) winInb[j] <- thisInb

      indInbs <- (theseROH[,3] - theseROH[,2])/(ends[j] - starts[j])
      uniqInds <- unique(theseROH[,1])
      for(k in 1:length(uniqInds)){
        chromIndWinInbs[match(uniqInds[k],allInds),j] <- sum(theseROH[match(uniqInds[k],theseROH[,1]),3] - theseROH[match(uniqInds[k],theseROH[,1]),2])/((ends[j] - starts[j])*length(match(uniqInds[k],theseROH[,1])))
      }
    }
  }
  indWinInbs <- cbind(indWinInbs,chromIndWinInbs)
  out <- cbind(rep(chroms[i], length(starts)),starts,ends,rohCounts,winInb)
  winMat <- rbind(winMat,out)
  print(i)
}

winMat[winMat[,4] == "0",5] <- 0


################################################################################
################################################################################
################################################################################
# Plot ROH density for each population along with recombination rate at the top
# minimum ROH length = 1 Mb
################################################################################
################################################################################
################################################################################
# analyze the genomic distribution of ROH in killer whales

ibdTracts <- read.table("allROH_18August2022",header=TRUE)

numChroms <- c(1:5,7:22)

################################################################################
## get the roh data
################################################################################

Froh <- read.table("Froh_7August2022",header=TRUE)


# get the Froh data
popKey <- read.table("popKey.txt",header=TRUE)


#------------------------------------------------------------------------
# get ROH density for all populations and for each population separately
#------------------------------------------------------------------------
minMb <- 1
popKey <- read.table("popKey.txt",header=TRUE)
ibdTracts2 <- ibdTracts[ibdTracts$lengths > minMb,]
srkws <- popKey[which(popKey[,1] == "SRKW"),2]
arkws <- popKey[which(popKey[,1] == "ARKW"),2]
tkws <- popKey[which(popKey[,1] == "TKW"),2]
offshores <- popKey[which(popKey[,1] == "offshore"),2]
nrkws <- popKey[which(popKey[,1] == "NRKW"),2]
minROHSize <- 1000000
winMat <- NULL
allInds <- unique(ibdTracts[,1])
chroms <- unique(ibdTracts[,2])
winSize <- 100000

################ all pops combined
indWinInbs <- NULL
for(i in 1:length(chroms)){
  thisDat <- ibdTracts[which(ibdTracts[,2] == chroms[i] & ibdTracts$length >= minROHSize),]
  chrLeng <-  chrLengs[i]
  starts <- seq(0,chrLeng,winSize)
  ends <- starts + winSize
  ends[length(ends)] <- chrLeng
  rohCounts <- rep(NA,length(starts))
  winInb <- rep(NA,length(starts))
  chromIndWinInbs <- matrix(0,nrow=length(allInds),ncol=length(starts))
  for(j in 1:length(starts)){
    rohCounts[j] <- sum(starts[j] < thisDat[,4] & ends[j] > thisDat[,3])
    if(rohCounts[j] > 0){
      theseROH <- thisDat[starts[j] < thisDat[,4] & ends[j] > thisDat[,3],c(1,3,4)]
      if(sum(theseROH[,2] < starts[j]) > 0) theseROH[which(theseROH[,2] < starts[j]),2] <- starts[j]
      if(sum(theseROH[,3] > ends[j]) > 0) theseROH[which(theseROH[,3] > ends[j]),3] <- ends[j]
      thisInb <- sum((theseROH[,3] - theseROH[,2])/(ends[j] - starts[j]))/length(unique(ibdTracts[,1]))
      if(thisInb <= 1) winInb[j] <- thisInb
      indInbs <- (theseROH[,3] - theseROH[,2])/(ends[j] - starts[j])
      uniqInds <- unique(theseROH[,1])
      for(k in 1:length(uniqInds)){
        chromIndWinInbs[match(uniqInds[k],allInds),j] <- sum(theseROH[match(uniqInds[k],theseROH[,1]),3] - theseROH[match(uniqInds[k],theseROH[,1]),2])/((ends[j] - starts[j])*length(match(uniqInds[k],theseROH[,1])))
      }
    }
  }
  indWinInbs <- cbind(indWinInbs,chromIndWinInbs)
  out <- cbind(rep(chroms[i], length(starts)),starts,ends,rohCounts,winInb)
  winMat <- rbind(winMat,out)
  print(i)
}
winMat[winMat[,4] == "0",5] <- 0

################ SRKWs
winMat_srkw <- NULL
indWinInbs_srkw <- NULL
for(i in 1:length(chroms)){
  thisDat <- ibdTracts[which(ibdTracts[,2] == chroms[i] & (ibdTracts[,4] - ibdTracts[,3]) >= minROHSize & ibdTracts[,1] %in% srkws),]
  chrLeng <- chrLengs[i]

  starts <- seq(0,chrLeng,winSize)
  ends <- starts + winSize
  ends[length(ends)] <- chrLeng
  rohCounts <- rep(NA,length(starts))
  winInb <- rep(NA,length(starts))
  chromIndWinInbs <- matrix(0,nrow=length(allInds),ncol=length(starts))
  for(j in 1:length(starts)){
    rohCounts[j] <- sum(starts[j] < thisDat[,4] & ends[j] > thisDat[,3])
    if(rohCounts[j] > 0){
      theseROH <- thisDat[starts[j] < thisDat[,4] & ends[j] > thisDat[,3],c(1,3,4)]
      if(sum(theseROH[,2] < starts[j]) > 0) theseROH[which(theseROH[,2] < starts[j]),2] <- starts[j]
      if(sum(theseROH[,3] > ends[j]) > 0) theseROH[which(theseROH[,3] > ends[j]),3] <- ends[j]
      thisInb <- sum((theseROH[,3] - theseROH[,2])/(ends[j] - starts[j]))/length(unique(ibdTracts[,1]))
      if(thisInb <= 1) winInb[j] <- thisInb
      indInbs <- (theseROH[,3] - theseROH[,2])/(ends[j] - starts[j])
      uniqInds <- unique(theseROH[,1])
      for(k in 1:length(uniqInds)){
        chromIndWinInbs[match(uniqInds[k],allInds),j] <- sum(theseROH[match(uniqInds[k],theseROH[,1]),3] - theseROH[match(uniqInds[k],theseROH[,1]),2])/((ends[j] - starts[j])*length(match(uniqInds[k],theseROH[,1])))
      }
    }
  }
  indWinInbs_srkw <- cbind(indWinInbs_srkw,chromIndWinInbs)
  out <- cbind(rep(chroms[i], length(starts)),starts,ends,rohCounts,winInb)
  winMat_srkw <- rbind(winMat_srkw,out)
  print(i)
}
winMat_srkw[winMat_srkw[,4] == "0",5] <- 0

################ ARKWs
winMat_arkw <- NULL
indWinInbs_arkw <- NULL
for(i in 1:length(chroms)){
  thisDat <- ibdTracts[which(ibdTracts[,2] == chroms[i] & (ibdTracts[,4] - ibdTracts[,3]) >= minROHSize & ibdTracts[,1] %in% arkws),]
  chrLeng <- chrLengs[i]
  starts <- seq(0,chrLeng,winSize)
  ends <- starts + winSize
  ends[length(ends)] <- chrLeng
  rohCounts <- rep(NA,length(starts))
  winInb <- rep(NA,length(starts))
  chromIndWinInbs <- matrix(0,nrow=length(allInds),ncol=length(starts))
  for(j in 1:length(starts)){
    rohCounts[j] <- sum(starts[j] < thisDat[,4] & ends[j] > thisDat[,3])
    if(rohCounts[j] > 0){
      theseROH <- thisDat[starts[j] < thisDat[,4] & ends[j] > thisDat[,3],c(1,3,4)]
      if(sum(theseROH[,2] < starts[j]) > 0) theseROH[which(theseROH[,2] < starts[j]),2] <- starts[j]
      if(sum(theseROH[,3] > ends[j]) > 0) theseROH[which(theseROH[,3] > ends[j]),3] <- ends[j]
      thisInb <- sum((theseROH[,3] - theseROH[,2])/(ends[j] - starts[j]))/length(unique(ibdTracts[,1]))
      if(thisInb <= 1) winInb[j] <- thisInb
      indInbs <- (theseROH[,3] - theseROH[,2])/(ends[j] - starts[j])
      uniqInds <- unique(theseROH[,1])
      for(k in 1:length(uniqInds)){
        chromIndWinInbs[match(uniqInds[k],allInds),j] <- sum(theseROH[match(uniqInds[k],theseROH[,1]),3] - theseROH[match(uniqInds[k],theseROH[,1]),2])/((ends[j] - starts[j])*length(match(uniqInds[k],theseROH[,1])))
      }
    }
  }
  indWinInbs_arkw <- cbind(indWinInbs_arkw,chromIndWinInbs)
  out <- cbind(rep(chroms[i], length(starts)),starts,ends,rohCounts,winInb)
  winMat_arkw <- rbind(winMat_arkw,out)
  print(i)
}
winMat_arkw[winMat_arkw[,4] == "0",5] <- 0

################ offshores
winMat_offshore <- NULL
indWinInbs_offshore <- NULL
for(i in 1:length(chroms)){
  thisDat <- ibdTracts[which(ibdTracts[,2] == chroms[i] & (ibdTracts[,4] - ibdTracts[,3]) >= minROHSize & ibdTracts[,1] %in% offshores),]
  chrLeng <- chrLengs[i]

  starts <- seq(0,chrLeng,winSize)
  ends <- starts + winSize
  ends[length(ends)] <- chrLeng
  rohCounts <- rep(NA,length(starts))
  winInb <- rep(NA,length(starts))
  chromIndWinInbs <- matrix(0,nrow=length(allInds),ncol=length(starts))
  for(j in 1:length(starts)){
    rohCounts[j] <- sum(starts[j] < thisDat[,4] & ends[j] > thisDat[,3])
    if(rohCounts[j] > 0){
      theseROH <- thisDat[starts[j] < thisDat[,4] & ends[j] > thisDat[,3],c(1,3,4)]
      if(sum(theseROH[,2] < starts[j]) > 0) theseROH[which(theseROH[,2] < starts[j]),2] <- starts[j]
      if(sum(theseROH[,3] > ends[j]) > 0) theseROH[which(theseROH[,3] > ends[j]),3] <- ends[j]
      thisInb <- sum((theseROH[,3] - theseROH[,2])/(ends[j] - starts[j]))/length(unique(ibdTracts[,1]))
      if(thisInb <= 1) winInb[j] <- thisInb
      indInbs <- (theseROH[,3] - theseROH[,2])/(ends[j] - starts[j])
      uniqInds <- unique(theseROH[,1])
      for(k in 1:length(uniqInds)){
        chromIndWinInbs[match(uniqInds[k],allInds),j] <- sum(theseROH[match(uniqInds[k],theseROH[,1]),3] - theseROH[match(uniqInds[k],theseROH[,1]),2])/((ends[j] - starts[j])*length(match(uniqInds[k],theseROH[,1])))
      }
    }
  }
  indWinInbs_offshore <- cbind(indWinInbs_offshore,chromIndWinInbs)
  out <- cbind(rep(chroms[i], length(starts)),starts,ends,rohCounts,winInb)
  winMat_offshore <- rbind(winMat_offshore,out)
  print(i)
}
winMat_offshore[winMat_offshore[,4] == "0",5] <- 0


################ TKWs
winMat_tkw <- NULL
indWinInbs_tkw <- NULL
for(i in 1:length(chroms)){
  thisDat <- ibdTracts[which(ibdTracts[,2] == chroms[i] & (ibdTracts[,4] - ibdTracts[,3]) >= minROHSize & ibdTracts[,1] %in% tkws),]
  chrLeng <- chrLengs[i]

  starts <- seq(0,chrLeng,winSize)
  ends <- starts + winSize
  ends[length(ends)] <- chrLeng
  rohCounts <- rep(NA,length(starts))
  winInb <- rep(NA,length(starts))
  chromIndWinInbs <- matrix(0,nrow=length(allInds),ncol=length(starts))
  for(j in 1:length(starts)){
    rohCounts[j] <- sum(starts[j] < thisDat[,4] & ends[j] > thisDat[,3])
    if(rohCounts[j] > 0){
      theseROH <- thisDat[starts[j] < thisDat[,4] & ends[j] > thisDat[,3],c(1,3,4)]
      if(sum(theseROH[,2] < starts[j]) > 0) theseROH[which(theseROH[,2] < starts[j]),2] <- starts[j]
      if(sum(theseROH[,3] > ends[j]) > 0) theseROH[which(theseROH[,3] > ends[j]),3] <- ends[j]

      thisInb <- sum((theseROH[,3] - theseROH[,2])/(ends[j] - starts[j]))/length(unique(ibdTracts[,1]))
      if(thisInb <= 1) winInb[j] <- thisInb

      indInbs <- (theseROH[,3] - theseROH[,2])/(ends[j] - starts[j])
      uniqInds <- unique(theseROH[,1])
      for(k in 1:length(uniqInds)){
        chromIndWinInbs[match(uniqInds[k],allInds),j] <- sum(theseROH[match(uniqInds[k],theseROH[,1]),3] - theseROH[match(uniqInds[k],theseROH[,1]),2])/((ends[j] - starts[j])*length(match(uniqInds[k],theseROH[,1])))
      }
    }
  }
  indWinInbs_tkw <- cbind(indWinInbs_tkw,chromIndWinInbs)
  out <- cbind(rep(chroms[i], length(starts)),starts,ends,rohCounts,winInb)
  winMat_tkw <- rbind(winMat_tkw,out)
  print(i)
}
winMat_tkw[winMat_tkw[,4] == "0",5] <- 0


#############################################################################################################
# plot ROH and ROH density for each population across the whole genome with density for each population
#############################################################################################################

minMb <- 1000000


ibdTracts2 <- ibdTracts[ibdTracts[,5] > minMb,]
srkws <- popKey[which(popKey[,1] == "SRKW"),2]
srkwFroh <- froh[match(srkws,froh[,1]),3]
srkws <- srkws[order(srkwFroh)]

arkws <- popKey[which(popKey[,1] == "ARKW"),2]
arkwFroh <- froh[match(arkws,froh[,1]),3]
arkws <- arkws[order(arkwFroh)]

offshores <- popKey[which(popKey[,1] == "offshore"),2]
offshoreFroh <- froh[match(offshores,froh[,1]),3]
offshores <- offshores[order(offshoreFroh)]

tkws <- popKey[which(popKey[,1] == "TKW"),2]
tkwFroh <- froh[match(tkws,froh[,1]),3]
tkws <- tkws[order(tkwFroh)]

nrkws <- popKey[which(popKey[,1] == "NRKW"),2]
nrkwFroh <- froh[match(nrkws,froh[,1]),3]
nrkws <- nrkws[order(nrkwFroh)]




par(mar=c(5,7,2,2))
#####################
# initialize the plot
#####################
gap <- 10000000
plot(c(0,sum(chrLengs/1000000) + (gap*20)/1000000),c(0,149 +157),type="n",
     ylab="",xlab="Chromosome Position",axes=FALSE,cex.lab=1.5)

#######----------------------------
# SRKWs first
#######----------------------------
yIter <- 1
theseTracts <- ibdTracts2[which(ibdTracts2[,1] %in% srkws),]
plotIndivs <- srkws[-which(srkws == "J32.fetus")]
plotIndFroh <- froh[match(plotIndivs,froh[,1]),3]
plotIndivs <- plotIndivs[order(plotIndFroh)]
par(xpd=TRUE)



chroms <- unique(ibdTracts[,2])
colVec <- rep(c("black","darkgray"),50)
xIter <- 0
for(i in 1:length(chroms)){
  rect(xleft=xIter,xright=xIter + chrLengs[i]/1000000,ybottom=yIter,ytop=149 +157 ,col="gray97",border="gray97")
  text(x=mean(c(xIter,xIter + chrLengs[i]/1000000)),y=-8,labels=c(1:5,7:22)[i])
  xIter <- xIter + gap/1000000 + chrLengs[i]/1000000
}



for(i in 1:length(plotIndivs)){
  thisDat <- theseTracts[theseTracts == plotIndivs[i],]
  xIter <- 0
  for(j in 1:length(chroms)){
    if(chroms[j] %in% thisDat[,2]){
      chrRows <- which(thisDat[,2] == chroms[j])
      for(k in 1:length(chrRows)){
        thisTract <- thisDat[chrRows[k],3:4]
        rect(xleft = xIter + thisTract[1]/1000000,xright=xIter + thisTract[2]/1000000,ybottom=yIter-0.45,ytop=yIter+0.45,col="#a6cef9",border=NA)
      }
    }
    xIter <- xIter + gap/1000000 + chrLengs[j]/1000000
  }
  yIter <- yIter + 1
}

# add roh density

xIter <- 0
densPlotMin <- yIter
densPlotMax <- yIter + 20
for(i in 1:length(chroms)){
  chrDat <- winMat_srkw[winMat_srkw[,1] == chroms[i],]
  lines( (as.numeric(chrDat[,3]) - 0.5*winSize)/1000000 + xIter,densPlotMin+(as.numeric(chrDat[,5])/0.6)*20 ,col="#a6cef9",lwd=1)
  xIter <- xIter + chrLengs[i]/1000000 + gap/1000000
}
axis(side=2,at=c(yIter,yIter+20),labels=c(0,0.6),pos=-10,cex.axis=0.5)
text(x=-175,y=mean(c(yIter,yIter + 20)),labels ="ROH",cex=0.75,col="#a6cef9",srt=90)
text(x=-125,y=mean(c(yIter,yIter + 20)),labels ="Density",cex=0.75,col="#a6cef9",srt=90)

yIter <- yIter + 20
rect(xleft=0,xright=sum(chrLengs/1000000) + (gap*20)/1000000,ybottom=yIter+0.5,ytop=yIter + 4,col="white",border="white")
yIter <- yIter + 5



text(x=-225,y=51,labels="Southern Residents", col="#a6cef9",cex=1.0)


#######----------------------------
# ARKWs
#######----------------------------
theseTracts <- ibdTracts2[which(ibdTracts2[,1] %in% arkws),]
plotIndivs <- arkws
plotIndFroh <- froh[match(plotIndivs,froh[,1]),5]
plotIndivs <- plotIndivs[order(plotIndFroh)]
par(xpd=TRUE)
gap <- 10000000
chroms <- unique(ibdTracts[,2])
colVec <- rep(c("blue","orange"),50)

for(i in 1:length(arkws)){
  if(plotIndivs[i] %in% theseTracts[,1]){
    thisDat <- theseTracts[theseTracts == plotIndivs[i],]
    xIter <- 0
    for(j in 1:length(chroms)){
      if(chroms[j] %in% thisDat[,2]){
        chrRows <- which(thisDat[,2] == chroms[j])
        for(k in 1:length(chrRows)){
          thisTract <- thisDat[chrRows[k],3:4]
          rect(xleft = xIter + thisTract[1]/1000000,xright=xIter + thisTract[2]/1000000,ybottom=yIter-0.45,ytop=yIter+0.45,col="#1f78b4",border=NA)
        }
      }
      xIter <- xIter + gap/1000000 + chrLengs[j]/1000000
    }
  }
  yIter <- yIter + 1
}



# add roh density
xIter <- 0
densPlotMin <- yIter
densPlotMax <- yIter + 20
yRef <- max(as.numeric(winMat[,4]))
for(i in 1:length(chroms)){
  chrDat <- winMat_arkw[winMat_arkw[,1] == chroms[i],]
  lines( (as.numeric(chrDat[,3]) - 0.5*winSize)/1000000 + xIter,densPlotMin+(as.numeric(chrDat[,5])/0.13)*20 ,col="#1f78b4",lwd=1)
  xIter <- xIter + chrLengs[i]/1000000 + gap/1000000
}
axis(side=2,at=c(yIter,yIter+20),labels=c(0,0.13),pos=-10,cex.axis=0.5)
text(x=-175,y=mean(c(yIter,yIter + 20)),labels ="ROH",cex=0.75,col="#1f78b4",srt=90)
text(x=-125,y=mean(c(yIter,yIter + 20)),labels ="Density",cex=0.75,col="#1f78b4",srt=90)

yIter <- yIter + 20
rect(xleft=0,xright=sum(chrLengs/1000000) + (gap*20)/1000000,ybottom=yIter+0.5,ytop=yIter + 4,col="white",border="white")
yIter <- yIter + 5



text(x=-225,y=137.5,labels="Alaska Residents", col="#1f78b4",cex=1.0)

#######----------------------------
# TKWs
#######----------------------------
theseTracts <- ibdTracts2[which(ibdTracts2[,1] %in% tkws),]
plotIndivs <- tkws
plotIndFroh <- froh[match(plotIndivs,froh[,1]),5]
plotIndivs <- plotIndivs[order(plotIndFroh)]
par(xpd=TRUE)
gap <- 10000000


chroms <- unique(ibdTracts[,2])
colVec <- rep(c("blue","orange"),50)



for(i in 1:length(plotIndivs)){
  if(plotIndivs[i] %in% theseTracts[,1]){
    thisDat <- theseTracts[theseTracts == plotIndivs[i],]
    xIter <- 0
    for(j in 1:length(chroms)){
      if(chroms[j] %in% thisDat[,2]){
        chrRows <- which(thisDat[,2] == chroms[j])
        for(k in 1:length(chrRows)){
          thisTract <- thisDat[chrRows[k],3:4]
          rect(xleft = xIter + thisTract[1]/1000000,xright=xIter + thisTract[2]/1000000,ybottom=yIter-0.45,ytop=yIter+0.45,col="#b2df8a",border=NA)
        }
      }
      xIter <- xIter + gap/1000000 + chrLengs[j]/1000000
    }
  }
  yIter <- yIter + 1
}



# add roh density
xIter <- 0
densPlotMin <- yIter
densPlotMax <- yIter + 20
for(i in 1:length(chroms)){
  chrDat <- winMat_tkw[winMat_tkw[,1] == chroms[i],]
  lines( (as.numeric(chrDat[,3]) - 0.5*winSize)/1000000 + xIter,densPlotMin+(as.numeric(chrDat[,5])/0.06)*20 ,col="#b2df8a",lwd=1)
  xIter <- xIter + chrLengs[i]/1000000 + gap/1000000
}
axis(side=2,at=c(yIter,yIter+20),labels=c(0,0.06),pos=-10,cex.axis=0.5)
text(x=-175,y=mean(c(yIter,yIter + 20)),labels ="ROH",cex=0.75,col="#b2df8a",srt=90)
text(x=-125,y=mean(c(yIter,yIter + 20)),labels ="Density",cex=0.75,col="#b2df8a",srt=90)

yIter <- yIter + 20
rect(xleft=0,xright=sum(chrLengs/1000000) + (gap*20)/1000000,ybottom=yIter+0.5,ytop=yIter + 4,col="white",border="white")
yIter <- yIter + 5




text(x=-225,y=183,labels="Transients", col="#b2df8a",cex=1.0)


#######----------------------------
# offshores
#######----------------------------
theseTracts <- ibdTracts2[which(ibdTracts2[,1] %in% offshores),]
plotIndivs <- offshores
plotIndFroh <- froh[match(plotIndivs,froh[,1]),5]
plotIndivs <- plotIndivs[order(plotIndFroh)]
par(xpd=TRUE)
gap <- 10000000


chroms <- unique(ibdTracts[,2])
colVec <- rep(c("blue","orange"),50)




for(i in 1:length(plotIndivs)){
  thisDat <- theseTracts[theseTracts == plotIndivs[i],]
  print(nrow(thisDat))
  xIter <- 0
  for(j in 1:length(chroms)){
    if(chroms[j] %in% thisDat[,2]){
      chrRows <- which(thisDat[,2] == chroms[j])
      for(k in 1:length(chrRows)){
        thisTract <- thisDat[chrRows[k],3:4]
        rect(xleft = xIter + thisTract[1]/1000000,xright=xIter + thisTract[2]/1000000,ybottom=yIter-0.45,ytop=yIter+0.45,col="#33a02c",border=NA)
      }
    }
    xIter <- xIter + gap/1000000 + chrLengs[j]/1000000
  }
  yIter <- yIter + 1
}

# add roh density
xIter <- 0
densPlotMin <- yIter
densPlotMax <- yIter + 20
for(i in 1:length(chroms)){
  chrDat <- winMat_offshore[winMat_offshore[,1] == chroms[i],]
  lines( (as.numeric(chrDat[,3]) - 0.5*winSize)/1000000 + xIter,densPlotMin+(as.numeric(chrDat[,5])/0.05)*20 ,col="#33a02c",lwd=1)
  xIter <- xIter + chrLengs[i]/1000000 + gap/1000000
}
axis(side=2,at=c(yIter,yIter+20),labels=c(0,0.05),pos=-10,cex.axis=0.5)
text(x=-175,y=mean(c(yIter,yIter + 20)),labels ="ROH",cex=0.75,col="#33a02c",srt=90)
text(x=-125,y=mean(c(yIter,yIter + 20)),labels ="Density",cex=0.75,col="#33a02c",srt=90)

yIter <- yIter + 20
rect(xleft=0,xright=sum(chrLengs/1000000) + (gap*20)/1000000,ybottom=yIter+0.5,ytop=yIter + 4,col="white",border="white")
yIter <- yIter + 5

text(x=-225,y=218,labels="Offshore", col="#33a02c",cex=1.0)

################################################
# add roh density above 1 Mb across all populations
################################################
xIter <- 0
densPlotMin <- yIter
densPlotMax <- yIter + 20

for(i in 1:length(chroms)){
  chrDat <- winMat[winMat[,1] == chroms[i],]
  lines( (as.numeric(chrDat[,3]) - 0.5*winSize)/1000000 + xIter,densPlotMin+(as.numeric(chrDat[,5])/0.7)*20 ,col="black",lwd=1)
  xIter <- xIter + chrLengs[i]/1000000 + gap/1000000
}

axis(side=2,at=c(yIter,yIter+20),labels=c(0,0.3),pos=-10,cex.axis=0.5)
text(x=-225,y=mean(c(yIter,yIter + 20)),labels ="Total",cex=0.75,col="black",srt=90)
text(x=-175,y=mean(c(yIter,yIter + 20)),labels ="ROH",cex=0.75,col="black",srt=90)
text(x=-125,y=mean(c(yIter,yIter + 20)),labels ="Density",cex=0.75,col="black",srt=90)

yIter <- yIter + 20
rect(xleft=0,xright=sum(chrLengs/1000000) + (gap*20)/1000000,ybottom=yIter+0.5,ytop=yIter + 4,col="white",border="white")
yIter <- yIter + 5



#-------------------------------------------------------------
# add recombination rates at the top
#-------------------------------------------------------------

setwd("~/Documents/orca/analyses_31May2022/recombinationRate")
# plot rec rates across the whole genome
chroms <- c(1:5,7:22)
fileNums <- 1:10
rates <- NULL     # matrix to store recombination rates
map <- NULL # store the genetic map of the genome
files <- paste("orcaChr1FastEPRR_4_repeatMassked_out.txt")

for(i in 1:length(chroms)){
  chromRates <- NULL
  for(j in 1:length(fileNums)){
    thisFile <- read.table(paste("orcaChr",chroms[i],"FastEPRR_",j,"_repeatMassked_out.txt",sep=""),header=TRUE)
    winCents <- rowMeans(thisFile[,1:2])
    if(j > 1) winCents <- winCents + max(chromRates[,1])
    outRates <- cbind(winCents,thisFile)
    chromRates <- rbind(chromRates,outRates)
  }
  chromRates <- cbind(rep(chroms[i],nrow(chromRates)),chromRates)
  rates <- rbind(rates,chromRates)

  physMb <- chromRates[,2]/1000000
  mapLengs <- (chromRates[,3]/(4*Ne))*physMb

  mapPos <- cumsum(mapLengs)
  chromRates <- cbind(chromRates,mapLengs,mapPos)

  outMap <- chromRates
  map <- rbind(map,outMap)
}

map[map[,5] == 0,5] <- 0.00001

xIter <- 0
smoothVals <- NULL
for(i in 1:length(chroms)){
  chrDat <- map[map[,1] == chroms[i] & map[,2] <= chrLengs[i],]

  mod <- loess(Rho ~ winCents,span = 0.2,data=chrDat)
  smooth <- predict(mod)
  smoothVals <- c(smoothVals,smooth)
  if(sum(smooth < 0) > 0) smooth[smooth < 0] <- 0
  lines(x=xIter + as.numeric(chrDat[,2])/1000000,y=yIter + (smooth/max(smooth))*20,col="black")
  xIter <- xIter + chrLengs[i]/1000000 + gap/1000000
}


axis(side=2,at=c(yIter,yIter+20),labels = c(0,7.9),pos=-10,cex.axis=0.5)
text(x=-150,y=mean(c(yIter,yIter + 20)),labels =expression(rho),cex=1,col="black",srt=90)

yIter <- yIter + 20
rect(xleft=0,xright=sum(chrLengs/1000000) + (gap*20)/1000000,ybottom=yIter+2,ytop=yIter + 50,col="white",border="white")
yIter <- yIter + 5





###############################################################################
###############################################################################
###############################################################################
# same figure with minimum 10 Mb ROH
###############################################################################
###############################################################################
###############################################################################


# analyze the genomic distribution of ROH in killer whales

setwd("~/Documents/orca/analyses_31May2022/ROH3")
ibdTracts <- read.table("allROH_18August2022",header=TRUE)

numChroms <- c(1:5,7:22)

################################################################################
## get the roh data
################################################################################

Froh <- read.table("Froh_7August2022",header=TRUE)
froh <- Froh

popKey <- read.table("popKey.txt",header=TRUE)
minROHSize <- 1000000
winMat <- NULL
allInds <- unique(ibdTracts[,1])
chroms <- unique(ibdTracts[,2])
winSize <- 1000000
chrLengs <- rep(NA,length(chroms))
indWinInbs <- NULL
for(i in 1:length(chroms)){
  thisDat <- ibdTracts[which(ibdTracts[,2] == chroms[i] & (ibdTracts[,4] - ibdTracts[,3]) >= minROHSize),]
  chrLeng <- max(thisDat[,4])
  chrLengs[i] <- chrLeng
  starts <- seq(0,chrLeng,winSize)
  ends <- starts + winSize
  ends[length(ends)] <- chrLeng
  rohCounts <- rep(NA,length(starts))
  winInb <- rep(NA,length(starts))
  chromIndWinInbs <- matrix(0,nrow=length(allInds),ncol=length(starts))
  
  for(j in 1:length(starts)){
    rohCounts[j] <- sum(starts[j] < thisDat[,4] & ends[j] > thisDat[,3])
    
    if(rohCounts[j] > 0){
      theseROH <- thisDat[starts[j] < thisDat[,4] & ends[j] > thisDat[,3],c(1,3,4)]
      if(sum(theseROH[,2] < starts[j]) > 0) theseROH[which(theseROH[,2] < starts[j]),2] <- starts[j]
      if(sum(theseROH[,3] > ends[j]) > 0) theseROH[which(theseROH[,3] > ends[j]),3] <- ends[j]
      
      thisInb <- sum((theseROH[,3] - theseROH[,2])/(ends[j] - starts[j]))/length(unique(ibdTracts[,1]))
      if(thisInb <= 1) winInb[j] <- thisInb
      
      indInbs <- (theseROH[,3] - theseROH[,2])/(ends[j] - starts[j])
      uniqInds <- unique(theseROH[,1])
      for(k in 1:length(uniqInds)){
        chromIndWinInbs[match(uniqInds[k],allInds),j] <- sum(theseROH[match(uniqInds[k],theseROH[,1]),3] - theseROH[match(uniqInds[k],theseROH[,1]),2])/((ends[j] - starts[j])*length(match(uniqInds[k],theseROH[,1])))
      }
    }
  }
  indWinInbs <- cbind(indWinInbs,chromIndWinInbs)
  out <- cbind(rep(chroms[i], length(starts)),starts,ends,rohCounts,winInb)
  winMat <- rbind(winMat,out)
  print(i)
}

winMat[winMat[,4] == "0",5] <- 0


################################################################################
################################################################################
################################################################################
# Plot ROH density for each population along with recombination rate at the top
# minimum ROH length = 1 Mb
################################################################################
################################################################################
################################################################################
# analyze the genomic distribution of ROH in killer whales

ibdTracts <- read.table("allROH_18August2022",header=TRUE)

numChroms <- c(1:5,7:22)

################################################################################
## get the roh data
################################################################################

Froh <- read.table("Froh_7August2022",header=TRUE)


# get the Froh data
popKey <- read.table("popKey.txt",header=TRUE)


#------------------------------------------------------------------------
# get ROH density for all populations and for each population separately
#------------------------------------------------------------------------
minMb <- 10
popKey <- read.table("popKey.txt",header=TRUE)
ibdTracts2 <- ibdTracts[ibdTracts$lengths > minMb,]
srkws <- popKey[which(popKey[,1] == "SRKW"),2]
arkws <- popKey[which(popKey[,1] == "ARKW"),2]
tkws <- popKey[which(popKey[,1] == "TKW"),2]
offshores <- popKey[which(popKey[,1] == "offshore"),2]
nrkws <- popKey[which(popKey[,1] == "NRKW"),2]
minROHSize <- 10000000
winMat <- NULL
allInds <- unique(ibdTracts[,1])
chroms <- unique(ibdTracts[,2])
winSize <- 100000

################ all pops combined
indWinInbs <- NULL
for(i in 1:length(chroms)){
  thisDat <- ibdTracts[which(ibdTracts[,2] == chroms[i] & ibdTracts$length >= minROHSize),]
  chrLeng <-  chrLengs[i]
  starts <- seq(0,chrLeng,winSize)
  ends <- starts + winSize
  ends[length(ends)] <- chrLeng
  rohCounts <- rep(NA,length(starts))
  winInb <- rep(NA,length(starts))
  chromIndWinInbs <- matrix(0,nrow=length(allInds),ncol=length(starts))
  for(j in 1:length(starts)){
    rohCounts[j] <- sum(starts[j] < thisDat[,4] & ends[j] > thisDat[,3])
    if(rohCounts[j] > 0){
      theseROH <- thisDat[starts[j] < thisDat[,4] & ends[j] > thisDat[,3],c(1,3,4)]
      if(sum(theseROH[,2] < starts[j]) > 0) theseROH[which(theseROH[,2] < starts[j]),2] <- starts[j]
      if(sum(theseROH[,3] > ends[j]) > 0) theseROH[which(theseROH[,3] > ends[j]),3] <- ends[j]
      thisInb <- sum((theseROH[,3] - theseROH[,2])/(ends[j] - starts[j]))/length(unique(ibdTracts[,1]))
      if(thisInb <= 1) winInb[j] <- thisInb
      indInbs <- (theseROH[,3] - theseROH[,2])/(ends[j] - starts[j])
      uniqInds <- unique(theseROH[,1])
      for(k in 1:length(uniqInds)){
        chromIndWinInbs[match(uniqInds[k],allInds),j] <- sum(theseROH[match(uniqInds[k],theseROH[,1]),3] - theseROH[match(uniqInds[k],theseROH[,1]),2])/((ends[j] - starts[j])*length(match(uniqInds[k],theseROH[,1])))
      }
    }
  }
  indWinInbs <- cbind(indWinInbs,chromIndWinInbs)
  out <- cbind(rep(chroms[i], length(starts)),starts,ends,rohCounts,winInb)
  winMat <- rbind(winMat,out)
  print(i)
}
winMat[winMat[,4] == "0",5] <- 0

################ SRKWs
winMat_srkw <- NULL
indWinInbs_srkw <- NULL
for(i in 1:length(chroms)){
  thisDat <- ibdTracts[which(ibdTracts[,2] == chroms[i] & (ibdTracts[,4] - ibdTracts[,3]) >= minROHSize & ibdTracts[,1] %in% srkws),]
  chrLeng <- chrLengs[i]
  
  starts <- seq(0,chrLeng,winSize)
  ends <- starts + winSize
  ends[length(ends)] <- chrLeng
  rohCounts <- rep(NA,length(starts))
  winInb <- rep(NA,length(starts))
  chromIndWinInbs <- matrix(0,nrow=length(allInds),ncol=length(starts))
  for(j in 1:length(starts)){
    rohCounts[j] <- sum(starts[j] < thisDat[,4] & ends[j] > thisDat[,3])
    if(rohCounts[j] > 0){
      theseROH <- thisDat[starts[j] < thisDat[,4] & ends[j] > thisDat[,3],c(1,3,4)]
      if(sum(theseROH[,2] < starts[j]) > 0) theseROH[which(theseROH[,2] < starts[j]),2] <- starts[j]
      if(sum(theseROH[,3] > ends[j]) > 0) theseROH[which(theseROH[,3] > ends[j]),3] <- ends[j]
      thisInb <- sum((theseROH[,3] - theseROH[,2])/(ends[j] - starts[j]))/length(unique(ibdTracts[,1]))
      if(thisInb <= 1) winInb[j] <- thisInb
      indInbs <- (theseROH[,3] - theseROH[,2])/(ends[j] - starts[j])
      uniqInds <- unique(theseROH[,1])
      for(k in 1:length(uniqInds)){
        chromIndWinInbs[match(uniqInds[k],allInds),j] <- sum(theseROH[match(uniqInds[k],theseROH[,1]),3] - theseROH[match(uniqInds[k],theseROH[,1]),2])/((ends[j] - starts[j])*length(match(uniqInds[k],theseROH[,1])))
      }
    }
  }
  indWinInbs_srkw <- cbind(indWinInbs_srkw,chromIndWinInbs)
  out <- cbind(rep(chroms[i], length(starts)),starts,ends,rohCounts,winInb)
  winMat_srkw <- rbind(winMat_srkw,out)
  print(i)
}
winMat_srkw[winMat_srkw[,4] == "0",5] <- 0

################ ARKWs
winMat_arkw <- NULL
indWinInbs_arkw <- NULL
for(i in 1:length(chroms)){
  thisDat <- ibdTracts[which(ibdTracts[,2] == chroms[i] & (ibdTracts[,4] - ibdTracts[,3]) >= minROHSize & ibdTracts[,1] %in% arkws),]
  chrLeng <- chrLengs[i]
  starts <- seq(0,chrLeng,winSize)
  ends <- starts + winSize
  ends[length(ends)] <- chrLeng
  rohCounts <- rep(NA,length(starts))
  winInb <- rep(NA,length(starts))
  chromIndWinInbs <- matrix(0,nrow=length(allInds),ncol=length(starts))
  for(j in 1:length(starts)){
    rohCounts[j] <- sum(starts[j] < thisDat[,4] & ends[j] > thisDat[,3])
    if(rohCounts[j] > 0){
      theseROH <- thisDat[starts[j] < thisDat[,4] & ends[j] > thisDat[,3],c(1,3,4)]
      if(sum(theseROH[,2] < starts[j]) > 0) theseROH[which(theseROH[,2] < starts[j]),2] <- starts[j]
      if(sum(theseROH[,3] > ends[j]) > 0) theseROH[which(theseROH[,3] > ends[j]),3] <- ends[j]
      thisInb <- sum((theseROH[,3] - theseROH[,2])/(ends[j] - starts[j]))/length(unique(ibdTracts[,1]))
      if(thisInb <= 1) winInb[j] <- thisInb
      indInbs <- (theseROH[,3] - theseROH[,2])/(ends[j] - starts[j])
      uniqInds <- unique(theseROH[,1])
      for(k in 1:length(uniqInds)){
        chromIndWinInbs[match(uniqInds[k],allInds),j] <- sum(theseROH[match(uniqInds[k],theseROH[,1]),3] - theseROH[match(uniqInds[k],theseROH[,1]),2])/((ends[j] - starts[j])*length(match(uniqInds[k],theseROH[,1])))
      }
    }
  }
  indWinInbs_arkw <- cbind(indWinInbs_arkw,chromIndWinInbs)
  out <- cbind(rep(chroms[i], length(starts)),starts,ends,rohCounts,winInb)
  winMat_arkw <- rbind(winMat_arkw,out)
  print(i)
}
winMat_arkw[winMat_arkw[,4] == "0",5] <- 0

################ offshores
winMat_offshore <- NULL
indWinInbs_offshore <- NULL
for(i in 1:length(chroms)){
  thisDat <- ibdTracts[which(ibdTracts[,2] == chroms[i] & (ibdTracts[,4] - ibdTracts[,3]) >= minROHSize & ibdTracts[,1] %in% offshores),]
  chrLeng <- chrLengs[i]
  
  starts <- seq(0,chrLeng,winSize)
  ends <- starts + winSize
  ends[length(ends)] <- chrLeng
  rohCounts <- rep(NA,length(starts))
  winInb <- rep(NA,length(starts))
  chromIndWinInbs <- matrix(0,nrow=length(allInds),ncol=length(starts))
  for(j in 1:length(starts)){
    rohCounts[j] <- sum(starts[j] < thisDat[,4] & ends[j] > thisDat[,3])
    if(rohCounts[j] > 0){
      theseROH <- thisDat[starts[j] < thisDat[,4] & ends[j] > thisDat[,3],c(1,3,4)]
      if(sum(theseROH[,2] < starts[j]) > 0) theseROH[which(theseROH[,2] < starts[j]),2] <- starts[j]
      if(sum(theseROH[,3] > ends[j]) > 0) theseROH[which(theseROH[,3] > ends[j]),3] <- ends[j]
      thisInb <- sum((theseROH[,3] - theseROH[,2])/(ends[j] - starts[j]))/length(unique(ibdTracts[,1]))
      if(thisInb <= 1) winInb[j] <- thisInb
      indInbs <- (theseROH[,3] - theseROH[,2])/(ends[j] - starts[j])
      uniqInds <- unique(theseROH[,1])
      for(k in 1:length(uniqInds)){
        chromIndWinInbs[match(uniqInds[k],allInds),j] <- sum(theseROH[match(uniqInds[k],theseROH[,1]),3] - theseROH[match(uniqInds[k],theseROH[,1]),2])/((ends[j] - starts[j])*length(match(uniqInds[k],theseROH[,1])))
      }
    }
  }
  indWinInbs_offshore <- cbind(indWinInbs_offshore,chromIndWinInbs)
  out <- cbind(rep(chroms[i], length(starts)),starts,ends,rohCounts,winInb)
  winMat_offshore <- rbind(winMat_offshore,out)
  print(i)
}
winMat_offshore[winMat_offshore[,4] == "0",5] <- 0


################ TKWs
winMat_tkw <- NULL
indWinInbs_tkw <- NULL
for(i in 1:length(chroms)){
  thisDat <- ibdTracts[which(ibdTracts[,2] == chroms[i] & (ibdTracts[,4] - ibdTracts[,3]) >= minROHSize & ibdTracts[,1] %in% tkws),]
  chrLeng <- chrLengs[i]
  
  starts <- seq(0,chrLeng,winSize)
  ends <- starts + winSize
  ends[length(ends)] <- chrLeng
  rohCounts <- rep(NA,length(starts))
  winInb <- rep(NA,length(starts))
  chromIndWinInbs <- matrix(0,nrow=length(allInds),ncol=length(starts))
  for(j in 1:length(starts)){
    rohCounts[j] <- sum(starts[j] < thisDat[,4] & ends[j] > thisDat[,3])
    if(rohCounts[j] > 0){
      theseROH <- thisDat[starts[j] < thisDat[,4] & ends[j] > thisDat[,3],c(1,3,4)]
      if(sum(theseROH[,2] < starts[j]) > 0) theseROH[which(theseROH[,2] < starts[j]),2] <- starts[j]
      if(sum(theseROH[,3] > ends[j]) > 0) theseROH[which(theseROH[,3] > ends[j]),3] <- ends[j]
      
      thisInb <- sum((theseROH[,3] - theseROH[,2])/(ends[j] - starts[j]))/length(unique(ibdTracts[,1]))
      if(thisInb <= 1) winInb[j] <- thisInb
      
      indInbs <- (theseROH[,3] - theseROH[,2])/(ends[j] - starts[j])
      uniqInds <- unique(theseROH[,1])
      for(k in 1:length(uniqInds)){
        chromIndWinInbs[match(uniqInds[k],allInds),j] <- sum(theseROH[match(uniqInds[k],theseROH[,1]),3] - theseROH[match(uniqInds[k],theseROH[,1]),2])/((ends[j] - starts[j])*length(match(uniqInds[k],theseROH[,1])))
      }
    }
  }
  indWinInbs_tkw <- cbind(indWinInbs_tkw,chromIndWinInbs)
  out <- cbind(rep(chroms[i], length(starts)),starts,ends,rohCounts,winInb)
  winMat_tkw <- rbind(winMat_tkw,out)
  print(i)
}
winMat_tkw[winMat_tkw[,4] == "0",5] <- 0


#############################################################################################################
# plot ROH and ROH density for each population across the whole genome with density for each population
#############################################################################################################

minMb <- 10000000


ibdTracts2 <- ibdTracts[ibdTracts[,5] > minMb,]
srkws <- popKey[which(popKey[,1] == "SRKW"),2]
srkwFroh <- froh[match(srkws,froh[,1]),3]
srkws <- srkws[order(srkwFroh)]

arkws <- popKey[which(popKey[,1] == "ARKW"),2]
arkwFroh <- froh[match(arkws,froh[,1]),3]
arkws <- arkws[order(arkwFroh)]

offshores <- popKey[which(popKey[,1] == "offshore"),2]
offshoreFroh <- froh[match(offshores,froh[,1]),3]
offshores <- offshores[order(offshoreFroh)]

tkws <- popKey[which(popKey[,1] == "TKW"),2]
tkwFroh <- froh[match(tkws,froh[,1]),3]
tkws <- tkws[order(tkwFroh)]

nrkws <- popKey[which(popKey[,1] == "NRKW"),2]
nrkwFroh <- froh[match(nrkws,froh[,1]),3]
nrkws <- nrkws[order(nrkwFroh)]




par(mar=c(5,7,2,2))
#####################
# initialize the plot
#####################
gap <- 10000000
plot(c(0,sum(chrLengs/1000000) + (gap*20)/1000000),c(0,149 +157),type="n",
     ylab="",xlab="Chromosome Position",axes=FALSE,cex.lab=1.5)

#######----------------------------
# SRKWs first
#######----------------------------
yIter <- 1
theseTracts <- ibdTracts2[which(ibdTracts2[,1] %in% srkws),]
plotIndivs <- srkws[-which(srkws == "J32.fetus")]
plotIndFroh <- froh[match(plotIndivs,froh[,1]),3]
plotIndivs <- plotIndivs[order(plotIndFroh)]
par(xpd=TRUE)



chroms <- unique(ibdTracts[,2])
colVec <- rep(c("black","darkgray"),50)
xIter <- 0
for(i in 1:length(chroms)){
  rect(xleft=xIter,xright=xIter + chrLengs[i]/1000000,ybottom=yIter,ytop=149 +157 ,col="gray97",border="gray97")
  text(x=mean(c(xIter,xIter + chrLengs[i]/1000000)),y=-8,labels=c(1:5,7:22)[i])
  xIter <- xIter + gap/1000000 + chrLengs[i]/1000000
}



for(i in 1:length(plotIndivs)){
  thisDat <- theseTracts[theseTracts == plotIndivs[i],]
  xIter <- 0
  for(j in 1:length(chroms)){
    if(chroms[j] %in% thisDat[,2]){
      chrRows <- which(thisDat[,2] == chroms[j])
      for(k in 1:length(chrRows)){
        thisTract <- thisDat[chrRows[k],3:4]
        rect(xleft = xIter + thisTract[1]/1000000,xright=xIter + thisTract[2]/1000000,ybottom=yIter-0.45,ytop=yIter+0.45,col="#a6cef9",border=NA)
      }
    }
    xIter <- xIter + gap/1000000 + chrLengs[j]/1000000
  }
  yIter <- yIter + 1
}

# add roh density

xIter <- 0
densPlotMin <- yIter
densPlotMax <- yIter + 20
for(i in 1:length(chroms)){
  chrDat <- winMat_srkw[winMat_srkw[,1] == chroms[i],]
  lines( (as.numeric(chrDat[,3]) - 0.5*winSize)/1000000 + xIter,densPlotMin+(as.numeric(chrDat[,5])/0.15)*20 ,col="#a6cef9",lwd=1)
  xIter <- xIter + chrLengs[i]/1000000 + gap/1000000
}
axis(side=2,at=c(yIter,yIter+20),labels=c(0,0.15),pos=-10,cex.axis=0.5)
text(x=-175,y=mean(c(yIter,yIter + 20)),labels ="ROH",cex=0.75,col="#a6cef9",srt=90)
text(x=-125,y=mean(c(yIter,yIter + 20)),labels ="Density",cex=0.75,col="#a6cef9",srt=90)

yIter <- yIter + 20
rect(xleft=0,xright=sum(chrLengs/1000000) + (gap*20)/1000000,ybottom=yIter+0.5,ytop=yIter + 4,col="white",border="white")
yIter <- yIter + 5



text(x=-225,y=51,labels="Southern Residents", col="#a6cef9",cex=1.0)


#######----------------------------
# ARKWs
#######----------------------------
theseTracts <- ibdTracts2[which(ibdTracts2[,1] %in% arkws),]
plotIndivs <- arkws
plotIndFroh <- froh[match(plotIndivs,froh[,1]),5]
plotIndivs <- plotIndivs[order(plotIndFroh)]
par(xpd=TRUE)
gap <- 10000000
chroms <- unique(ibdTracts[,2])
colVec <- rep(c("blue","orange"),50)

for(i in 1:length(arkws)){
  if(plotIndivs[i] %in% theseTracts[,1]){
    thisDat <- theseTracts[theseTracts == plotIndivs[i],]
    xIter <- 0
    for(j in 1:length(chroms)){
      if(chroms[j] %in% thisDat[,2]){
        chrRows <- which(thisDat[,2] == chroms[j])
        for(k in 1:length(chrRows)){
          thisTract <- thisDat[chrRows[k],3:4]
          rect(xleft = xIter + thisTract[1]/1000000,xright=xIter + thisTract[2]/1000000,ybottom=yIter-0.45,ytop=yIter+0.45,col="#1f78b4",border=NA)
        }
      }
      xIter <- xIter + gap/1000000 + chrLengs[j]/1000000
    }
  }
  yIter <- yIter + 1
}



# add roh density
xIter <- 0
densPlotMin <- yIter
densPlotMax <- yIter + 20
yRef <- max(as.numeric(winMat[,4]))
for(i in 1:length(chroms)){
  chrDat <- winMat_arkw[winMat_arkw[,1] == chroms[i],]
  lines( (as.numeric(chrDat[,3]) - 0.5*winSize)/1000000 + xIter,densPlotMin+(as.numeric(chrDat[,5])/0.014)*20 ,col="#1f78b4",lwd=1)
  xIter <- xIter + chrLengs[i]/1000000 + gap/1000000
}
axis(side=2,at=c(yIter,yIter+20),labels=c(0,0.014),pos=-10,cex.axis=0.5)
text(x=-175,y=mean(c(yIter,yIter + 20)),labels ="ROH",cex=0.75,col="#1f78b4",srt=90)
text(x=-125,y=mean(c(yIter,yIter + 20)),labels ="Density",cex=0.75,col="#1f78b4",srt=90)

yIter <- yIter + 20
rect(xleft=0,xright=sum(chrLengs/1000000) + (gap*20)/1000000,ybottom=yIter+0.5,ytop=yIter + 4,col="white",border="white")
yIter <- yIter + 5



text(x=-225,y=137.5,labels="Alaska Residents", col="#1f78b4",cex=1.0)

#######----------------------------
# TKWs
#######----------------------------
theseTracts <- ibdTracts2[which(ibdTracts2[,1] %in% tkws),]
plotIndivs <- tkws
plotIndFroh <- froh[match(plotIndivs,froh[,1]),5]
plotIndivs <- plotIndivs[order(plotIndFroh)]
par(xpd=TRUE)
gap <- 10000000


chroms <- unique(ibdTracts[,2])
colVec <- rep(c("blue","orange"),50)



for(i in 1:length(plotIndivs)){
  if(plotIndivs[i] %in% theseTracts[,1]){
    thisDat <- theseTracts[theseTracts == plotIndivs[i],]
    xIter <- 0
    for(j in 1:length(chroms)){
      if(chroms[j] %in% thisDat[,2]){
        chrRows <- which(thisDat[,2] == chroms[j])
        for(k in 1:length(chrRows)){
          thisTract <- thisDat[chrRows[k],3:4]
          rect(xleft = xIter + thisTract[1]/1000000,xright=xIter + thisTract[2]/1000000,ybottom=yIter-0.45,ytop=yIter+0.45,col="#b2df8a",border=NA)
        }
      }
      xIter <- xIter + gap/1000000 + chrLengs[j]/1000000
    }
  }
  yIter <- yIter + 1
}



# add roh density
xIter <- 0
densPlotMin <- yIter
densPlotMax <- yIter + 20
for(i in 1:length(chroms)){
  chrDat <- winMat_tkw[winMat_tkw[,1] == chroms[i],]
  lines( (as.numeric(chrDat[,3]) - 0.5*winSize)/1000000 + xIter,densPlotMin+(as.numeric(chrDat[,5])/0.03)*20 ,col="#b2df8a",lwd=1)
  xIter <- xIter + chrLengs[i]/1000000 + gap/1000000
}
axis(side=2,at=c(yIter,yIter+20),labels=c(0,0.03),pos=-10,cex.axis=0.5)
text(x=-175,y=mean(c(yIter,yIter + 20)),labels ="ROH",cex=0.75,col="#b2df8a",srt=90)
text(x=-125,y=mean(c(yIter,yIter + 20)),labels ="Density",cex=0.75,col="#b2df8a",srt=90)

yIter <- yIter + 20
rect(xleft=0,xright=sum(chrLengs/1000000) + (gap*20)/1000000,ybottom=yIter+0.5,ytop=yIter + 4,col="white",border="white")
yIter <- yIter + 5




text(x=-225,y=183,labels="Transients", col="#b2df8a",cex=1.0)


#######----------------------------
# offshores
#######----------------------------
theseTracts <- ibdTracts2[which(ibdTracts2[,1] %in% offshores),]
plotIndivs <- offshores
plotIndFroh <- froh[match(plotIndivs,froh[,1]),5]
plotIndivs <- plotIndivs[order(plotIndFroh)]
par(xpd=TRUE)
gap <- 10000000


chroms <- unique(ibdTracts[,2])
colVec <- rep(c("blue","orange"),50)




for(i in 1:length(plotIndivs)){
  thisDat <- theseTracts[theseTracts == plotIndivs[i],]
  print(nrow(thisDat))
  xIter <- 0
  for(j in 1:length(chroms)){
    if(chroms[j] %in% thisDat[,2]){
      chrRows <- which(thisDat[,2] == chroms[j])
      for(k in 1:length(chrRows)){
        thisTract <- thisDat[chrRows[k],3:4]
        rect(xleft = xIter + thisTract[1]/1000000,xright=xIter + thisTract[2]/1000000,ybottom=yIter-0.45,ytop=yIter+0.45,col="#33a02c",border=NA)
      }
    }
    xIter <- xIter + gap/1000000 + chrLengs[j]/1000000
  }
  yIter <- yIter + 1
}

# add roh density
xIter <- 0
densPlotMin <- yIter
densPlotMax <- yIter + 20
for(i in 1:length(chroms)){
  chrDat <- winMat_offshore[winMat_offshore[,1] == chroms[i],]
  lines( (as.numeric(chrDat[,3]) - 0.5*winSize)/1000000 + xIter,densPlotMin+(as.numeric(chrDat[,5])/0.014)*20 ,col="#33a02c",lwd=1)
  xIter <- xIter + chrLengs[i]/1000000 + gap/1000000
}
axis(side=2,at=c(yIter,yIter+20),labels=c(0,0.014),pos=-10,cex.axis=0.5)
text(x=-175,y=mean(c(yIter,yIter + 20)),labels ="ROH",cex=0.75,col="#33a02c",srt=90)
text(x=-125,y=mean(c(yIter,yIter + 20)),labels ="Density",cex=0.75,col="#33a02c",srt=90)

yIter <- yIter + 20
rect(xleft=0,xright=sum(chrLengs/1000000) + (gap*20)/1000000,ybottom=yIter+0.5,ytop=yIter + 4,col="white",border="white")
yIter <- yIter + 5

text(x=-225,y=218,labels="Offshore", col="#33a02c",cex=1.0)

################################################
# add roh density above 1 Mb across all populations
################################################
xIter <- 0
densPlotMin <- yIter
densPlotMax <- yIter + 20

for(i in 1:length(chroms)){
  chrDat <- winMat[winMat[,1] == chroms[i],]
  lines( (as.numeric(chrDat[,3]) - 0.5*winSize)/1000000 + xIter,densPlotMin+(as.numeric(chrDat[,5])/0.17)*20 ,col="black",lwd=1)
  xIter <- xIter + chrLengs[i]/1000000 + gap/1000000
}

axis(side=2,at=c(yIter,yIter+20),labels=c(0,0.17),pos=-10,cex.axis=0.5)
text(x=-225,y=mean(c(yIter,yIter + 20)),labels ="Total",cex=0.75,col="black",srt=90)
text(x=-175,y=mean(c(yIter,yIter + 20)),labels ="ROH",cex=0.75,col="black",srt=90)
text(x=-125,y=mean(c(yIter,yIter + 20)),labels ="Density",cex=0.75,col="black",srt=90)

yIter <- yIter + 20
rect(xleft=0,xright=sum(chrLengs/1000000) + (gap*20)/1000000,ybottom=yIter+0.5,ytop=yIter + 4,col="white",border="white")
yIter <- yIter + 5



#-------------------------------------------------------------
# add recombination rates at the top
#-------------------------------------------------------------

setwd("~/Documents/orca/analyses_31May2022/recombinationRate")
# plot rec rates across the whole genome
chroms <- c(1:5,7:22)
fileNums <- 1:10
rates <- NULL     # matrix to store recombination rates
map <- NULL # store the genetic map of the genome
files <- paste("orcaChr1FastEPRR_4_repeatMassked_out.txt")

for(i in 1:length(chroms)){
  chromRates <- NULL
  for(j in 1:length(fileNums)){
    thisFile <- read.table(paste("orcaChr",chroms[i],"FastEPRR_",j,"_repeatMassked_out.txt",sep=""),header=TRUE)
    winCents <- rowMeans(thisFile[,1:2])
    if(j > 1) winCents <- winCents + max(chromRates[,1])
    outRates <- cbind(winCents,thisFile)
    chromRates <- rbind(chromRates,outRates)
  }
  chromRates <- cbind(rep(chroms[i],nrow(chromRates)),chromRates)
  rates <- rbind(rates,chromRates)
  
  physMb <- chromRates[,2]/1000000
  mapLengs <- (chromRates[,3]/(4*Ne))*physMb
  
  mapPos <- cumsum(mapLengs)
  chromRates <- cbind(chromRates,mapLengs,mapPos)
  
  outMap <- chromRates
  map <- rbind(map,outMap)
}

map[map[,5] == 0,5] <- 0.00001

xIter <- 0
smoothVals <- NULL
for(i in 1:length(chroms)){
  chrDat <- map[map[,1] == chroms[i] & map[,2] <= chrLengs[i],]
  mod <- loess(Rho ~ winCents,span = 0.2,data=chrDat)
  smooth <- predict(mod)
  smoothVals <- c(smoothVals,smooth)
  if(sum(smooth < 0) > 0) smooth[smooth < 0] <- 0
  lines(x=xIter + as.numeric(chrDat[,2])/1000000,y=yIter + (smooth/max(smooth))*20,col="black")
  xIter <- xIter + chrLengs[i]/1000000 + gap/1000000
}


axis(side=2,at=c(yIter,yIter+20),labels = c(0,7.9),pos=-10,cex.axis=0.5)
text(x=-150,y=mean(c(yIter,yIter + 20)),labels =expression(rho),cex=1,col="black",srt=90)

yIter <- yIter + 20
rect(xleft=0,xright=sum(chrLengs/1000000) + (gap*20)/1000000,ybottom=yIter+2,ytop=yIter + 50,col="white",border="white")
yIter <- yIter + 5


