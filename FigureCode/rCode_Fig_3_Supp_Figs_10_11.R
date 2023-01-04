################################################################################
# plot deleterious homozygous count versus Froh for the different populations
################################################################################
setwd("~/Documents/orca/analyses_31May2022/RXY")
Froh <- read.table("Froh_7August2022",header=TRUE)
popKey <- read.table("popKey.txt",header=TRUE)
froh <- Froh
het <- read.table("orcaHeterozygosity_17August20221",header=TRUE)
colVec <-  c("#a6cef9","#1f78b4","#33a02c","#b2df8a")
popID <- c("SRKW","ARKW","Offshore","TKW")
fixedSeg <- read.table("numDelFixedAndSegregating_allPops",header=TRUE)

froh[which(froh[,1] == "NGOS09.11"),1] <- "NGOS09-11"
#######################################################################
# plot homozygous load versus Froh
#######################################################################
tkw_delHoms <- read.table("deleteriousHomozygousCount_tkw",header=TRUE)
arkw_delHoms <- read.table("deleteriousHomozygousCount_arkw",header=TRUE)
srkw_delHoms <- read.table("deleteriousHomozygousCount_srkw",header=TRUE)
offshore_delHoms <- read.table("deleteriousHomozygousCount_offshore",header=TRUE)
nrkw_delHoms <- read.table("deleteriousHomozygousCount_nrkw",header=TRUE)


arkw_delHoms <- cbind(arkw_delHoms,froh[match(arkw_delHoms[,1],froh[,1]),3:4])
srkw_delHoms <- cbind(srkw_delHoms,froh[match(srkw_delHoms[,1],froh[,1]),3:4])
tkw_delHoms <- cbind(tkw_delHoms,froh[match(tkw_delHoms[,1],froh[,1]),3:4])
nrkw_delHoms <- cbind(nrkw_delHoms,froh[match(nrkw_delHoms[,1],froh[,1]),3:4])
offshore_delHoms <- cbind(offshore_delHoms,froh[match(offshore_delHoms[,1],froh[,1]),3:4])

arkw_delHoms <- cbind(arkw_delHoms,het[match(arkw_delHoms[,1],het[,1]),2])
srkw_delHoms <- cbind(srkw_delHoms,het[match(srkw_delHoms[,1],het[,1]),2])
tkw_delHoms <- cbind(tkw_delHoms,het[match(tkw_delHoms[,1],het[,1]),2])
nrkw_delHoms <- cbind(nrkw_delHoms,het[match(nrkw_delHoms[,1],het[,1]),2])
offshore_delHoms <- cbind(offshore_delHoms,het[match(offshore_delHoms[,1],het[,1]),2])

colnames(arkw_delHoms)[3:5] <- c("Froh_1Mb","Froh_10Mb","het")
colnames(tkw_delHoms)[3:5] <- c("Froh_1Mb","Froh_10Mb","het")
colnames(srkw_delHoms)[3:5] <- c("Froh_1Mb","Froh_10Mb","het")
colnames(offshore_delHoms)[3:5] <- c("Froh_1Mb","Froh_10Mb","het")
colnames(nrkw_delHoms)[3:5] <- c("Froh_1Mb","Froh_10Mb","het")

allHoms <- rbind(tkw_delHoms,arkw_delHoms,srkw_delHoms,offshore_delHoms,nrkw_delHoms)
yLims <- c(min(as.numeric(allHoms[,2])),max(as.numeric(allHoms[,2])))
xLims <- c(min(froh[,3]),max(froh[,3]))

#------------------------------------------
# add Rxy
#------------------------------------------

par(mar=c(5,6,2,2),mfrow=c(1,3))
missRxyEsts <- read.table("missRxyEsts")
lofRxyEsts <- read.table("lofRxyEsts")
missRxySDs <- read.table("missRxyCIs")
lofRxySDs <- read.table("lofRxyCIs")
par(xpd=TRUE)

#LOF
plot(c(0.5,3.5),c(0.9,1.2),type="n",axes=FALSE,xlab="Population Pair",
     ylab=expression(italic(""*R*"")["X/Y"]),cex.lab=1.8)
text(x=-.25,y=1.22,labels="A",cex=2)
starts <- seq(1.5,3.5,1)
for (i in 1:3){
  rect(xleft=starts[i] - 0.99,xright=starts[i] - 0.01,ybottom=0.9,ytop=1.2,col="gray98",border=FALSE)
}
lines(c(0.5,3.5),c(1,1),lty="dashed")
points((1:3)-0.1,lofRxyEsts[,1],pch=16,col="darkred",cex=2)
points((1:3)+0.1,missRxyEsts[,1],pch=16,col="darkgray",cex=2)

for(i in 1:3){
  arrows(x0=i-0.1,x1=i-0.1,y0=lofRxyEsts[i,1],y1=lofRxySDs[i,1],col="darkred",angle=90,length=0.1)
  arrows(x0=i-0.1,x1=i-0.1,y0=lofRxyEsts[i,1],y1=lofRxySDs[i,2],col="darkred",angle=90,length=0.1)
  arrows(x0=i+0.1,x1=i+0.1,y0=missRxyEsts[i,1],y1=missRxySDs[i,1],col="darkgray",angle=90,length=0.1)
  arrows(x0=i+0.1,x1=i+0.1,y0=missRxyEsts[i,1],y1=missRxySDs[i,2],col="darkgray",angle=90,length=0.1)
}

axis(side=2,at=seq(0.9,1.2,0.05))
axis(side=1,at=seq(0.5,3.5,1),labels=rep("",4))
text(x=1,y=0.87,labels="SRKW / ARKW")
text(x=2,y=0.87,labels="SRKW / TKW")
text(x=3,y=0.87,labels="ARKW / TKW")

legend(x=0.5,y=1.15,legend=c("Loss of function","Missense"),col=c("darkred","darkgray"),pch=16,cex=1.3,
       xjust=FALSE,yjust=FALSE,bty="n")

#------------------------------------------------
# make the plot of homozygous load versus Froh
#------------------------------------------------
plot(c(0,0.5),c(2800,4800),type="n",ylab="Homozygous mutation load",
     xlab=expression(italic(""*F*"")["ROH,1Mb"]),cex.lab=1.8,axes=FALSE)
axis(side=1,at=seq(0,0.5,0.1))
axis(side=2,at=seq(2800,4800,500))

library(scales)
text(x=-0.12,y=4900,labels="B",cex=2)
points(srkw_delHoms[,3],srkw_delHoms[,2],bg=alpha("#a6cee3",alpha=0.8),cex=2,pch=21)
points(arkw_delHoms[,3],arkw_delHoms[,2],bg=alpha("#1f78b4",alpha=0.8),cex=2,pch=21)
points(tkw_delHoms[,3],tkw_delHoms[,2],bg=alpha("#b2df8a",alpha=0.8),cex=2,pch=21)
points(nrkw_delHoms[,3],nrkw_delHoms[,2],bg=alpha( "#fb9a99",alpha=0.8),cex=2,pch=21)
points(offshore_delHoms[,3],offshore_delHoms[,2],bg=alpha("#33a02c",alpha=0.8),cex=2,pch=21)
legend(x=0,y=4250,xjust=FALSE,yjust=FALSE,bty="n",legend=c("Southern Residents","Northern Residents","Alaska Residents","Offshore","Transients"),
       pch=16,col= c("#a6cef9","#fb9a99","#1f78b4","#33a02c","#b2df8a"),cex=1.1)
legend(x=0,y=4250,xjust=FALSE,yjust=FALSE,bty="n",legend=c("Southern Residents","Northern Residents","Alaska Residents","Offshore","Transients"),
       pch=1,cex=1.1)

# regression models
# TKW
modDat <- tkw_delHoms
mod <- lm(modDat[,2]~modDat[,3])
int <- mod$coefficients[1]
slop <- mod$coefficients[2]
xRange <- c(min(modDat[,3],na.rm=TRUE),max(modDat[,3],na.rm=TRUE))
yRange <- c(int+slop*xRange[1],int+slop*xRange[2])
lines(x=xRange,y=yRange,col="orange",lwd=5)
lines(x=xRange,y=yRange,col="#b2df8a",lwd=3)
lines(x=xRange,y=yRange,col="black",lwd=3,lty="dashed")

# offshore
modDat <- offshore_delHoms
mod <- lm(modDat[,2]~modDat[,3])
# model is not statistically significant for the offshores

# ARKW
modDat <- arkw_delHoms
mod <- lm(modDat[,2]~modDat[,3])
int <- mod$coefficients[1]
slop <- mod$coefficients[2]
xRange <- c(min(modDat[,3],na.rm=TRUE),max(modDat[,3],na.rm=TRUE))
yRange <- c(int+slop*xRange[1],int+slop*xRange[2])
lines(x=xRange,y=yRange,col="orange",lwd=5)
lines(x=xRange,y=yRange,col="#1f78b4",lwd=3)
lines(x=xRange,y=yRange,col="black",lwd=3,lty="dashed")

# SRKW
modDat <- srkw_delHoms
mod <- lm(modDat[,2]~modDat[,3])
int <- mod$coefficients[1]
slop <- mod$coefficients[2]
xRange <- c(min(modDat[,3],na.rm=TRUE),max(modDat[,3],na.rm=TRUE))
yRange <- c(int+slop*xRange[1],int+slop*xRange[2])
lines(x=xRange,y=yRange,col="orange",lwd=5)
lines(x=xRange,y=yRange,col="#a6cef9",lwd=3)
lines(x=xRange,y=yRange,col="black",lwd=3,lty="dashed")

#-----------------------------------------------------------------------
# plot the homozygous load due to segregating alleles, then the total
# all in one panel
#-----------------------------------------------------------------------
# homozygous load
yLims <- c(0,5000)
plot(c(0,4),yLims,type="n",axes=FALSE,ylab="Homozygous mutation load",
     xlab="Population",cex.lab=1.6)

homLoadList <- list()
homLoadList[[1]] <- srkw_delHoms[,2]
homLoadList[[2]] <- arkw_delHoms[,2]
homLoadList[[3]] <- offshore_delHoms[,2]
homLoadList[[4]] <- tkw_delHoms[,2]

axis(side=1,at=seq(0,4),labels=rep("",5))
text(x=(1:4)-0.5,y=min(yLims)-rep(0.1*(yLims[2] - yLims[1]),5),labels=popID)
axis(side=2,at=seq(0,5000,1000))

# segregating load
homLoadList <- list()
homLoadList[[1]] <- srkw_delHoms[,2]
homLoadList[[2]] <- arkw_delHoms[,2]
homLoadList[[3]] <- offshore_delHoms[,2]
homLoadList[[4]] <- tkw_delHoms[,2]

for(i in 1:4){
  thisDat <- homLoadList[[i]]
  xVals <- runif(length(thisDat),min=i-0.36,max=i-0.24)
  yVals1 <- fixedSeg[match(popID[i],fixedSeg[,1]),3]
  yVals2 <- mean(homLoadList[[i]])
  rect(xleft=i-0.95,xright=i-0.05,ybottom=yVals1,ytop=yVals2,col=alpha(colVec[i],alpha=0.5))
  rect(xleft=i-0.95,xright=i-0.05,ybottom=0,ytop=yVals1,col=alpha(colVec[i],alpha=0.5))
  rect(xleft=i-0.95,xright=i-0.05,ybottom=0,ytop=yVals1,density=20,angle=45,border=TRUE)
}
legend(x=1.5,y=4300,xjust=FALSE,yjust=FALSE,col=c("darkgray","darkgray"),angle=c(45,45),
       density = c(0,40),legend=c("Polymorphic loci","Fixed loci"),bty="n",cex=1.2)

# points for all homozygous load
for(i in 1:4){
  thisDat <- homLoadList[[i]]
  xVals <- runif(length(thisDat),min=i-0.75,max=i-0.25)
  points(xVals,thisDat,bg=alpha(colVec[i],alpha=0.4),pch=21,cex=1.5)
  arrows(x0=i-0.5, y0=mean(thisDat), x1 = i-0.5, y1 = mean(thisDat)+sd(thisDat), length = 0.1, angle = 90,col="orange",lwd=2)
  arrows(x0=i-0.5, y0=mean(thisDat), x1 = i-0.5, y1 = mean(thisDat)-sd(thisDat), length = 0.1, angle = 90,col="orange",lwd=2)
  points(x=i-0.5,y=mean(thisDat),col="orange",cex=1.5,pch=16)
  print(length(thisDat))
}
text(x=-1,y=5280,labels="C",cex=2)







#-----------------------------------------------------------------------
########################################################################
# deleterious homozygous count for moderate and LOF mutations separately
########################################################################
#-----------------------------------------------------------------------

#######################################################################
# plot homozygous load versus Froh
#######################################################################
tkw_delHoms <- read.table("lofHomozygousCount_tkw",header=TRUE)
arkw_delHoms <- read.table("lofHomozygousCount_arkw",header=TRUE)
srkw_delHoms <- read.table("lofHomozygousCount_srkw",header=TRUE)
offshore_delHoms <- read.table("lofHomozygousCount_offshore",header=TRUE)
nrkw_delHoms <- read.table("lofHomozygousCount_nrkw",header=TRUE)

arkw_delHoms <- cbind(arkw_delHoms,froh[match(arkw_delHoms[,1],froh[,1]),3:4])
srkw_delHoms <- cbind(srkw_delHoms,froh[match(srkw_delHoms[,1],froh[,1]),3:4])
tkw_delHoms <- cbind(tkw_delHoms,froh[match(tkw_delHoms[,1],froh[,1]),3:4])
nrkw_delHoms <- cbind(nrkw_delHoms,froh[match(nrkw_delHoms[,1],froh[,1]),3:4])
offshore_delHoms <- cbind(offshore_delHoms,froh[match(offshore_delHoms[,1],froh[,1]),3:4])

arkw_delHoms <- cbind(arkw_delHoms,het[match(arkw_delHoms[,1],het[,1]),2])
srkw_delHoms <- cbind(srkw_delHoms,het[match(srkw_delHoms[,1],het[,1]),2])
tkw_delHoms <- cbind(tkw_delHoms,het[match(tkw_delHoms[,1],het[,1]),2])
nrkw_delHoms <- cbind(nrkw_delHoms,het[match(nrkw_delHoms[,1],het[,1]),2])
offshore_delHoms <- cbind(offshore_delHoms,het[match(offshore_delHoms[,1],het[,1]),2])

colnames(arkw_delHoms)[3:5] <- c("Froh_1Mb","Froh_10Mb","het")
colnames(tkw_delHoms)[3:5] <- c("Froh_1Mb","Froh_10Mb","het")
colnames(srkw_delHoms)[3:5] <- c("Froh_1Mb","Froh_10Mb","het")
colnames(offshore_delHoms)[3:5] <- c("Froh_1Mb","Froh_10Mb","het")
colnames(nrkw_delHoms)[3:5] <- c("Froh_1Mb","Froh_10Mb","het")

allHoms <- rbind(tkw_delHoms,arkw_delHoms,srkw_delHoms,offshore_delHoms,nrkw_delHoms)
yLims <- c(min(as.numeric(allHoms[,2])),max(as.numeric(allHoms[,2])))
xLims <- c(min(froh[,3]),max(froh[,3]))

par(mfrow=c(1,1),mar=c(5,6,2,2))
#------------------------------------------------
# make the plot of homozygous load versus Froh
#------------------------------------------------
plot(c(0,0.5),c(60,160),type="n",ylab="Homozygous mutation load",
     xlab=expression(italic(""*F*"")["ROH,1Mb"]),cex.lab=1.8,axes=FALSE)
axis(side=1,at=seq(0,0.5,0.1))
axis(side=2,at=seq(60,160,20))

library(scales)
text(x=-0.12,y=4900,labels="B",cex=2)
points(srkw_delHoms[,3],srkw_delHoms[,2],bg=alpha("#a6cee3",alpha=0.8),cex=2,pch=21)
points(arkw_delHoms[,3],arkw_delHoms[,2],bg=alpha("#1f78b4",alpha=0.8),cex=2,pch=21)
points(tkw_delHoms[,3],tkw_delHoms[,2],bg=alpha("#b2df8a",alpha=0.8),cex=2,pch=21)
points(nrkw_delHoms[,3],nrkw_delHoms[,2],bg=alpha( "#fb9a99",alpha=0.8),cex=2,pch=21)
points(offshore_delHoms[,3],offshore_delHoms[,2],bg=alpha("#33a02c",alpha=0.8),cex=2,pch=21)
legend(x=0,y=4250,xjust=FALSE,yjust=FALSE,bty="n",legend=c("Southern Residents","Northern Residents","Alaska Residents","Offshore","Transients"),
       pch=16,col= c("#a6cef9","#fb9a99","#1f78b4","#33a02c","#b2df8a"),cex=1.1)
legend(x=0,y=4250,xjust=FALSE,yjust=FALSE,bty="n",legend=c("Southern Residents","Northern Residents","Alaska Residents","Offshore","Transients"),
       pch=1,cex=1.1)

# regression models
# TKW
modDat <- tkw_delHoms
mod <- lm(modDat[,2]~modDat[,3])
int <- mod$coefficients[1]
slop <- mod$coefficients[2]
xRange <- c(min(modDat[,3],na.rm=TRUE),max(modDat[,3],na.rm=TRUE))
yRange <- c(int+slop*xRange[1],int+slop*xRange[2])
lines(x=xRange,y=yRange,col="orange",lwd=5)
lines(x=xRange,y=yRange,col="#b2df8a",lwd=3)
lines(x=xRange,y=yRange,col="black",lwd=3,lty="dashed")

# offshore
modDat <- offshore_delHoms
mod <- lm(modDat[,2]~modDat[,3])
# model is not statistically significant for the offshores

# ARKW
modDat <- arkw_delHoms
mod <- lm(modDat[,2]~modDat[,3])
int <- mod$coefficients[1]
slop <- mod$coefficients[2]
xRange <- c(min(modDat[,3],na.rm=TRUE),max(modDat[,3],na.rm=TRUE))
yRange <- c(int+slop*xRange[1],int+slop*xRange[2])
lines(x=xRange,y=yRange,col="orange",lwd=5)
lines(x=xRange,y=yRange,col="#1f78b4",lwd=3)
lines(x=xRange,y=yRange,col="black",lwd=3,lty="dashed")

# SRKW
modDat <- srkw_delHoms
mod <- lm(modDat[,2]~modDat[,3])
int <- mod$coefficients[1]
slop <- mod$coefficients[2]
xRange <- c(min(modDat[,3],na.rm=TRUE),max(modDat[,3],na.rm=TRUE))
yRange <- c(int+slop*xRange[1],int+slop*xRange[2])
lines(x=xRange,y=yRange,col="orange",lwd=5)
lines(x=xRange,y=yRange,col="#a6cef9",lwd=3)
lines(x=xRange,y=yRange,col="black",lwd=3,lty="dashed")


legend(x=0.25,y=70,xjust=FALSE,yjust=FALSE,bty="n",legend=c("Southern Residents","Northern Residents","Alaska Residents","Offshore","Transients"),
       pch=16,col= c("#a6cef9","#fb9a99","#1f78b4","#33a02c","#b2df8a"),cex=1.1)





#-----------------------------------------------------------------------------
##############################################################################
# deleterious homozygous count for moderate and missense mutations separately
##############################################################################
#-----------------------------------------------------------------------------

#######################################################################
# plot homozygous load versus Froh
#######################################################################
tkw_delHoms <- read.table("moderateHomozygousCount_tkw",header=TRUE)
arkw_delHoms <- read.table("moderateHomozygousCount_arkw",header=TRUE)
srkw_delHoms <- read.table("moderateHomozygousCount_srkw",header=TRUE)
offshore_delHoms <- read.table("moderateHomozygousCount_offshore",header=TRUE)
nrkw_delHoms <- read.table("moderateHomozygousCount_nrkw",header=TRUE)

arkw_delHoms <- cbind(arkw_delHoms,froh[match(arkw_delHoms[,1],froh[,1]),3:4])
srkw_delHoms <- cbind(srkw_delHoms,froh[match(srkw_delHoms[,1],froh[,1]),3:4])
tkw_delHoms <- cbind(tkw_delHoms,froh[match(tkw_delHoms[,1],froh[,1]),3:4])
nrkw_delHoms <- cbind(nrkw_delHoms,froh[match(nrkw_delHoms[,1],froh[,1]),3:4])
offshore_delHoms <- cbind(offshore_delHoms,froh[match(offshore_delHoms[,1],froh[,1]),3:4])

arkw_delHoms <- cbind(arkw_delHoms,het[match(arkw_delHoms[,1],het[,1]),2])
srkw_delHoms <- cbind(srkw_delHoms,het[match(srkw_delHoms[,1],het[,1]),2])
tkw_delHoms <- cbind(tkw_delHoms,het[match(tkw_delHoms[,1],het[,1]),2])
nrkw_delHoms <- cbind(nrkw_delHoms,het[match(nrkw_delHoms[,1],het[,1]),2])
offshore_delHoms <- cbind(offshore_delHoms,het[match(offshore_delHoms[,1],het[,1]),2])

colnames(arkw_delHoms)[3:5] <- c("Froh_1Mb","Froh_10Mb","het")
colnames(tkw_delHoms)[3:5] <- c("Froh_1Mb","Froh_10Mb","het")
colnames(srkw_delHoms)[3:5] <- c("Froh_1Mb","Froh_10Mb","het")
colnames(offshore_delHoms)[3:5] <- c("Froh_1Mb","Froh_10Mb","het")
colnames(nrkw_delHoms)[3:5] <- c("Froh_1Mb","Froh_10Mb","het")

allHoms <- rbind(tkw_delHoms,arkw_delHoms,srkw_delHoms,offshore_delHoms,nrkw_delHoms)
yLims <- c(min(as.numeric(allHoms[,2])),max(as.numeric(allHoms[,2])))
xLims <- c(min(froh[,3]),max(froh[,3]))

par(mfrow=c(1,1),mar=c(5,6,2,2))
#------------------------------------------------
# make the plot of homozygous load versus Froh
#------------------------------------------------
plot(c(0,0.5),c(2700,4500),type="n",ylab="Homozygous mutation load",
     xlab=expression(italic(""*F*"")["ROH,1Mb"]),cex.lab=1.8,axes=FALSE)
axis(side=1,at=seq(0,0.5,0.1))
axis(side=2,at=seq(2700,4500,200))

library(scales)
text(x=-0.12,y=4900,labels="B",cex=2)
points(srkw_delHoms[,3],srkw_delHoms[,2],bg=alpha("#a6cee3",alpha=0.8),cex=2,pch=21)
points(arkw_delHoms[,3],arkw_delHoms[,2],bg=alpha("#1f78b4",alpha=0.8),cex=2,pch=21)
points(tkw_delHoms[,3],tkw_delHoms[,2],bg=alpha("#b2df8a",alpha=0.8),cex=2,pch=21)
points(nrkw_delHoms[,3],nrkw_delHoms[,2],bg=alpha( "#fb9a99",alpha=0.8),cex=2,pch=21)
points(offshore_delHoms[,3],offshore_delHoms[,2],bg=alpha("#33a02c",alpha=0.8),cex=2,pch=21)
legend(x=0,y=4100,xjust=FALSE,yjust=FALSE,bty="n",legend=c("Southern Residents","Northern Residents","Alaska Residents","Offshore","Transients"),
       pch=16,col= c("#a6cef9","#fb9a99","#1f78b4","#33a02c","#b2df8a"),cex=1.1)
legend(x=0,y=4100,xjust=FALSE,yjust=FALSE,bty="n",legend=c("Southern Residents","Northern Residents","Alaska Residents","Offshore","Transients"),
       pch=1,cex=1.1)

# regression models
# TKW
modDat <- tkw_delHoms
mod <- lm(modDat[,2]~modDat[,3])
int <- mod$coefficients[1]
slop <- mod$coefficients[2]
xRange <- c(min(modDat[,3],na.rm=TRUE),max(modDat[,3],na.rm=TRUE))
yRange <- c(int+slop*xRange[1],int+slop*xRange[2])
lines(x=xRange,y=yRange,col="orange",lwd=5)
lines(x=xRange,y=yRange,col="#b2df8a",lwd=3)
lines(x=xRange,y=yRange,col="black",lwd=3,lty="dashed")

# offshore
modDat <- offshore_delHoms
mod <- lm(modDat[,2]~modDat[,3])
# model is not statistically significant for the offshores

# ARKW
modDat <- arkw_delHoms
mod <- lm(modDat[,2]~modDat[,3])
int <- mod$coefficients[1]
slop <- mod$coefficients[2]
xRange <- c(min(modDat[,3],na.rm=TRUE),max(modDat[,3],na.rm=TRUE))
yRange <- c(int+slop*xRange[1],int+slop*xRange[2])
lines(x=xRange,y=yRange,col="orange",lwd=5)
lines(x=xRange,y=yRange,col="#1f78b4",lwd=3)
lines(x=xRange,y=yRange,col="black",lwd=3,lty="dashed")

# SRKW
modDat <- srkw_delHoms
mod <- lm(modDat[,2]~modDat[,3])
int <- mod$coefficients[1]
slop <- mod$coefficients[2]
xRange <- c(min(modDat[,3],na.rm=TRUE),max(modDat[,3],na.rm=TRUE))
yRange <- c(int+slop*xRange[1],int+slop*xRange[2])
lines(x=xRange,y=yRange,col="orange",lwd=5)
lines(x=xRange,y=yRange,col="#a6cef9",lwd=3)
lines(x=xRange,y=yRange,col="black",lwd=3,lty="dashed")


#-----------------------------------------------------------------------
########################################################################
# deleterious homozygous count for moderate and LOF mutations separately
########################################################################
#-----------------------------------------------------------------------

#######################################################################
# plot homozygous load versus Froh
#######################################################################
tkw_delHoms <- read.table("lofHomozygousCount_tkw",header=TRUE)
arkw_delHoms <- read.table("lofHomozygousCount_arkw",header=TRUE)
srkw_delHoms <- read.table("lofHomozygousCount_srkw",header=TRUE)
offshore_delHoms <- read.table("lofHomozygousCount_offshore",header=TRUE)
nrkw_delHoms <- read.table("lofHomozygousCount_nrkw",header=TRUE)

arkw_delHoms <- cbind(arkw_delHoms,froh[match(arkw_delHoms[,1],froh[,1]),3:4])
srkw_delHoms <- cbind(srkw_delHoms,froh[match(srkw_delHoms[,1],froh[,1]),3:4])
tkw_delHoms <- cbind(tkw_delHoms,froh[match(tkw_delHoms[,1],froh[,1]),3:4])
nrkw_delHoms <- cbind(nrkw_delHoms,froh[match(nrkw_delHoms[,1],froh[,1]),3:4])
offshore_delHoms <- cbind(offshore_delHoms,froh[match(offshore_delHoms[,1],froh[,1]),3:4])

arkw_delHoms <- cbind(arkw_delHoms,het[match(arkw_delHoms[,1],het[,1]),2])
srkw_delHoms <- cbind(srkw_delHoms,het[match(srkw_delHoms[,1],het[,1]),2])
tkw_delHoms <- cbind(tkw_delHoms,het[match(tkw_delHoms[,1],het[,1]),2])
nrkw_delHoms <- cbind(nrkw_delHoms,het[match(nrkw_delHoms[,1],het[,1]),2])
offshore_delHoms <- cbind(offshore_delHoms,het[match(offshore_delHoms[,1],het[,1]),2])

colnames(arkw_delHoms)[3:5] <- c("Froh_1Mb","Froh_10Mb","het")
colnames(tkw_delHoms)[3:5] <- c("Froh_1Mb","Froh_10Mb","het")
colnames(srkw_delHoms)[3:5] <- c("Froh_1Mb","Froh_10Mb","het")
colnames(offshore_delHoms)[3:5] <- c("Froh_1Mb","Froh_10Mb","het")
colnames(nrkw_delHoms)[3:5] <- c("Froh_1Mb","Froh_10Mb","het")

allHoms <- rbind(tkw_delHoms,arkw_delHoms,srkw_delHoms,offshore_delHoms,nrkw_delHoms)
yLims <- c(min(as.numeric(allHoms[,2])),max(as.numeric(allHoms[,2])))
xLims <- c(min(froh[,3]),max(froh[,3]))

par(mfrow=c(1,1),mar=c(5,6,2,2))
#------------------------------------------------
# make the plot of homozygous load versus Froh
#------------------------------------------------
plot(c(0,0.5),c(70,170),type="n",ylab="Homozygous mutation load",
     xlab=expression(italic(""*F*"")["ROH,1Mb"]),cex.lab=1.8,axes=FALSE)
axis(side=1,at=seq(0,0.5,0.1))
axis(side=2,at=seq(70,170,20))

library(scales)
text(x=-0.12,y=4900,labels="B",cex=2)
points(srkw_delHoms[,3],srkw_delHoms[,2],bg=alpha("#a6cee3",alpha=0.8),cex=2,pch=21)
points(arkw_delHoms[,3],arkw_delHoms[,2],bg=alpha("#1f78b4",alpha=0.8),cex=2,pch=21)
points(tkw_delHoms[,3],tkw_delHoms[,2],bg=alpha("#b2df8a",alpha=0.8),cex=2,pch=21)
points(nrkw_delHoms[,3],nrkw_delHoms[,2],bg=alpha( "#fb9a99",alpha=0.8),cex=2,pch=21)
points(offshore_delHoms[,3],offshore_delHoms[,2],bg=alpha("#33a02c",alpha=0.8),cex=2,pch=21)
legend(x=0,y=4250,xjust=FALSE,yjust=FALSE,bty="n",legend=c("Southern Residents","Northern Residents","Alaska Residents","Offshore","Transients"),
       pch=16,col= c("#a6cef9","#fb9a99","#1f78b4","#33a02c","#b2df8a"),cex=1.1)
legend(x=0,y=4250,xjust=FALSE,yjust=FALSE,bty="n",legend=c("Southern Residents","Northern Residents","Alaska Residents","Offshore","Transients"),
       pch=1,cex=1.1)

# regression models
# TKW
modDat <- tkw_delHoms
mod <- lm(modDat[,2]~modDat[,3])
int <- mod$coefficients[1]
slop <- mod$coefficients[2]
xRange <- c(min(modDat[,3],na.rm=TRUE),max(modDat[,3],na.rm=TRUE))
yRange <- c(int+slop*xRange[1],int+slop*xRange[2])
lines(x=xRange,y=yRange,col="orange",lwd=5)
lines(x=xRange,y=yRange,col="#b2df8a",lwd=3)
lines(x=xRange,y=yRange,col="black",lwd=3,lty="dashed")

# offshore
modDat <- offshore_delHoms
mod <- lm(modDat[,2]~modDat[,3])
# model is not statistically significant for the offshores

# ARKW
modDat <- arkw_delHoms
mod <- lm(modDat[,2]~modDat[,3])
int <- mod$coefficients[1]
slop <- mod$coefficients[2]
xRange <- c(min(modDat[,3],na.rm=TRUE),max(modDat[,3],na.rm=TRUE))
yRange <- c(int+slop*xRange[1],int+slop*xRange[2])
lines(x=xRange,y=yRange,col="orange",lwd=5)
lines(x=xRange,y=yRange,col="#1f78b4",lwd=3)
lines(x=xRange,y=yRange,col="black",lwd=3,lty="dashed")

# SRKW
modDat <- srkw_delHoms
mod <- lm(modDat[,2]~modDat[,3])
int <- mod$coefficients[1]
slop <- mod$coefficients[2]
xRange <- c(min(modDat[,3],na.rm=TRUE),max(modDat[,3],na.rm=TRUE))
yRange <- c(int+slop*xRange[1],int+slop*xRange[2])
lines(x=xRange,y=yRange,col="orange",lwd=5)
lines(x=xRange,y=yRange,col="#a6cef9",lwd=3)
lines(x=xRange,y=yRange,col="black",lwd=3,lty="dashed")


legend(x=0.25,y=70,xjust=FALSE,yjust=FALSE,bty="n",legend=c("Southern Residents","Northern Residents","Alaska Residents","Offshore","Transients"),
       pch=16,col= c("#a6cef9","#fb9a99","#1f78b4","#33a02c","#b2df8a"),cex=1.1)






#######################################################################
# plot homozygous load versus Hsnp
#######################################################################
tkw_delHoms <- read.table("moderateHomozygousCount_tkw",header=TRUE)
arkw_delHoms <- read.table("moderateHomozygousCount_arkw",header=TRUE)
srkw_delHoms <- read.table("moderateHomozygousCount_srkw",header=TRUE)
offshore_delHoms <- read.table("moderateHomozygousCount_offshore",header=TRUE)
nrkw_delHoms <- read.table("moderateHomozygousCount_nrkw",header=TRUE)

arkw_delHoms <- cbind(arkw_delHoms,froh[match(arkw_delHoms[,1],froh[,1]),3:4])
srkw_delHoms <- cbind(srkw_delHoms,froh[match(srkw_delHoms[,1],froh[,1]),3:4])
tkw_delHoms <- cbind(tkw_delHoms,froh[match(tkw_delHoms[,1],froh[,1]),3:4])
nrkw_delHoms <- cbind(nrkw_delHoms,froh[match(nrkw_delHoms[,1],froh[,1]),3:4])
offshore_delHoms <- cbind(offshore_delHoms,froh[match(offshore_delHoms[,1],froh[,1]),3:4])

arkw_delHoms <- cbind(arkw_delHoms,het[match(arkw_delHoms[,1],het[,1]),2])
srkw_delHoms <- cbind(srkw_delHoms,het[match(srkw_delHoms[,1],het[,1]),2])
tkw_delHoms <- cbind(tkw_delHoms,het[match(tkw_delHoms[,1],het[,1]),2])
nrkw_delHoms <- cbind(nrkw_delHoms,het[match(nrkw_delHoms[,1],het[,1]),2])
offshore_delHoms <- cbind(offshore_delHoms,het[match(offshore_delHoms[,1],het[,1]),2])

colnames(arkw_delHoms)[3:5] <- c("Froh_1Mb","Froh_10Mb","het")
colnames(tkw_delHoms)[3:5] <- c("Froh_1Mb","Froh_10Mb","het")
colnames(srkw_delHoms)[3:5] <- c("Froh_1Mb","Froh_10Mb","het")
colnames(offshore_delHoms)[3:5] <- c("Froh_1Mb","Froh_10Mb","het")
colnames(nrkw_delHoms)[3:5] <- c("Froh_1Mb","Froh_10Mb","het")

allHoms <- rbind(tkw_delHoms,arkw_delHoms,srkw_delHoms,offshore_delHoms,nrkw_delHoms)
yLims <- c(min(as.numeric(allHoms[,2])),max(as.numeric(allHoms[,2])))
xLims <- c(0.05,0.25)

par(mfrow=c(1,1),mar=c(5,6,2,2))
#------------------------------------------------
# make the plot of homozygous load versus Froh
#------------------------------------------------
plot(xLims,c(2700,4600),type="n",ylab="Homozygous mutation load",
     xlab=expression(italic(""*H*"")["SNP"]),cex.lab=1.8,axes=FALSE)
axis(side=1,at=seq(0.05,0.25,0.05))
axis(side=2,at=seq(2700,4700,400))

library(scales)
text(x=-0.12,y=4900,labels="B",cex=2)
points(srkw_delHoms[,5],srkw_delHoms[,2],bg=alpha("#a6cee3",alpha=0.8),cex=2,pch=21)
points(arkw_delHoms[,5],arkw_delHoms[,2],bg=alpha("#1f78b4",alpha=0.8),cex=2,pch=21)
points(tkw_delHoms[,5],tkw_delHoms[,2],bg=alpha("#b2df8a",alpha=0.8),cex=2,pch=21)
points(nrkw_delHoms[,5],nrkw_delHoms[,2],bg=alpha( "#fb9a99",alpha=0.8),cex=2,pch=21)
points(offshore_delHoms[,5],offshore_delHoms[,2],bg=alpha("#33a02c",alpha=0.8),cex=2,pch=21)
legend(x=0.05,y=2700,xjust=FALSE,yjust=FALSE,bty="n",legend=c("Southern Residents","Northern Residents","Alaska Residents","Offshore","Transients"),
       pch=16,col= c("#a6cef9","#fb9a99","#1f78b4","#33a02c","#b2df8a"),cex=1.1)
legend(x=0.05,y=2700,xjust=FALSE,yjust=FALSE,bty="n",legend=c("Southern Residents","Northern Residents","Alaska Residents","Offshore","Transients"),
       pch=1,cex=1.1)

# regression models
# TKW
modDat <- tkw_delHoms
mod <- lm(modDat[,2]~modDat[,5])
int <- mod$coefficients[1]
slop <- mod$coefficients[2]
xRange <- c(min(modDat[,5],na.rm=TRUE),max(modDat[,5],na.rm=TRUE))
yRange <- c(int+slop*xRange[1],int+slop*xRange[2])
lines(x=xRange,y=yRange,col="orange",lwd=5)
lines(x=xRange,y=yRange,col="#b2df8a",lwd=3)
lines(x=xRange,y=yRange,col="black",lwd=3,lty="dashed")

# offshore
modDat <- offshore_delHoms
mod <- lm(modDat[,2]~modDat[,3])
# model is not statistically significant for the offshores

# ARKW
modDat <- arkw_delHoms
mod <- lm(modDat[,2]~modDat[,5])
int <- mod$coefficients[1]
slop <- mod$coefficients[2]
xRange <- c(min(modDat[,5],na.rm=TRUE),max(modDat[,5],na.rm=TRUE))
yRange <- c(int+slop*xRange[1],int+slop*xRange[2])
lines(x=xRange,y=yRange,col="orange",lwd=5)
lines(x=xRange,y=yRange,col="#1f78b4",lwd=3)
lines(x=xRange,y=yRange,col="black",lwd=3,lty="dashed")

# SRKW
modDat <- srkw_delHoms
mod <- lm(modDat[,2]~modDat[,5])
int <- mod$coefficients[1]
slop <- mod$coefficients[2]
xRange <- c(min(modDat[,5],na.rm=TRUE),max(modDat[,5],na.rm=TRUE))
yRange <- c(int+slop*xRange[1],int+slop*xRange[2])
lines(x=xRange,y=yRange,col="orange",lwd=5)
lines(x=xRange,y=yRange,col="#a6cef9",lwd=3)
lines(x=xRange,y=yRange,col="black",lwd=3,lty="dashed")
