################################################################################
# plot population distributions of Froh first
################################################################################

setwd("~/Documents/orca/analyses_31May2022/ROH3")
Froh <- read.table("Froh_7August2022",header=TRUE)

setwd("/Users/martin.kardos/Documents/orca/orca")


popKey <- read.table("popKey.txt",header=TRUE)
popVec <- popKey[match(Froh[,1],popKey[,2]),1]
froh <- Froh

pops <- rev(c("SRKW", "NRKW" ,"ARKW" ,"offshore",  "TKW"))
# make a barplot of inbreeding for all populations

popMeans <- NULL
for(i in 1:length(pops)){
  theseMeans <- colMeans(froh[froh[,2] == pops[i],3:4])
  popMeans <- rbind(popMeans, theseMeans)
}
popMeans <- cbind(pops,popMeans)
#############################################
#bar plots only for 1Mb and 10Mb ROH
#############################################
theseFrohs <- froh[,c(1,2,3)]

theseFrohs <- theseFrohs[order(theseFrohs[,3]),]

colVec <- rep(NA,nrow(froh))
colVec[theseFrohs[,2] == "ARKW"] <-  "#1f78b4"
colVec[theseFrohs[,2] == "SRKW"] <-  "#a6cef9"
colVec[theseFrohs[,2] == "TKW"] <-  "#b2df8a"
colVec[theseFrohs[,2] == "offshore"] <-  "#33a02c"
colVec[theseFrohs[,2] == "NRKW"] <-  "#fb9a99"

meanCols <- rev(c("#a6cef9","#fb9a99","#1f78b4","#33a02c", "#b2df8a"))

par(mfrow=c(2,4),mar=c(5,6,2,4),xpd=TRUE)

# skip two to make room for the map
plot(c(0,1),c(0,1), axes=FALSE,xlab="",ylab="",type="n")
plot(c(0,1),c(0,1), axes=FALSE,xlab="",ylab="",type="n")

#############################################
#bar plot with minimum ROH length of 1 mb
#############################################

plot(c(0,nrow(theseFrohs)),c(0,0.5),ylab=expression(paste(italic(""*F*"")["ROH,1Mb"],sep="")),type="n",axes=FALSE,xlab="",cex.lab=1.5)
axis(side=2,at=seq(0,0.5,0.1))

for(i in 1:nrow(theseFrohs)){
  rect(xleft=i-1,xright=i,ybottom=0,ytop=theseFrohs[i,3],col=colVec[i],border=colVec[i])
}


########## add CIs
FrohMat <- matrix(NA,ncol=length(pops),nrow=1000)
for(i in 1:1000){
  fVec <- rep(NA,5)
  for(j in 1:5){
    fVec[j] <- mean(sample(froh[froh[,2] == pops[j],3],sum(froh[,2] == pops[j]),replace=TRUE))
  }
  FrohMat[i,] <- fVec
}

FrohEsts <- colMeans(FrohMat)
CIMat <- matrix(NA,ncol=2,nrow=5)

for(i in 1:5){
  CIMat [i,] <- quantile(FrohMat[,i],probs=c(0.025,0.975))
}

plotX <- c(160,162,164,164,166)
for(i in c(1,2,3,5)){
  points(x=plotX[i],y=FrohEsts[i],cex=1.5,col=meanCols[i],pch=16)
  points(x=plotX[i],y=FrohEsts[i],cex=1.54,col=meanCols[i],pch=16)
  arrows(x0=plotX[i],x1=plotX[i],y0=FrohEsts[i],y1=CIMat[i,1],angle=90,col=meanCols[i],length=0.02)
  arrows(x0=plotX[i],x1=plotX[i],y0=FrohEsts[i],y1=CIMat[i,2],angle=90,col=meanCols[i],length=0.02)
}
text(x=75,y=-0.05,labels="147 Individuals",cex=1.3)
text(x=-40,y=0.545,labels="B",cex=2)


#############################################
#bar plot with minimum ROH length of 10 mb
#############################################

theseFrohs <- froh[,c(1,2,4)]

theseFrohs <- theseFrohs[order(theseFrohs[,3]),]

colVec <- rep(NA,nrow(froh))
colVec[theseFrohs[,2] == "ARKW"] <-  "#1f78b4"
colVec[theseFrohs[,2] == "SRKW"] <-  "#a6cef9"
colVec[theseFrohs[,2] == "TKW"] <-  "#b2df8a"
colVec[theseFrohs[,2] == "offshore"] <-  "#33a02c"
colVec[theseFrohs[,2] == "NRKW"] <-  "#fb9a99"

plot(c(0,nrow(theseFrohs)),c(0,0.15),ylab=expression(paste(italic(""*F*"")["ROH,10Mb"],sep="")),type="n",axes=FALSE,xlab="",cex.lab=1.5)
axis(side=2,at=seq(0,0.15,0.03))

for(i in 1:nrow(theseFrohs)){
  rect(xleft=i-1,xright=i,ybottom=0,ytop=theseFrohs[i,3],col=colVec[i],border=colVec[i])
}


text(x=75,y=-0.015,labels="147 Individuals",cex=1.3)
text(x=-40,y=0.1625,labels="C",cex=2)




########## add CIs
FrohMat <- matrix(NA,ncol=length(pops),nrow=1000)
for(i in 1:1000){
  fVec <- rep(NA,5)
  for(j in 1:5){
    fVec[j] <- mean(sample(froh[froh[,2] == pops[j],4],sum(froh[,2] == pops[j]),replace=TRUE))
  }
  FrohMat[i,] <- fVec
}

FrohEsts <- colMeans(FrohMat)
CIMat <- matrix(NA,ncol=2,nrow=5)

for(i in 1:5){
  CIMat [i,] <- quantile(FrohMat[,i],probs=c(0.025,0.975))
}


for(i in c(1,2,3,5)){
  points(x=plotX[i],y=FrohEsts[i],cex=1.5,col=meanCols[i],pch=16)
  points(x=plotX[i],y=FrohEsts[i],cex=1.54,col=meanCols[i],pch=16)
  arrows(x0=plotX[i],x1=plotX[i],y0=FrohEsts[i],y1=CIMat[i,1],angle=90,col=meanCols[i],length=0.02)
  arrows(x0=plotX[i],x1=plotX[i],y0=FrohEsts[i],y1=CIMat[i,2],angle=90,col=meanCols[i],length=0.02)
}



legend(x=0,y=0.07,xjust=FALSE,yjust=FALSE,bty="n",legend=c("Southern Residents","Northern Residents","Alaska Residents","Offshore","Transients"),
       pch=15,col= c("#a6cef9","#fb9a99","#1f78b4","#33a02c","#b2df8a"),cex=1.2)


# skip two to make room for map and NJ tree
plot(c(0,1),c(0,1), axes=FALSE,xlab="",ylab="",type="n")
plot(c(0,1),c(0,1), axes=FALSE,xlab="",ylab="",type="n")

#--------------------------------------------------------
#########################################################
# plot Ne
#########################################################
#--------------------------------------------------------

# plot the different GONE reps for the SRKWs
par(xpd=FALSE)
library(scales)
library(matrixStats)



plot(c(0,150),c(0,8000),type="n",xlab="Generations back in time",ylab=expression(paste("Hissrotical ",italic(""*N*"")[e],sep="")),
     cex.lab=1.5)

##################################################################################
# Alaska residents
##################################################################################
setwd( "/Users/martin.kardos/Documents/orca/analyses_31May2022/Ne/GONE/lowMissingIndivs/arkw")
files <- paste("outfileLD_",1:500,"_GONE_Nebest",sep="")
NeMat <- NULL
for(i in 1:500){
  dat <- read.table(files[i],skip=2)
  NeMat <- cbind(NeMat,dat[,2])
}

# get CI
NeCI <- matrix(NA,nrow=500,ncol=2)
NeCI2 <- matrix(NA,nrow=500,ncol=2)

for(i in 1:500){
  NeCI[i,] <- quantile(NeMat[i,],probs=c(0.025,0.975))
  NeCI2[i,] <- quantile(NeMat[i,],probs=c(0.25,0.75))
}


lines(1:500,rowMedians(NeMat[1:500,]),col="#1f78b4",lwd=4)
polygon(x=c(1:500,rev(1:500)),y=c(NeCI[1:500,1],rev(NeCI[1:500,2])),col=adjustcolor("#1f78b4",alpha.f=0.2),border=NA)
polygon(x=c(1:500,rev(1:500)),y=c(NeCI2[1:500,1],rev(NeCI2[1:500,2])),col=adjustcolor("#1f78b4",alpha.f=0.3),border=NA)


##################################################################################
# Transients
##################################################################################
setwd( "/Users/martin.kardos/Documents/orca/analyses_31May2022/Ne/GONE/lowMissingIndivs/tkw")

files <- paste("outfileLD_",1:500,"_GONE_Nebest",sep="")
NeMat <- NULL
for(i in 1:500){
  dat <- read.table(files[i],skip=2)
  NeMat <- cbind(NeMat,dat[,2])
}


# get CI
# get CI
NeCI <- matrix(NA,nrow=500,ncol=2)
NeCI2 <- matrix(NA,nrow=500,ncol=2)

for(i in 1:500){
  NeCI[i,] <- quantile(NeMat[i,],probs=c(0.025,0.975))
  NeCI2[i,] <- quantile(NeMat[i,],probs=c(0.25,0.75))
}
lines(1:500,rowMedians(NeMat[1:500,]),col="#b2df8a",lwd=4)
polygon(x=c(1:500,rev(1:500)),y=c(NeCI[1:500,1],rev(NeCI[1:500,2])),col=adjustcolor("#b2df8a",alpha.f=0.4),border=NA,angle=135,density=30)
polygon(x=c(1:500,rev(1:500)),y=c(NeCI2[1:500,1],rev(NeCI2[1:500,2])),col=adjustcolor("#b2df8a",alpha.f=0.3),border=NA)


polygon(x=c(1:500,rev(1:500)),y=c(NeCI[1:500,1],rev(NeCI[1:500,2])),col=adjustcolor("#b2df8a",alpha.f=0.2),border=NA)
polygon(x=c(1:500,rev(1:500)),y=c(NeCI2[1:500,1],rev(NeCI2[1:500,2])),col=adjustcolor("#b2df8a",alpha.f=0.3),border=NA)





##################################################################################
# southern residents
##################################################################################
setwd( "/Users/martin.kardos/Documents/orca/analyses_31May2022/Ne/GONE/lowMissingIndivs/srkw")
files <- paste("outfileLD_",1:500,"_GONE_Nebest",sep="")
NeMat <- NULL
for(i in 1:500){
  dat <- read.table(files[i],skip=2)
  NeMat <- cbind(NeMat,dat[,2])
}

# get CI
# get CI
NeCI <- matrix(NA,nrow=500,ncol=2)
NeCI2 <- matrix(NA,nrow=500,ncol=2)

for(i in 1:500){
  NeCI[i,] <- quantile(NeMat[i,],probs=c(0.025,0.975))
  NeCI2[i,] <- quantile(NeMat[i,],probs=c(0.25,0.75))
}
lines(1:500,rowMedians(NeMat[1:500,]),col="#a6cef9",lwd=4)

polygon(x=c(1:500,rev(1:500)),y=c(NeCI[1:500,1],rev(NeCI[1:500,2])),col=adjustcolor("#a6cef9",alpha.f=0.4),border=NA,angle=45,density=30)
polygon(x=c(1:500,rev(1:500)),y=c(NeCI2[1:500,1],rev(NeCI2[1:500,2])),col=adjustcolor("#a6cef9",alpha.f=0.3),border=NA)

polygon(x=c(1:500,rev(1:500)),y=c(NeCI[1:500,1],rev(NeCI[1:500,2])),col=adjustcolor("#a6cef9",alpha.f=0.2),border=NA)
polygon(x=c(1:500,rev(1:500)),y=c(NeCI2[1:500,1],rev(NeCI2[1:500,2])),col=adjustcolor("#a6cef9",alpha.f=0.3),border=NA)
par(xpd=TRUE)
text(x=-40,y=8800,labels="D",cex=2)

sr_NeVec <- c(27.4, 27.4,rowMedians(NeMat))
srNe <- NULL
srNe <- cbind(srNe,sr_NeVec)
write.table(srNe,file="SRKW_pastNe",quote=FALSE,row.names=FALSE)
srHistNe <- NULL
srHistNe <- cbind(1:nrow(NeMat),rowMedians(NeMat))
write.table(srHistNe,file="SRKWHistoricalNe",quote=FALSE,row.names=FALSE)
###############################################################################
# plot LDNe CIs
###############################################################################

SRNe <- c(27.4,22.7,33.3)
ARNe <- c(38.9,22.1,107.4)
TRNe <- c(86.0,54.8,188.4)

NeMat <- rbind(SRNe,ARNe,TRNe)
plot(c(0.5,1.5,2.5),NeMat[,1],xlim=c(0,3),ylim=c(0,250),col=c("#a6cef9","#1f78b4","#b2df8a"),pch=16,cex=1.6,xlab="Population",cex.lab=1.3,
     ylab=expression(paste("Contemporary ",italic(""*N*"")[e])),axes=FALSE)
text(x=-0.8,y=275,labels="E",cex=2)

arrows(x0=0.5,x1=0.5,y0=NeMat[1,1],y1=NeMat[1,2],col="#a6cef9",angle=90,length=0.1)
arrows(x0=0.5,x1=0.5,y0=NeMat[1,1],y1=NeMat[1,3],col="#a6cef9",angle=90,length=0.1)

arrows(x0=1.5,x1=1.5,y0=NeMat[2,1],y1=NeMat[2,2],col="#1f78b4",angle=90,length=0.1)
arrows(x0=1.5,x1=1.5,y0=NeMat[2,1],y1=NeMat[2,3],col="#1f78b4",angle=90,length=0.1)

arrows(x0=2.5,x1=2.5,y0=NeMat[3,1],y1=NeMat[3,2],col="#b2df8a",angle=90,length=0.1)
arrows(x0=2.5,x1=2.5,y0=NeMat[3,1],y1=NeMat[3,3],col="#b2df8a",angle=90,length=0.1)

axis(side=2,at=seq(0,250,50))
axis(side=1,at=c(0,1,2,3),labels=rep("",4))
text(x=0.5,y=-35,labels="SRKW")
text(x=1.5,y=-35,labels="ARKW")
text(x=2.5,y=-35,labels="TKW")