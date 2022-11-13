
#######################################
# plot very recent Ne
#######################################
setwd("/Users/martin.kardos/Documents/orca/orca")

# plot the different GONE reps for the SRKWs
par(xpd=FALSE,mfrow=c(1,1))
library(scales)
library(matrixStats)



plot(c(0,50),c(0,1500),type="n",xlab="Generations back in time",ylab=expression(paste("Hissrotical ",italic(""*N*"")[e],sep="")),
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
text(x=-40,y=5500,labels="D",cex=2)
