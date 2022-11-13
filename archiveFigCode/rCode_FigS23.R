# make a multi-panel figure showing population size through time for different mutation models
setwd("~/Documents/orca/analyses_31May2022/Projections")


shapes <- rep(0.2,14)
scales <- c(rep(0.2,6),rep(0.1,6),0.05,0.2)
propLethals <- c(rep(c(0,0.01,0.02),4),0,0)
betas <- c(50,50,50,13,13,13,50,50,50,13,13,13,50,13)
betas2 <- as.character(betas)

infiles <- c(paste("popSize_shape",shapes,"_scale",scales,"_propLethal",propLethals,"_beta",betas,sep="")[1:13],
             "popSize_shape0.2_scale0.2_propLethal0_beta13_noMutationInFuture",
             "popSize_BforYearlySurvival0.12")

par(mfrow=c(4,4),mar=c(4,4,2,1))
noInbDepDat <- read.table("popSize_noInbDep",header=TRUE)[1:200,1:100]
noInbDepDat <- cbind(rep(74,nrow(noInbDepDat)),noInbDepDat)

for(i in 1:length(infiles)){
  popSize <- read.table(infiles[i],header=TRUE)
  # plot results
  library(scales)
  #newPopSize <- popSize[which(popSize[,10] >0),]                 # get rid of the extra rows where this are not data
  newPopSize <- popSize[1:200,1:100]                 # get rid of the extra rows where this are not data
  
  newPopSize <- cbind(rep(74,nrow(newPopSize)),newPopSize)
  
  if(i %in% c(1,5,9)) plot(c(0,100),c(0,150),xlab="",ylab="Population Size",type="n",cex.main=1,main=paste("Model ",i,sep=""))
  if(i %in% c(2:4,6:8,10:12)) plot(c(0,100),c(0,150),xlab="",ylab="",type="n",cex.main=1,main=paste("Model ",i,sep=""))
  if(i == 13) plot(c(0,100),c(0,150),xlab="Year",ylab="Population Size",type="n",cex.main=1,main=paste("Model ",i,sep=""))
  if(i > 13)plot(c(0,100),c(0,150),xlab="Year",ylab="",type="n",cex.main=1,main=paste("Model ",i,sep=""))
  
  #---------------------------------------------------------------------
  # plot population sizes for simulations without inbreeding depression
  #---------------------------------------------------------------------
  for(j in 1:nrow(noInbDepDat)){
    lines(0:100,as.numeric(noInbDepDat[j,]),col=alpha("darkgray",alpha=0.02))
  }
  medians <- apply(noInbDepDat,2,median)
  quant95  <- matrix(NA,nrow=2,ncol=ncol(noInbDepDat))
  quant50  <- matrix(NA,nrow=2,ncol=ncol(noInbDepDat))
  
  for(j in 1:ncol(noInbDepDat)){
    quant95[,j] <- quantile(noInbDepDat[,j],probs=c(0.025,0.975))
    quant50[,j] <- quantile(noInbDepDat[,j],probs=c(0.25,0.75))
  }
  polygon(c(0:100,rev(0:100)),c(quant95[1,1:101],rev(quant95[2,1:101])),col=alpha("darkgray",alpha=0.15),border=NA)
  polygon(c(0:100,rev(0:100)),c(quant50[1,1:101],rev(quant50[2,1:101])),col=alpha("darkgray",alpha=0.35),border=NA)
  lines(0:100,medians[1:101],col="darkgray",lwd=3)
  
  #---------------------------------------------------------------------
  # plot population sizes for simulations with inbreeding depression
  #---------------------------------------------------------------------
  for(j in 1:nrow(newPopSize)){
    lines(0:100,as.numeric(newPopSize[j,1:101]),col=alpha("darkred",alpha=0.02))
  }

  medians <- apply(newPopSize,2,median)
  quant95  <- matrix(NA,nrow=2,ncol=ncol(newPopSize))
  quant50  <- matrix(NA,nrow=2,ncol=ncol(newPopSize))
  
  for(j in 1:ncol(newPopSize)){
    quant95[,j] <- quantile(newPopSize[,j],probs=c(0.025,0.975))
    quant50[,j] <- quantile(newPopSize[,j],probs=c(0.25,0.75))
  }
  polygon(c(0:100,rev(0:100)),c(quant95[1,1:101],rev(quant95[2,1:101])),col=alpha("darkred",alpha=0.1),border=NA)
  polygon(c(0:100,rev(0:100)),c(quant50[1,1:101],rev(quant50[2,1:101])),col=alpha("darkred",alpha=0.3),border=NA)
  lines(0:100,medians[1:101],col="darkred",lwd=3)
  
  if(i < 15){
    legend(x=-18,y=120,xjust=FALSE,yjust=FALSE,bty="n",legend = c(paste("Shape=",shapes[i],sep=""),
   paste("Scale=",scales[i],sep="")))
    
    legend(x=30,y=134,xjust=FALSE,yjust=FALSE,bty="n",legend =paste("P(lethal)=",propLethals[i],sep=""))
    legend(x=30,y=120,xjust=FALSE,yjust=FALSE,bty="n",legend = expression(beta))
    legend(x=36,y=120,xjust=FALSE,yjust=FALSE,bty="n",legend = paste("=",betas2[i],sep=""))
    
  }
  if(i == 15){
    legend(x=-10,y=130,xjust=FALSE,yjust=FALSE,bty="n",cex=0.8,
           legend = expression(paste(italic(""*S*""), " = ",italic(""*S*"")[0],e^-italic(""*BF*""))))
  }
}
