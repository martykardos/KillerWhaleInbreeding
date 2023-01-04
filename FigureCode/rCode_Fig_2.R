setwd("/Users/martin.kardos/Documents/orca/analyses_31May2022/demograhpy")

fig_s18 <- as.data.frame(readRDS("fig_s18.rds"))
male1Mb <- fig_s18[which(fig_s18$name == "Male 1Mb"),]
male1Mb <- male1Mb[order(male1Mb$draw),]


male10Mb <- fig_s18[which(fig_s18$name == "Male 10Mb"),]
male10Mb <- male10Mb[order(male10Mb$draw),]

female1Mb <- fig_s18[which(fig_s18$name == "Female 1Mb"),]
female1Mb <- female1Mb[order(female1Mb$draw),]

female10Mb <- fig_s18[which(fig_s18$name == "Female 10Mb"),]
female10Mb <- female10Mb[order(female10Mb$draw),]

xLims <- c(min(fig_s18$Froh_raw),max(fig_s18$Froh_raw))
yLims <- c(0,1)

library(scales)
par(mfrow=c(2,1),mar=c(5,5,2,1),xpd=FALSE)

# female survival to 40 years versus FROH10Mb
plot(female10Mb$Froh_raw,female10Mb$pred_invlogit,ylim=yLims,xlab=expression(italic(""*F*"")["ROH,10Mb"]),
     ylab="Female survival to 40 years",type="n",cex.lab=1.2)
par(xpd=TRUE)
text(x=-0.056,y=1.1,labels="A",cex=1.5)
par(xpd=FALSE)
frohs <- sort(unique(female10Mb$Froh_raw))
means <- rep(NA,length(frohs))
for(i in 1:length(means)){
  means[i] <- median(female10Mb$pred_invlogit[which(female10Mb$Froh_raw == frohs[i])])
}

survPreds <- matrix(NA,nrow=5000,ncol=length(unique(frohs)))
for(i in 1:5000){
  thisRep <- female10Mb[which(female10Mb$draw == i),]
  survPreds[i,] <- thisRep$pred_invlogit
  lines(thisRep$Froh_raw,thisRep$pred_invlogit,col=alpha("blue",alpha=0.008),lwd=0.5)
}

quantSurv95  <- matrix(NA,nrow=2,ncol=length(frohs))
quantSurv50  <- matrix(NA,nrow=2,ncol=length(frohs))

for(i in 1:length(frohs)){
  quantSurv95[,i] <- quantile(survPreds[,i],probs=c(0.025,0.975))
  quantSurv50[,i] <- quantile(survPreds[,i],probs=c(0.25,0.75))
}
polygon(c(frohs,rev(frohs)),c(quantSurv95[1,],rev(quantSurv95[2,])),col=alpha("blue3",alpha=0.05),border=NA)
polygon(c(frohs,rev(frohs)),c(quantSurv50[1,],rev(quantSurv50[2,])),col=alpha("blue3",alpha=0.2),border=NA)
lines(frohs,means,col="blue3",lwd=3)



#-------------------------------------------------------------
# plot projections
#-------------------------------------------------------------

# make a multi-panel figure showing population size through time for different mutation models
setwd("~/Documents/orca/analyses_31May2022/Projections")


shapes <- 0.2
scales <- 0.1
propLethals <- 0.02
betas <- 13
betas2 <- as.character(betas)

infiles <- paste("popSize_shape",shapes,"_scale",scales,"_propLethal",propLethals,"_beta",betas,sep="")

noInbDepDat <- read.table("popSize_noInbDep",header=TRUE)[1:200,1:100]
noInbDepDat <- cbind(rep(74,nrow(noInbDepDat)),noInbDepDat)

for(i in 1:length(infiles)){
  popSize <- read.table(infiles[i],header=TRUE)
  # plot results
  library(scales)
  newPopSize <- popSize[1:200,1:100]                 # get rid of the extra rows where this are not data
  
  newPopSize <- cbind(rep(74,nrow(newPopSize)),newPopSize)
  plot(c(0,100),c(0,200),xlab="Year",ylab="Population Size",type="n",cex.lab=1.3)
  
  par(xpd=TRUE)
  text(x=-40,y=220,labels="B",cex=1.5)
  par(xpd=FALSE)
  #---------------------------------------------------------------------
  # plot population sizes for simulations without inbreeding depression
  #---------------------------------------------------------------------
  for(j in 1:nrow(noInbDepDat)){
    lines(0:100,as.numeric(noInbDepDat[j,]),col=alpha("darkgray",alpha=0.02))
  }
  medians <- apply(noInbDepDat,2,median)
  quant95  <- matrix(NA,nrow=2,ncol=length(medians))
  quant50  <- matrix(NA,nrow=2,ncol=length(medians))
  
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
  quant95  <- matrix(NA,nrow=2,ncol=length(medians))
  quant50  <- matrix(NA,nrow=2,ncol=length(medians))
  
  for(j in 1:ncol(newPopSize)){
    quant95[,j] <- quantile(newPopSize[,j],probs=c(0.025,0.975))
    quant50[,j] <- quantile(newPopSize[,j],probs=c(0.25,0.75))
  }
  polygon(c(0:100,rev(0:100)),c(quant95[1,1:101],rev(quant95[2,1:101])),col=alpha("darkred",alpha=0.1),border=NA)
  polygon(c(0:100,rev(0:100)),c(quant50[1,1:101],rev(quant50[2,1:101])),col=alpha("darkred",alpha=0.3),border=NA)
  lines(0:100,medians[1:101],col="darkred",lwd=3)
  legend(x=0,y=170,xjust=FALSE,yjust=FALSE,legend=c("No inbreeding depression",
    "Inbreeding depression"),bty="n",lwd=2,col=c("darkgray","darkred"),cex=0.8)
}
