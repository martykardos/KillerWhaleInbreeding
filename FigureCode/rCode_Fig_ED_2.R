# plot male annual survival versus Froh
setwd("/Users/martin.kardos/Documents/orca/analyses_31May2022/demograhpy")
fig_s17 <- as.data.frame(readRDS("fig_s17.rds"))

#################################################
# Fig S17
#################################################
par(mfrow=c(2,2),mar=c(5,5,3,1),xpd=FALSE)
library(scales)
yLims = c(0.9,1)

# Male 1Mb
thisDat <- fig_s17[which(fig_s17$name == "Male 1Mb"),]
thisDat <- thisDat[order(thisDat$draws),]
frohs <- sort(unique(thisDat$Froh_raw[which(thisDat$draws == 1)]))
survPreds <- matrix(NA,nrow=5000,ncol=length(frohs))
for(i in 1:5000){                 # show 5000 MCMC draws
  thisRep <- thisDat[thisDat$draws == i,]
  thisRep <- thisRep[order(thisRep$Froh_raw),]
  thisRep <- thisRep[-which(duplicated(thisRep$Froh_raw)),]
  survPreds[i,] <- thisRep$pred_invlogit
}
means <- apply(survPreds,2,median)# plot median estimate
par(xpd=FALSE)

plot(frohs,survPreds[1,],ylim=yLims,xlab="",ylab="Male yearly survival",type="n",cex.lab=1.5)
for(i in 1:5000){     # show 5000 MCMC draws
  lines(frohs,survPreds[i,],col=alpha("blue",alpha=0.008),lwd=0.5)
}

# plot central 50% and 95% intervals
quantSurv95  <- matrix(NA,nrow=2,ncol=length(frohs))
quantSurv50  <- matrix(NA,nrow=2,ncol=length(frohs))
for(i in 1:length(frohs)){
  quantSurv95[,i] <- quantile(survPreds[,i],probs=c(0.025,0.975))
  quantSurv50[,i] <- quantile(survPreds[,i],probs=c(0.25,0.75))
}
polygon(c(frohs,rev(frohs)),c(quantSurv95[1,],rev(quantSurv95[2,])),col=alpha("blue3",alpha=0.05),border=NA)
polygon(c(frohs,rev(frohs)),c(quantSurv50[1,],rev(quantSurv50[2,])),col=alpha("blue3",alpha=0.2),border=NA)
lines(frohs,means,col="blue3",lwd=3)

# Male 10Mb
yLims = c(0.75,1)
par(xpd=FALSE)
thisDat <- fig_s17[which(fig_s17$name == "Male 10Mb"),]
thisDat <- thisDat[order(thisDat$draws),]
frohs <- (sort(unique(thisDat$Froh_raw[which(thisDat$draws == 1)])-min(thisDat$Froh_raw[which(thisDat$draws == 1)]))/max(unique(thisDat$Froh_raw[which(thisDat$draws == 1)])-min(thisDat$Froh_raw[which(thisDat$draws == 1)]))) * 0.145716671
survPreds <- matrix(NA,nrow=5000,ncol=length(frohs))
for(i in 1:5000){# show 5000 MCMC draws
  thisRep <- thisDat[thisDat$draws == i,]
  thisRep <- thisRep[order(thisRep$Froh_raw),]
  thisRep <- thisRep[-which(duplicated(thisRep$Froh_raw)),]
  survPreds[i,] <- thisRep$pred_invlogit
}
means <- apply(survPreds,2,median)# plot median estimate

par(xpd=FALSE)
plot(frohs,survPreds[1,],ylim=yLims,xlab="",ylab="",type="n")
for(i in 1:5000){       # show 5000 MCMC draws
  lines(frohs,survPreds[i,],col=alpha("blue",alpha=0.008),lwd=0.5)
}

# plot central 50% and 95% intervals
quantSurv95  <- matrix(NA,nrow=2,ncol=length(frohs))
quantSurv50  <- matrix(NA,nrow=2,ncol=length(frohs))
for(i in 1:length(frohs)){
  quantSurv95[,i] <- quantile(survPreds[,i],probs=c(0.025,0.975))
  quantSurv50[,i] <- quantile(survPreds[,i],probs=c(0.25,0.75))
}
polygon(c(frohs,rev(frohs)),c(quantSurv95[1,],rev(quantSurv95[2,])),col=alpha("blue3",alpha=0.05),border=NA)
polygon(c(frohs,rev(frohs)),c(quantSurv50[1,],rev(quantSurv50[2,])),col=alpha("blue3",alpha=0.2),border=NA)
lines(frohs,means,col="blue3",lwd=3)


# Female 1Mb
yLims = c(0.9,1)
thisDat <- fig_s17[which(fig_s17$name == "Female 1Mb"),]
thisDat <- thisDat[order(thisDat$draws),]
frohs <- sort(unique(thisDat$Froh_raw[which(thisDat$draws == 1)]))
survPreds <- matrix(NA,nrow=5000,ncol=length(frohs))
for(i in 1:5000){# show 5000 MCMC draws
  thisRep <- thisDat[thisDat$draws == i,]
  thisRep <- thisRep[order(thisRep$Froh_raw),]
  thisRep <- thisRep[-which(duplicated(thisRep$Froh_raw)),]
  survPreds[i,] <- thisRep$pred_invlogit
}
means <- apply(survPreds,2,median)# plot median estimate

plot(frohs,survPreds[1,],ylim=yLims,xlab=expression(italic(""*F*"")["ROH,1Mb"]),ylab="Female yearly survival",type="n",cex.lab=1.5)
for(i in 1:5000){      # show 5000 MCMC draws
  lines(frohs,survPreds[i,],col=alpha("blue",alpha=0.006),lwd=0.5)
}

# plot central 50% and 95% intervals
quantSurv95  <- matrix(NA,nrow=2,ncol=length(frohs))
quantSurv50  <- matrix(NA,nrow=2,ncol=length(frohs))
for(i in 1:length(frohs)){
  quantSurv95[,i] <- quantile(survPreds[,i],probs=c(0.025,0.975))
  quantSurv50[,i] <- quantile(survPreds[,i],probs=c(0.25,0.75))
}
polygon(c(frohs,rev(frohs)),c(quantSurv95[1,],rev(quantSurv95[2,])),col=alpha("blue3",alpha=0.05),border=NA)
polygon(c(frohs,rev(frohs)),c(quantSurv50[1,],rev(quantSurv50[2,])),col=alpha("blue3",alpha=0.2),border=NA)
lines(frohs,means,col="blue3",lwd=3)


# Female 10Mb
yLims = c(0.75,1)
thisDat <- fig_s17[which(fig_s17$name == "Female 10Mb"),]
thisDat <- thisDat[order(thisDat$draws),]
frohs <- (sort(unique(thisDat$Froh_raw[which(thisDat$draws == 1)])-min(thisDat$Froh_raw[which(thisDat$draws == 1)]))/max(unique(thisDat$Froh_raw[which(thisDat$draws == 1)])-min(thisDat$Froh_raw[which(thisDat$draws == 1)]))) * 0.145716671

survPreds <- matrix(NA,nrow=5000,ncol=length(frohs))
for(i in 1:5000){# show 5000 MCMC draws
  thisRep <- thisDat[thisDat$draws == i,]
  thisRep <- thisRep[order(thisRep$Froh_raw),]
  thisRep <- thisRep[-which(duplicated(thisRep$Froh_raw)),]
  survPreds[i,] <- thisRep$pred_invlogit
}
means <- apply(survPreds,2,median)        # plot median estimate
par(xpd=FALSE)
plot(frohs,survPreds[1,],ylim=yLims,xlab=expression(italic(""*F*"")["ROH,10Mb"]),ylab="",type="n",cex.lab=1.5)
for(i in 1:5000){     # show 5000 MCMC draws
  lines(frohs,survPreds[i,],col=alpha("blue",alpha=0.006),lwd=0.5)
}

# plot central 50% and 95% intervals
quantSurv95  <- matrix(NA,nrow=2,ncol=length(frohs))
quantSurv50  <- matrix(NA,nrow=2,ncol=length(frohs))
for(i in 1:length(frohs)){
  quantSurv95[,i] <- quantile(survPreds[,i],probs=c(0.025,0.975))
  quantSurv50[,i] <- quantile(survPreds[,i],probs=c(0.25,0.75))
}
polygon(c(frohs,rev(frohs)),c(quantSurv95[1,],rev(quantSurv95[2,])),col=alpha("blue3",alpha=0.05),border=NA)
polygon(c(frohs,rev(frohs)),c(quantSurv50[1,],rev(quantSurv50[2,])),col=alpha("blue3",alpha=0.2),border=NA)
lines(frohs,means,col="blue3",lwd=3)


