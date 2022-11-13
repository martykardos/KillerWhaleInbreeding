# plot probability of surviving to 40 versus Froh

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
par(mfrow=c(2,2),mar=c(5,5,2,1))


# male 1Mb
plot(male1Mb$Froh_raw,male1Mb$pred_invlogit,ylim=yLims,xlab="",
     ylab="Male survival to 40 years",type="n",cex.lab=1.3)

text(x=-0.035,y=1.1,labels="A",cex=2)
frohs <- sort(unique(male1Mb$Froh_raw))
means <- rep(NA,length(frohs))          # plot median across 5000 MCMC draws
for(i in 1:length(means)){
  means[i] <- median(male1Mb$pred_invlogit[which(male1Mb$Froh_raw == frohs[i])])
}

survPreds <- matrix(NA,nrow=5000,ncol=length(unique(frohs)))            # plot 5000 MCMC draws
for(i in 1:5000){
  thisRep <- male1Mb[which(male1Mb$draw == i),]
  survPreds[i,] <- thisRep$pred_invlogit
  lines(thisRep$Froh_raw,thisRep$pred_invlogit,col=alpha("blue",alpha=0.008),lwd=0.5)
}


#plot 50% and 95% central bands
quantSurv95  <- matrix(NA,nrow=2,ncol=length(frohs))
quantSurv50  <- matrix(NA,nrow=2,ncol=length(frohs))
for(i in 1:length(frohs)){
  quantSurv95[,i] <- quantile(survPreds[,i],probs=c(0.025,0.975))
  quantSurv50[,i] <- quantile(survPreds[,i],probs=c(0.25,0.75))
}
polygon(c(frohs,rev(frohs)),c(quantSurv95[1,],rev(quantSurv95[2,])),col=alpha("blue3",alpha=0.05),border=NA)
polygon(c(frohs,rev(frohs)),c(quantSurv50[1,],rev(quantSurv50[2,])),col=alpha("blue3",alpha=0.2),border=NA)
lines(frohs,means,col="blue3",lwd=3)




# male 10Mb
plot(male10Mb$Froh_raw,male10Mb$pred_invlogit,ylim=yLims,xlab="",
     ylab="",type="n",cex.lab=1.3)

text(x=-0.035,y=1.1,labels="A",cex=2)
frohs <- sort(unique(male10Mb$Froh_raw))
means <- rep(NA,length(frohs))       # plot median across 5000 MCMC draws
for(i in 1:length(means)){
  means[i] <- median(male10Mb$pred_invlogit[which(male10Mb$Froh_raw == frohs[i])])
}

survPreds <- matrix(NA,nrow=5000,ncol=length(unique(frohs)))          # plot 5000 MCMC draws
for(i in 1:5000){
  thisRep <- male10Mb[which(male10Mb$draw == i),]
  survPreds[i,] <- thisRep$pred_invlogit
  lines(thisRep$Froh_raw,thisRep$pred_invlogit,col=alpha("blue",alpha=0.008),lwd=0.5)
}

#plot 50% and 95% central bands
quantSurv95  <- matrix(NA,nrow=2,ncol=length(frohs))
quantSurv50  <- matrix(NA,nrow=2,ncol=length(frohs))
for(i in 1:length(frohs)){
  quantSurv95[,i] <- quantile(survPreds[,i],probs=c(0.025,0.975))
  quantSurv50[,i] <- quantile(survPreds[,i],probs=c(0.25,0.75))
}
polygon(c(frohs,rev(frohs)),c(quantSurv95[1,],rev(quantSurv95[2,])),col=alpha("blue3",alpha=0.05),border=NA)
polygon(c(frohs,rev(frohs)),c(quantSurv50[1,],rev(quantSurv50[2,])),col=alpha("blue3",alpha=0.2),border=NA)
lines(frohs,means,col="blue3",lwd=3)


# female 1Mb
plot(female1Mb$Froh_raw,female1Mb$pred_invlogit,ylim=yLims,xlab=expression(italic(""*F*"")["ROH,1Mb"]),
     ylab="Female survival to 40 years",type="n",cex.lab=1.3)

text(x=-0.035,y=1.1,labels="A",cex=2)
frohs <- sort(unique(female1Mb$Froh_raw))
means <- rep(NA,length(frohs))          # plot median across 5000 MCMC draws
for(i in 1:length(means)){
  means[i] <- median(female1Mb$pred_invlogit[which(female1Mb$Froh_raw == frohs[i])])
}

survPreds <- matrix(NA,nrow=5000,ncol=length(unique(frohs)))             # plot 5000 MCMC draws
for(i in 1:5000){
  thisRep <- female1Mb[which(female1Mb$draw == i),]
  survPreds[i,] <- thisRep$pred_invlogit
  lines(thisRep$Froh_raw,thisRep$pred_invlogit,col=alpha("blue",alpha=0.008),lwd=0.5)
}

#plot 50% and 95% central bands
quantSurv95  <- matrix(NA,nrow=2,ncol=length(frohs))
quantSurv50  <- matrix(NA,nrow=2,ncol=length(frohs))
for(i in 1:length(frohs)){
  quantSurv95[,i] <- quantile(survPreds[,i],probs=c(0.025,0.975))
  quantSurv50[,i] <- quantile(survPreds[,i],probs=c(0.25,0.75))
}
polygon(c(frohs,rev(frohs)),c(quantSurv95[1,],rev(quantSurv95[2,])),col=alpha("blue3",alpha=0.05),border=NA)
polygon(c(frohs,rev(frohs)),c(quantSurv50[1,],rev(quantSurv50[2,])),col=alpha("blue3",alpha=0.2),border=NA)
lines(frohs,means,col="blue3",lwd=3)





# female 10Mb
plot(female10Mb$Froh_raw,female10Mb$pred_invlogit,ylim=yLims,xlab=expression(italic(""*F*"")["ROH,10Mb"]),
     ylab="",type="n",cex.lab=1.3)

text(x=-0.035,y=1.1,labels="A",cex=2)
frohs <- sort(unique(female10Mb$Froh_raw))
means <- rep(NA,length(frohs))                # plot median across 5000 MCMC draws
for(i in 1:length(means)){
  means[i] <- median(female10Mb$pred_invlogit[which(female10Mb$Froh_raw == frohs[i])])
}

survPreds <- matrix(NA,nrow=5000,ncol=length(unique(frohs)))         # plot 5000 MCMC draws
for(i in 1:5000){
  thisRep <- female10Mb[which(female10Mb$draw == i),]
  survPreds[i,] <- thisRep$pred_invlogit
  lines(thisRep$Froh_raw,thisRep$pred_invlogit,col=alpha("blue",alpha=0.008),lwd=0.5)
}

#plot 50% and 95% central bands
quantSurv95  <- matrix(NA,nrow=2,ncol=length(frohs))
quantSurv50  <- matrix(NA,nrow=2,ncol=length(frohs))
for(i in 1:length(frohs)){
  quantSurv95[,i] <- quantile(survPreds[,i],probs=c(0.025,0.975))
  quantSurv50[,i] <- quantile(survPreds[,i],probs=c(0.25,0.75))
}
polygon(c(frohs,rev(frohs)),c(quantSurv95[1,],rev(quantSurv95[2,])),col=alpha("blue3",alpha=0.05),border=NA)
polygon(c(frohs,rev(frohs)),c(quantSurv50[1,],rev(quantSurv50[2,])),col=alpha("blue3",alpha=0.2),border=NA)
lines(frohs,means,col="blue3",lwd=3)

