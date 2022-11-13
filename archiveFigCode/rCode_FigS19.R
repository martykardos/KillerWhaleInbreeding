####################################################################
# plot survival versus heterozygosity
####################################################################
setwd("/Users/martin.kardos/Documents/orca/analyses_31May2022/demograhpy")
fig_sX <- as.data.frame(readRDS("20_year_survival_vs_het.rds"))

male <- fig_sX[which(fig_sX$sexF1M2 == 2),]
male <- male[order(male$draw),]


female <- fig_sX[which(fig_sX$sexF1M2 == 1),]
female <- female[order(female$draw),]

xLims <- c(0.19,0.318)
yLims <- c(0,1)

library(scales)
par(mfrow=c(1,2),mar=c(5,5,2,1))


# males
plot(male$Froh_raw,male$pred_invlogit,ylim=c(0.8,1),xlim=xLims,xlab=expression(italic(""*H*"")),
     ylab="Male annual survival",type="n",cex.lab=1.3)

text(x=-0.035,y=1.1,labels="A",cex=2)
frohs <- sort(unique(male$Froh_raw))
means <- rep(NA,length(frohs))
for(i in 1:length(means)){
  means[i] <- median(male$pred_invlogit[which(male$Froh_raw == frohs[i])])
}

survPreds <- matrix(NA,nrow=5000,ncol=23)
for(i in 1:5000){
  thisRep <- male[which(male$draw == i),]
  thisRep <- thisRep[-nrow(thisRep),]
  thisRep <- thisRep[order(thisRep$Froh_raw),]
  survPreds[i,] <- thisRep$pred_invlogit
}

# truncate the plot at the observed range of H values
newFrohs <- frohs[which(frohs < 0.318 & frohs > 0.194)]
newSurvPreds <- survPreds[,which(frohs < 0.318 & frohs > 0.194)]
for(i in 1:5000){
  lines(newFrohs,newSurvPreds[i,],col=alpha("blue",alpha=0.008),lwd=0.5)
}
means <- colMeans(newSurvPreds)
quantSurv95  <- matrix(NA,nrow=2,ncol=length(newFrohs))
quantSurv50  <- matrix(NA,nrow=2,ncol=length(newFrohs))

for(i in 1:length(newFrohs)){
  quantSurv95[,i] <- quantile(newSurvPreds[,i],probs=c(0.025,0.975))
  quantSurv50[,i] <- quantile(newSurvPreds[,i],probs=c(0.25,0.75))
}
polygon(c(newFrohs,rev(newFrohs)),c(quantSurv95[1,],rev(quantSurv95[2,])),col=alpha("blue3",alpha=0.05),border=NA)
polygon(c(newFrohs,rev(newFrohs)),c(quantSurv50[1,],rev(quantSurv50[2,])),col=alpha("blue3",alpha=0.2),border=NA)
lines(newFrohs,means,col="blue3",lwd=3)

rect(xleft=0.19,xright=0.1962,ybottom=0.8,ytop=1,border=NA,col="white")     # limit plot range to empirical range
rect(xleft=0.3143,xright=0.32,ybottom=0.8,ytop=1,border=NA,col="white")






# females
plot(female$Froh_raw,female$pred_invlogit,ylim=c(0.8,1),xlim=xLims,xlab=expression(italic(""*H*"")),
     ylab="Female annual survival",type="n",cex.lab=1.3)

text(x=-0.035,y=1.1,labels="A",cex=2)
frohs <- sort(unique(female$Froh_raw))
means <- rep(NA,length(frohs))
for(i in 1:length(means)){
  means[i] <- median(female$pred_invlogit[which(female$Froh_raw == frohs[i])])
}

survPreds <- matrix(NA,nrow=5000,ncol=23)
for(i in 1:5000){
  thisRep <- female[which(female$draw == i),]
  thisRep <- thisRep[-nrow(thisRep),]
  thisRep <- thisRep[order(thisRep$Froh_raw),]
  survPreds[i,] <- thisRep$pred_invlogit

}

# truncate the plot at the observed range of H values
newFrohs <- frohs[which(frohs < 0.318 & frohs > 0.194)]
newSurvPreds <- survPreds[,which(frohs < 0.318 & frohs > 0.194)]
for(i in 1:5000){
  lines(newFrohs,newSurvPreds[i,],col=alpha("blue",alpha=0.008),lwd=0.5)
}
means <- colMeans(newSurvPreds)
quantSurv95  <- matrix(NA,nrow=2,ncol=length(newFrohs))
quantSurv50  <- matrix(NA,nrow=2,ncol=length(newFrohs))

for(i in 1:length(newFrohs)){
  quantSurv95[,i] <- quantile(newSurvPreds[,i],probs=c(0.025,0.975))
  quantSurv50[,i] <- quantile(newSurvPreds[,i],probs=c(0.25,0.75))
}
polygon(c(newFrohs,rev(newFrohs)),c(quantSurv95[1,],rev(quantSurv95[2,])),col=alpha("blue3",alpha=0.05),border=NA)
polygon(c(newFrohs,rev(newFrohs)),c(quantSurv50[1,],rev(quantSurv50[2,])),col=alpha("blue3",alpha=0.2),border=NA)
lines(newFrohs,means,col="blue3",lwd=3)

rect(xleft=0.19,xright=0.1962,ybottom=0.8,ytop=1,border=NA,col="white")      # limit plot range to empirical range
rect(xleft=0.3143,xright=0.32,ybottom=0.8,ytop=1,border=NA,col="white")









####################################################################
# plot survival versus Fh
####################################################################
setwd("/Users/martin.kardos/Documents/orca/analyses_31May2022/demograhpy")
fig_sX <- as.data.frame(readRDS("20_year_survival_vs_het.rds"))

male <- fig_sX[which(fig_sX$sexF1M2 == 2),]
male <- male[order(male$draw),]


female <- fig_sX[which(fig_sX$sexF1M2 == 1),]
female <- female[order(female$draw),]

xLims <- c(0,0.38)
yLims <- c(0,1)

library(scales)
par(mfrow=c(1,2),mar=c(5,5,2,1))


# males

text(x=-0.035,y=1.1,labels="A",cex=2)
frohs <- sort(unique(male$Froh_raw))
means <- rep(NA,length(frohs))
for(i in 1:length(means)){
  means[i] <- median(male$pred_invlogit[which(male$Froh_raw == frohs[i])])
}

survPreds <- matrix(NA,nrow=5000,ncol=23)
for(i in 1:5000){
  thisRep <- male[which(male$draw == i),]
  thisRep <- thisRep[-nrow(thisRep),]
  thisRep <- thisRep[order(thisRep$Froh_raw),]
  survPreds[i,] <- thisRep$pred_invlogit
  #lines(thisRep$Froh_raw,thisRep$pred_invlogit,col=alpha("blue",alpha=0.008),lwd=0.5)
}

# truncate the plot at the observed range of H values
newFrohs <- frohs[which(frohs < 0.318 & frohs > 0.194)]
newFrohs <- (max(newFrohs)-newFrohs)/max(newFrohs)
newSurvPreds <- survPreds[,which(frohs < 0.318 & frohs > 0.194)]


plot(male$Froh_raw,male$pred_invlogit,ylim=c(0.8,1),xlim=xLims,xlab=expression(italic(""*F*"")[h]),
     ylab="Male annual survival",type="n",cex.lab=1.3)

for(i in 1:5000){
  lines(newFrohs,newSurvPreds[i,],col=alpha("blue",alpha=0.008),lwd=0.5)
}
means <- colMeans(newSurvPreds)
quantSurv95  <- matrix(NA,nrow=2,ncol=length(newFrohs))
quantSurv50  <- matrix(NA,nrow=2,ncol=length(newFrohs))

for(i in 1:length(newFrohs)){
  quantSurv95[,i] <- quantile(newSurvPreds[,i],probs=c(0.025,0.975))
  quantSurv50[,i] <- quantile(newSurvPreds[,i],probs=c(0.25,0.75))
}
polygon(c(newFrohs,rev(newFrohs)),c(quantSurv95[1,],rev(quantSurv95[2,])),col=alpha("blue3",alpha=0.05),border=NA)
polygon(c(newFrohs,rev(newFrohs)),c(quantSurv50[1,],rev(quantSurv50[2,])),col=alpha("blue3",alpha=0.2),border=NA)
lines(newFrohs,means,col="blue3",lwd=3)






# females


text(x=-0.035,y=1.1,labels="A",cex=2)
frohs <- sort(unique(female$Froh_raw))
means <- rep(NA,length(frohs))
for(i in 1:length(means)){
  means[i] <- median(female$pred_invlogit[which(female$Froh_raw == frohs[i])])
}

survPreds <- matrix(NA,nrow=5000,ncol=23)
for(i in 1:5000){
  thisRep <- female[which(female$draw == i),]
  thisRep <- thisRep[-nrow(thisRep),]
  thisRep <- thisRep[order(thisRep$Froh_raw),]
  survPreds[i,] <- thisRep$pred_invlogit
  
}

# truncate the plot at the observed range of H values
newFrohs <- frohs[which(frohs < 0.318 & frohs > 0.194)]
newFrohs <- (max(newFrohs)-newFrohs)/max(newFrohs)
newSurvPreds <- survPreds[,which(frohs < 0.318 & frohs > 0.194)]

plot(female$Froh_raw,female$pred_invlogit,ylim=c(0.8,1),xlim=xLims,xlab=expression(italic(""*F*"")[h]),
     ylab="Female annual survival",type="n",cex.lab=1.3)
for(i in 1:5000){
  lines(newFrohs,newSurvPreds[i,],col=alpha("blue",alpha=0.008),lwd=0.5)
}
means <- colMeans(newSurvPreds)
quantSurv95  <- matrix(NA,nrow=2,ncol=length(newFrohs))
quantSurv50  <- matrix(NA,nrow=2,ncol=length(newFrohs))

for(i in 1:length(newFrohs)){
  quantSurv95[,i] <- quantile(newSurvPreds[,i],probs=c(0.025,0.975))
  quantSurv50[,i] <- quantile(newSurvPreds[,i],probs=c(0.25,0.75))
}
polygon(c(newFrohs,rev(newFrohs)),c(quantSurv95[1,],rev(quantSurv95[2,])),col=alpha("blue3",alpha=0.05),border=NA)
polygon(c(newFrohs,rev(newFrohs)),c(quantSurv50[1,],rev(quantSurv50[2,])),col=alpha("blue3",alpha=0.2),border=NA)
lines(newFrohs,means,col="blue3",lwd=3)



