####################################################################
# plot survival versus heterozygosity
####################################################################
setwd("/Users/martin.kardos/Documents/orca/analyses_31May2022/demograhpy")
fig_sX <- as.data.frame(readRDS("cumulative_survival_vs_het.rds"))
#fig_sX$Froh_raw [which(fig_sX$Froh_raw > 0.315 | fig_sX$Froh_raw < 0.195)] <- NA

male <- fig_sX[which(fig_sX$Sex == "Male"),]
male <- male[order(male$draw),]

female <- fig_sX[which(fig_sX$Sex == "Female"),]
female <- female[order(female$draw),]

xLims <- c(0.196,0.314)
yLims <- c(0,1)

library(scales)
par(mfrow=c(1,2),mar=c(5,5,2,1))


# males
plot(male$Froh_raw,male$pred_invlogit,ylim=yLims,xlim=xLims,xlab=expression(italic(""*H*"")[SNP]),
     ylab="Male survival to 40 years",type="n",cex.lab=1.3)

frohs <- sort(unique(male$Froh_raw))
means <- rep(NA,length(frohs))
for(i in 1:length(means)){
  means[i] <- median(male$pred_invlogit[which(male$Froh_raw == frohs[i])])
}

survPreds <- matrix(NA,nrow=5000,ncol=length(unique(frohs)))
for(i in 1:5000){
  thisRep <- male[which(male$draw == i),]
  survPreds[i,] <- thisRep$pred_invlogit
}

# truncate the plot at the observed range of H values
newFrohs <- frohs[which(frohs < 0.315 & frohs > 0.195)]
newSurvPreds <- survPreds[,which(frohs < 0.315 & frohs > 0.195)]
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
rect(xleft=0.194,xright=0.1962,ybottom=0.0,ytop=1,border=NA,col="white")
rect(xleft=0.3143,xright=0.316,ybottom=0.0,ytop=1,border=NA,col="white")


# male lethal eqs.

FhMax <- (0.315-0.195)/0.315
slope <- (log(0.7799869)-log(0.1875211))/FhMax
fSlope <- (0.7799869-0.1875211)/(FhMax)
hetSlope <- (0.7799869-0.1875211)/(0.314-0.196)

# females
plot(female$Froh_raw,female$pred_invlogit,ylim=yLims,xlim=xLims,xlab=expression(italic(""*H*"")[SNP]),
     ylab="Female survival to 40 years",type="n",cex.lab=1.3)

frohs <- sort(unique(female$Froh_raw))
means <- rep(NA,length(frohs))
for(i in 1:length(means)){
  means[i] <- median(female$pred_invlogit[which(female$Froh_raw == frohs[i])])
}

survPreds <- matrix(NA,nrow=5000,ncol=length(unique(frohs)))
for(i in 1:5000){
  thisRep <- female[which(female$draw == i),]
  survPreds[i,] <- thisRep$pred_invlogit
}

# truncate the plot at the observed range of H values
newFrohs <- frohs[which(frohs < 0.315 & frohs > 0.195)]
newSurvPreds <- survPreds[,which(frohs < 0.315 & frohs > 0.195)]
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
rect(xleft=0.194,xright=0.1962,ybottom=0.0,ytop=1,border=NA,col="white")
rect(xleft=0.3143,xright=0.316,ybottom=0.0,ytop=1,border=NA,col="white")

# female lethal eqs

FhMax <- (0.315-0.195)/0.315
slope <- (log(0.84312)-log(0.29705))/FhMax





