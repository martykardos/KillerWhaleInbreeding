setwd("/Users/martin.kardos/Documents/orca/analyses_31May2022/RXY")
################################################
# plot the sfs
################################################
sfsSRKW <- read.table("sfs_SRKW_downSamp_12",header=TRUE)
sfsARKW <- read.table("sfs_ARKW_downSamp_12",header=TRUE)
sfsTKW <- read.table("sfs_TKW_downSamp_12",header=TRUE)

par(mfrow=c(1,3))
library(scales)

# SRKW
plot(c(0,1),c(0,0.22),type="n",ylab="Proportion of Loci",xlab="Derived Allele Frequency",main="Southern Residents",cex.lab=1.3)
lines(sfsSRKW[,1]/(12*2),sfsSRKW[,2]/sum(sfsSRKW[,2]),col=alpha("orange",alpha=0.5),lwd=3)
lines(sfsSRKW[,1]/(12*2),(sfsSRKW[,3]+sfsSRKW[,4])/sum(sfsSRKW[,3]+sfsSRKW[,4]),col=alpha("blue",alpha=0.5),lwd=3)

legend(x=0,y=0.17,xjust=FALSE,yjust=FALSE,legend=c("Deleterious","Neutral"),lty="solid",
       lwd=2,col=c(alpha("blue",alpha=0.5),alpha("orange",alpha=0.5)),bty="n")


# ARKW
plot(c(0,1),c(0,0.22),type="n",ylab="",xlab="Derived Allele Frequency",main="Alaska Residents",cex.lab=1.3)
lines(sfsARKW[,1]/(12*2),sfsARKW[,2]/sum(sfsARKW[,2]),col=alpha("orange",alpha=0.5),lwd=3)
lines(sfsARKW[,1]/(12*2),(sfsARKW[,3]+sfsARKW[,4])/sum(sfsARKW[,3]+sfsARKW[,4]),col=alpha("blue",alpha=0.5),lwd=3)

# TKW
plot(c(0,1),c(0,0.22),type="n",ylab="",xlab="Derived Allele Frequency",main="Transients",cex.lab=1.3)
lines(sfsTKW[,1]/(12*2),sfsTKW[,2]/sum(sfsTKW[,2]),col=alpha("orange",alpha=0.5),lwd=3)
lines(sfsTKW[,1]/(12*2),(sfsTKW[,3]+sfsTKW[,4])/sum(sfsTKW[,3]+sfsTKW[,4]),col=alpha("blue",alpha=0.5),lwd=3)


