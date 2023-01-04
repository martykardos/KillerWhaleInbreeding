# code for fig S19
setwd("~/Documents/orca/analyses_31May2022/demograhpy")
library(tidyverse)
library(viridis)
fig_19 <- readRDS("~/Documents/orca/analyses_31May2022/demograhpy/fig_19.rds")
dat1Mb <- fig_19[which(fig_19$cov == "1Mb" & fig_19$model == "Absolute" & fig_19$Froh_raw >= 0.19 & fig_19$Froh_raw <= 0.45),]
plot(dat1Mb$Froh_raw,dat1Mb$lrs,xlim=c(0,0.5),type="n",ylim=c(1.5,2.6),ylab="Lifetime reproductive success",xlab=expression(italic(""*F*"")[ROH]),cex.lab=1.3)
lines(dat1Mb$Froh_raw,dat1Mb$lrs,col="blue",lwd=3)
dat10Mb <- fig_19[which(fig_19$cov == "10Mb" & fig_19$model == "Absolute" & fig_19$Froh_raw >= 0 & fig_19$Froh_raw <= 0.145),]
lines(dat10Mb$Froh_raw,dat10Mb$lrs,col="darkred",lwd=3)
legend(x=0.35,y=2.5,legend=c(expression(italic(""*F*"")[ROH1Mb]),expression(italic(""*F*"")[ROH10Mb])),lty="solid",lwd=2,col=c("blue","darkred"),bty="n")




