
###########################################################
# PLOT SURVIVAL PROBABILITY AGAINST INBREEDING
###########################################################
library(scales)
Shape <- 0.2
Scale <- 0.1
PropLethal <- 0
beta <- 50
SimSex <- "male"
setwd("~/Documents/orca/analyses_31May2022/Projections")
outSurvProbs <- read.table(paste("survTo40Probs_shape",Shape,"_scale",Scale,"_propLethal",PropLethal,"_beta",beta,"_sex_",SimSex,sep=""),header=TRUE)
outFh <- read.table(paste("Fh_shape",Shape,"_scale",Scale,"_propLethal",PropLethal,"_beta",beta,"_sex_",SimSex,sep=""),header=TRUE)
plot(c(0,0.36),c(0,1),xlab=expression(italic(""*F*"")[h]),ylab="Male survival so 40 years",type="n")
for(i in 1:nrow(outFh)){
  theseProbs <- outSurvProbs[i,]
  if(sum(theseProbs > 0.999) > 0) theseProbs[which(theseProbs > 0.999)] <- 0.999
  lines(lowess(outFh[i,],theseProbs,f=0.66),col=alpha("blue",alpha=0.05),lwd=1)
}

# make vectors of inbreeding and survival probability
allFhVec <- NULL
allSurvProbVec <- NULL
for(i in 1:nrow(outFh)){
  allFhVec <- c(allFhVec,outFh[i,])
  theseSurvProbs <- outSurvProbs[i,]
  if(sum(theseSurvProbs > 0.999) > 0)theseSurvProbs[which(theseSurvProbs > 0.999)] <- 0.999
  allSurvProbVec <- c(allSurvProbVec,theseSurvProbs)
}
lines(lowess(allFhVec,allSurvProbVec),lwd=5)




###########################################################
# PLOT SURVIVAL PROBABILITY AGAINST INBREEDING
###########################################################
par(mfrow=c(1,1),xpd=FALSE)
Shape <- 0.2
Scale <- 0.1
PropLethal <- 0
beta <- 50
SimSex <- "male"
setwd("~/Documents/orca/analyses_31May2022/Projections")
outSurvProbs <- read.table(paste("survTo40Probs_shape",Shape,"_scale",Scale,"_propLethal",PropLethal,"_beta",beta,"_sex_",SimSex,sep=""),header=TRUE)
outFh <- read.table(paste("het_shape",Shape,"_scale",Scale,"_propLethal",PropLethal,"_beta",beta,"_sex_",SimSex,sep=""),header=TRUE)
plot(c(0.19,0.31),c(0,1),xlab=expression(italic(""*H*"")[SNP]),ylab="Male survival so 40 years",type="n")
for(i in 1:nrow(outFh)){
  theseProbs <- outSurvProbs[i,]
  if(sum(theseProbs > 0.999) > 0) theseProbs[which(theseProbs > 0.999)] <- 0.999
  lines(lowess(outFh[i,],theseProbs,f=0.66),col=alpha("blue",alpha=0.05),lwd=1)
}

# make vectors of inbreeding and survival probability
allFhVec <- NULL
allSurvProbVec <- NULL
for(i in 1:nrow(outFh)){
  allFhVec <- c(allFhVec,outFh[i,])
  theseSurvProbs <- outSurvProbs[i,]
  if(sum(theseSurvProbs > 0.999) > 0)theseSurvProbs[which(theseSurvProbs > 0.999)] <- 0.999
  allSurvProbVec <- c(allSurvProbVec,theseSurvProbs)
}
lines(lowess(allFhVec,allSurvProbVec),lwd=5)