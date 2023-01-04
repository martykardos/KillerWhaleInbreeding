
#----------------------------------------------------------------------------------------------------------
# exclude individuals born before 1960 and plot males and females separately for each ROH length threshold
#----------------------------------------------------------------------------------------------------------
setwd("/Users/martin.kardos/Documents/orca/analyses_31May2022/ROH3")
deathDat <- read.table("ageAtDeath_inbreeding_fromMFord_29April2022",header=TRUE)
FrohDat <- read.table("Froh_7August2022",header=TRUE)
library(asbio)

# add froh to death data
deathDat <- cbind(deathDat,FrohDat[match(deathDat[,1],FrohDat[,1]),])
colnames(deathDat) <- c("id","age","sex","birthYear","deathYear","id2","pop","Froh_min1Mb","Froh_min10Mb")
deathDat <- deathDat[-which(deathDat$birthYear < 1960),]
par(mar=c(4,5,2,1),mfrow=c(2,2))

# 1Mb 
mod <- lm(deathDat$age~ deathDat$Froh_min1Mb + deathDat$sex)
summary(mod)
col<- rep("gray",nrow(deathDat))
col[which(deathDat$sex == "M")] <- "darkred"
library(scales)
yLims<- c(0,max(deathDat$age[which(deathDat$birt >= 1960)]))
xLims<- c(min(deathDat$Froh_min1Mb[which(deathDat$birt >= 1960)]),max(deathDat$Froh_min1Mb[which(deathDat$birt >= 1960)]))
plot(deathDat$Froh_min1Mb[which(deathDat$sex == "F")],deathDat$age[which(deathDat$sex == "F")],col=alpha("gray",alpha=0.6),
     pch=16,cex=2,ylab="Age at death",xlab=expression(italic(""*F*"")["ROH,1Mb"]),cex.lab=1.5,ylim=yLims,xlim=xLims)
x1 <- min(deathDat$Froh_min1Mb)
x2 <- max(deathDat$Froh_min1Mb)
int <- mod$coefficients[1]
slope <- mod$coefficients[2]
sex <- mod$coefficients[3]

# line for females
lines(c(x1,x2),c(int+slope*x1,int+slope*x2),col="gray",lwd=3)
par(xpd=TRUE)
text(x=mean(c(par("usr")[1],par("usr")[2])),y=59,labels=expression("Females"),cex=1.7)
par(xpd=FALSE)
plot(deathDat$Froh_min1Mb[which(deathDat$sex == "M")],deathDat$age[which(deathDat$sex == "M")],col=alpha("darkred",alpha=0.6),
     pch=16,cex=2,ylab="",xlab=expression(italic(""*F*"")["ROH,1Mb"]),cex.lab=1.5,ylim=yLims,xlim=xLims)

# line for males
lines(c(x1,x2),c(int+sex+slope*x1,int+sex+slope*x2),col="darkred",lwd=3)
par(xpd=TRUE)
text(x=mean(c(par("usr")[1],par("usr")[2])),y=59,labels=expression("Males"),cex=1.7)
par(xpd=FALSE)

# 10Mb 
mod <- lm(deathDat$age~ deathDat$Froh_min10Mb + deathDat$sex)
summary(mod)
col<- rep("gray",nrow(deathDat))
col[which(deathDat$sex == "M")] <- "darkred"
library(scales)
yLims<- c(0,max(deathDat$age[which(deathDat$birt >= 1960)]))
xLims<- c(min(deathDat$Froh_min10Mb[which(deathDat$birt >= 1960)]),max(deathDat$Froh_min10Mb[which(deathDat$birt >= 1960)]))
plot(deathDat$Froh_min10Mb[which(deathDat$sex == "F")],deathDat$age[which(deathDat$sex == "F")],col=alpha("gray",alpha=0.6),
     pch=16,cex=2,ylab="Age at death",xlab=expression(italic(""*F*"")["ROH,10Mb"]),cex.lab=1.5,ylim=yLims,xlim=xLims)
x1 <- min(deathDat$Froh_min10Mb)
x2 <- max(deathDat$Froh_min10Mb)
int <- mod$coefficients[1]
slope <- mod$coefficients[2]
sex <- mod$coefficients[3]

# line for females
lines(c(x1,x2),c(int+slope*x1,int+slope*x2),col="gray",lwd=3)
plot(deathDat$Froh_min10Mb[which(deathDat$sex == "M")],deathDat$age[which(deathDat$sex == "M")],col=alpha("darkred",alpha=0.6),
     pch=16,cex=2,ylab="",xlab=expression(italic(""*F*"")["ROH,10Mb"]),cex.lab=1.5,ylim=yLims,xlim=xLims)

# line for males
lines(c(x1,x2),c(int+sex+slope*x1,int+sex+slope*x2),col="darkred",lwd=3)



#----------------------------------------------------------------------------------------
#########################################################################################
# repeat this with genomic heterozygosity measured using called genotypes on the x axis
#########################################################################################
#----------------------------------------------------------------------------------------

setwd("~/Documents/orca/analyses_31May2022/RXY")
het <- read.table("srkw_indivHet",header=TRUE)
het <- het[match(deathDat[,1],het[,1]),]

deathDat <- cbind(deathDat,het[,c(3,7)])


# 10000 Loci 
mod <- lm(deathDat$age~ deathDat$X10000_loci + deathDat$sex)
summary(mod)
col<- rep("gray",nrow(deathDat))
col[which(deathDat$sex == "M")] <- "darkred"
library(scales)
yLims<- c(0,max(deathDat$age[which(deathDat$birt >= 1960)]))
xLims<- c(min(deathDat$X10000_loci),max(deathDat$X10000_loci))
plot(deathDat$X10000_loci[which(deathDat$sex == "F")],deathDat$age[which(deathDat$sex == "F")],col=alpha("gray",alpha=0.6),
     pch=16,cex=2,ylab="Age at death",xlab=expression(italic(""*H*"")["10,000 loci"]),cex.lab=1.5,ylim=yLims,xlim=xLims)
x1 <- min(deathDat$Froh_min1Mb)
x2 <- max(deathDat$Froh_min1Mb)
int <- mod$coefficients[1]
slope <- mod$coefficients[2]
sex <- mod$coefficients[3]

# line for females
lines(c(x1,x2),c(int+slope*x1,int+slope*x2),col="gray",lwd=3)
par(xpd=TRUE)
text(x=mean(c(par("usr")[1],par("usr")[2])),y=59,labels=expression("Females"),cex=1.7)
par(xpd=FALSE)
plot(deathDat$X10000_loci[which(deathDat$sex == "M")],deathDat$age[which(deathDat$sex == "M")],col=alpha("darkred",alpha=0.6),
     pch=16,cex=2,ylab="",xlab=expression(italic(""*H*"")["10,000 loci"]),cex.lab=1.5,ylim=yLims,xlim=xLims)

# line for males
lines(c(x1,x2),c(int+sex+slope*x1,int+sex+slope*x2),col="darkred",lwd=3)
par(xpd=TRUE)
text(x=mean(c(par("usr")[1],par("usr")[2])),y=59,labels=expression("Males"),cex=1.7)
par(xpd=FALSE)

# 10Mb 
mod <- lm(deathDat$age~ deathDat$allLoci + deathDat$sex)
summary(mod)
col<- rep("gray",nrow(deathDat))
col[which(deathDat$sex == "M")] <- "darkred"
library(scales)
yLims<- c(0,max(deathDat$age[which(deathDat$birt >= 1960)]))
xLims<- c(min(deathDat$allLoci[which(deathDat$birt >= 1960)]),max(deathDat$allLoci[which(deathDat$birt >= 1960)]))
plot(deathDat$allLoci[which(deathDat$sex == "F")],deathDat$age[which(deathDat$sex == "F")],col=alpha("gray",alpha=0.6),
     pch=16,cex=2,ylab="Age at death",xlab=expression(italic(""*H*"")["all loci"]),cex.lab=1.5,ylim=yLims,xlim=xLims)
x1 <- min(deathDat$allLoci)
x2 <- max(deathDat$allLoci)
int <- mod$coefficients[1]
slope <- mod$coefficients[2]
sex <- mod$coefficients[3]

# line for females
lines(c(x1,x2),c(int+slope*x1,int+slope*x2),col="gray",lwd=3)
plot(deathDat$allLoci[which(deathDat$sex == "M")],deathDat$age[which(deathDat$sex == "M")],col=alpha("darkred",alpha=0.6),
     pch=16,cex=2,ylab="",xlab=expression(italic(""*H*"")["all loci"]),cex.lab=1.5,ylim=yLims,xlim=xLims)

# line for males
lines(c(x1,x2),c(int+sex+slope*x1,int+sex+slope*x2),col="darkred",lwd=3)




# infer the total genomic inbreeding of the individual with the least heterozygosity in the population
maxF <- (max(het$allLoci)-min(het$allLoci))/max(het$allLoci)
print(maxF)



