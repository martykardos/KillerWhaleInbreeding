# make a histogram of read depth

setwd("~/Documents/orca/analyses_31May2022/ROH3")
library(data.table)
hardy <- as.matrix(fread("kw_151.snp.final.passOnly.biallelicOnly.noIndels.mac1_indivFiltered.hardy.hwe",header=TRUE))
depth <- as.matrix(fread("kw_151.snp.final.passOnly.biallelicOnly.noIndels.mac1_indivFiltered.siteDepth.ldepth",header=TRUE))

library(reshape2)
hets <- colsplit(hardy[,3],"/",names=c("hom1","het1","hom2"))[,2]
expHets <- colsplit(hardy[,4],"/",names=c("hom1","het1","hom2"))[,2]
excessHets <- which(hets > expHets & as.numeric(hardy[,6] < 0.01))
hetExPVal <- as.numeric(hardy[,8])
hetDefPVal <- as.numeric(hardy[,7])

depthVec <- as.numeric(depth[,3])/148
hist(depthVec[which(depthVec < 50)],breaks=100)
lines(c(17,17),c(0,800000))

# associate heterozygosity statistics with read depth
hardyLoci <- paste(hardy[,1],"_",hardy[,2],sep="")
depthLoci <- paste(depth[,1],"_",depth[,2],sep="")
hardyMatch <- match(depthLoci,hardyLoci)

outStats <- cbind(depth[,1],depth[,2],depthVec,hets[hardyMatch],expHets[hardyMatch],hetExPVal[hardyMatch],hetDefPVal[hardyMatch])
colnames(outStats) <- c("chrom","position","meanDepth","obsHet","expHet","hetExcessP","hetDefP")

par(mfrow=c(2,1))
hist(as.numeric(outStats[which(as.numeric(outStats[,6]) <0.05),3]))
median(as.numeric(outStats[which(as.numeric(outStats[,6]) <0.05),3]))
median(as.numeric(outStats[which(as.numeric(outStats[,6]) >0.05),3]))

removeLoci <- which(as.numeric(outStats[,3]) > 17 | as.numeric(outStats[,6]) < 0.01 | as.numeric(outStats[,3]) < 5)  # identify loci with high read depth OR an excess of heterozygotes to remove from analysis

keepPos <- trimws(outStats[-removeLoci,2])
keepChr <- outStats[-removeLoci,1]

outKeep <- cbind(keepChr,keepPos)
hist(as.numeric(outStats[which(as.numeric(outStats[,3]) < 20),3]),xlab="Read depth",ylab="Frequency",breaks = 50,main="")
lines(c(5,5),c(0,10e5),lty="dashed")
lines(c(17,17),c(0,10e5),lty="dashed")