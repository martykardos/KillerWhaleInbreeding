################################################################################
# plot population distributions of Froh first
################################################################################
setwd("/Users/martin.kardos/Documents/orca/orca")
setwd("~/Documents/orca/analyses_31May2022/ROH3")
Froh <- read.table("Froh_7August2022",header=TRUE)
popKey <- read.table("popKey.txt",header=TRUE)
popVec <- popKey[match(Froh[,1],popKey[,2]),1]
Froh <- Froh[which(popVec == "SRKW"),]
Froh <- Froh[order(Froh[,3]),]
Fp <- rep(0,nrow(Froh))
Fp[which (Froh[,1] == "J52")] <- 0.125
Fp[which (Froh[,1] == "J51")] <- 0.0625
Fp[which (Froh[,1] == "Sooke_calf")] <- 0.0625
Fp[which (Froh[,1] == "J42")] <- 0.25
par(mar=c(4,6,2,2))
plot(c(0,100),c(0,0.5),ylab=expression(italic(""*F*"")["ROH,1Mb"]),xlab="",type="n",cex.axis=1.2,cex.lab=1.5)
cols <- rep("lightblue",100)
cols[which(Fp == 0.0625)] <- "#fee0d2"
cols[which(Fp == 0.125)] <- "#fc9272"
cols[which(Fp == 0.25)] <- "#de2d26"
for(i in 1:100){
  rect(xleft=i-1,xright=i,ybottom=0,ytop=Froh[i,3],col=cols[i])  
}
legend(x=0,y=0.35,pch=15,col=c("lightblue","#fee0d2","#fc9272","#de2d26"),xjust=FALSE,yjust=FALSE,
       legend = c(
         expression(paste(italic(""*F*"")[P]," = 0",sep="")),
         expression(paste(italic(""*F*"")[P]," = 0.0625",sep="")),
         expression(paste(italic(""*F*"")[P]," = 0.12",sep="")),
         expression(paste(italic(""*F*"")[P]," = 0.25",sep=""))
       ),border="NULL")