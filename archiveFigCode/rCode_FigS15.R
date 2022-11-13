par(mfrow=c(2,1),mar=c(5,5,2,2),xpd=TRUE)


# shape = 02, scale= 0.1
sDist <- rgamma(100000,shape=0.2,scale=0.1)
sDist[1:2000] <- sample(seq(0.96,1,0.01),2000,replace=TRUE)

starts <- seq(0,0.95,0.05)
counts <- rep(NA,20)
for(i in 1:length(starts)){
  counts[i] <- sum(sDist > starts[i] & sDist <= (starts[i] + 0.05))
}
props <- counts/100000
h <- 0.5*(exp(-13*seq(0,1,0.0001)))

plot(c(-1,0),c(0,0.9),xlab=expression(paste("-",italic(""*s*""))),ylab="Frequency",type="n",cex.lab=1.5)
for(i in 1:length(starts)){
  rect(xright=0-(starts[i]),xleft=0-(starts[i]+0.05),ybottom=0,ytop=props[i],col="#a6bddb",border="#a6bddb")
}
lines(rev(-seq(0,1,0.0001)),rev(h),lwd=2,col="black")

legend(x=-1,y=0.8,legend="shape = 0.2, scale=0.1",bty="n")




# shape = 0.2, scale= 0.2
sDist <- rgamma(100000,shape=0.4,scale=0.2)
sDist[1:2000] <- sample(seq(0.96,1,0.01),2000,replace=TRUE)

starts <- seq(0,0.95,0.05)
counts <- rep(NA,20)
for(i in 1:length(starts)){
  counts[i] <- sum(sDist > starts[i] & sDist <= (starts[i] + 0.05))
}
props <- counts/100000
h <- 0.5*(exp(-13*seq(0,1,0.0001)))

plot(c(-1,0),c(0,0.6),xlab=expression(paste("-",italic(""*s*""))),ylab="Frequency",type="n",cex.lab=1.5)
text(x=-1.25,y=0.65,labels="B",cex=2)
for(i in 1:length(starts)){
  rect(xright=0-(starts[i]),xleft=0-(starts[i]+0.05),ybottom=0,ytop=props[i],col="#a6bddb",border="#a6bddb")
}
lines(rev(-seq(0,1,0.0001)),rev(h),lwd=2,col="black")

legend(x=-1,y=0.6,legend="shape = 0.2, scale=0.2",bty="n")
