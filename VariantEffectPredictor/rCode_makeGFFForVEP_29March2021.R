###### configure the gff file for killer whales to work in VEP

setwd("~/Documents/orca/genomeAnnotations")

gffIn <- as.matrix(read.table("KW.gene.gff",header=FALSE))
badGenes <- which(gffIn[,3] == "CDS" & gffIn[,8] == ".")     # remove loci with no phase information
addRows <- 68
for(i in 1:(length(badGenes)-1)){
  if(badGenes[i] != badGenes[i+1]-1) addRows <- c(addRows,badGenes[i+1]-1)
}
remRows <- sort(c(addRows,badGenes))

gffIn <- gffIn[-remRows,]

geneNames <- paste("KW0000",1:9,sep="")
geneNames <- c(geneNames,paste("KW000",10:99,sep=""))
geneNames <- c(geneNames,paste("KW00",100:999,sep=""))
geneNames <- c(geneNames,paste("KW0",1000:9999,sep=""))
geneNames <- c(geneNames,paste("KW",10000:19710,sep=""))

outGff <- NULL

for(i in 1:length(geneNames)){
  if(length(grep(geneNames[i],gffIn[,9])) != 0){
    thisInGff <- gffIn[grep(geneNames[i],gffIn[,9]),]          # grab the annotations from this gene

    geneLine <- c(thisInGff[1,1:2],"gene",thisInGff[1,4:8])    # add a line with a gene identifier
    geneLine <- c(geneLine,paste("ID=gene:",geneNames[i],"_GENE;biotype=protein_coding",sep=""))

    mrnaLine <- geneLine
    mrnaLine [3] <- "transcript"
    mrnaLine[9] <- paste("ID=transcript:",geneNames[i],"_MRNA;Parent=",paste("gene:",geneNames[i],"_GENE",sep=""),";biotype=protein_coding",sep="")
    if(nrow(thisInGff) > 2){
      exons <- thisInGff[2:nrow(thisInGff),]
      exons[,3] <- rep("exon",nrow(exons))
      exons[,8] <- rep(".",nrow(exons))
      exonNames <- paste("exon_",1:nrow(exons),"_",geneNames[i],sep="")
      exons[,9] <- paste("Name=",exonNames,";Parent=transcript:",rep(geneNames[i],nrow(exons)),"_MRNA",sep="")

      cds <- thisInGff[2:nrow(thisInGff),]
      cds[,3] <- rep("CDS",nrow(cds))
      cds[,9] <- rep(paste("ID=CDS:",geneNames[i],"P;Parent=transcript:",geneNames[i],"_MRNA",sep=""),nrow(exons))
    }

    if(nrow(thisInGff) == 2){
      exons <- thisInGff[2,]
      exonNames <- paste("exon_1_",geneNames[i],sep="")
      exons[3] <- "exon"
      exons[8] <- "."
      exons[9] <- paste("Name=",exonNames,";Parent=transcript:",geneNames[i],"_MRNA",sep="")

      cds <- thisInGff[2,]
      cds[3] <- "CDS"
      cds[9] <- paste("ID=CDS:",geneNames[i],"P;Parent=transcript:",geneNames[i],"_MRNA",sep="")
    }

    thisOutGff <- rbind(geneLine,mrnaLine,exons,cds)
    outGff <- rbind(outGff,thisOutGff)
    print(i)
  }

}

outGff <- outGff[which(outGff[,1] %in% chroms),]
write.table(outGff,file="KW.gene.appended.gff",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
chroms <- paste("chr",1:22,sep="")







