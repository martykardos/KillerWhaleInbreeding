library(data.table)
samps <- read.table("kw_151.snp.final.minGQ30.minGQ30_maxMiss20Perc_srkw_indvMissMax20Perc.chr9.depth.idepth",header=TRUE)[,1]
# make a vector of sample IDs. These will be the sequence IDs in the fasta file output
sampVec <- NULL
for(i in 1:length(samps)){
  sampVec <- c(sampVec,rep(paste(">",as.character(samps[i+1]),sep=""),2))
}
nInds <- length(samps)     # number of individuals
library(reshape2)
######################################################
# loop over chromosomes
######################################################
for(z in 21){
  vcf <- read.table(paste("kw_151.snp.final.minGQ30.minGQ30_maxMiss20Perc_srkw_indvMissMax20Perc_readDepCent95Perc_HWPValMin005_chr",z,".repeatMasked_phased.vcf",sep=""),header=FALSE)
  allNum <- nchar(as.character(vcf[,5]))
  if(sum(allNum > 1) > 0) vcf <- vcf[-which(allNum > 1),]
  genos1 <- matrix(NA,nrow=nInds,ncol=nrow(vcf))      # first haploid copy
  genos2 <- matrix(NA,nrow=nInds,ncol=nrow(vcf))      # second haploid copy
  print(ncol(genos1))
  print(nrow(genos1))
  
  for(i in 1:nInds){
    splitGenos <- colsplit(as.character(vcf[,i+9]),pattern="",names=c("all1","all2"))
    genos1[i,] <- splitGenos[,1]
    genos2[i,] <- colsplit(splitGenos[,2],pattern="",names=c("all1","all2"))[,2]
    print(i)
  }
  alleles <- vcf[,4:5]
  seqLeng <- max(as.numeric(vcf[,2]))
  snpLocs <- as.numeric(vcf[,2])
  starts <- seq(1,seqLeng,round(seqLeng/10))[1:10]
  ends <- (starts + round(seqLeng/10) -1)[1:10]
  
  #################################################
  # make the sequences in 10 batches
  #################################################
  for(j in 1:10){
    seqs <- rep(NA,nInds*4)   # vector for sequences from each haploid genome
    theseSNPLocs <- snpLocs[which(snpLocs >= starts[j] & snpLocs <= ends[j])]       # snp locations, genotypes, and alleles for loci in this batch
    theseGenos1 <- genos1[,which(snpLocs >= starts[j] & snpLocs <= ends[j])]
    theseGenos2 <- genos2[,which(snpLocs >= starts[j] & snpLocs <= ends[j])]
    theseAlleles <- alleles[which(snpLocs >= starts[j] & snpLocs <= ends[j]),]
    snpInds <- theseSNPLocs - starts[j] + 1                 # index the SNPs
    thisSeqLeng <- ends[j] - starts[j]      # length of the output sequence
    ticker <- 1   # iterate through seqs
    for(i in 1:nInds){
      thisSeq1 <- rep(NA,ncol(theseGenos1))
      thisSeq2 <- rep(NA,ncol(theseGenos2))
      thisSeq1 [which(theseGenos1[i,] == 0)] <- as.character(theseAlleles[which(theseGenos1[i,] == 0),1])
      thisSeq1 [which(theseGenos1[i,] == 1)] <- as.character(theseAlleles[which(theseGenos1[i,] == 1),2])
      thisSeq2 [which(theseGenos2[i,] == 0)] <- as.character(theseAlleles[which(theseGenos2[i,] == 0),1])
      thisSeq2 [which(theseGenos2[i,] == 1)] <- as.character(theseAlleles[which(theseGenos2[i,] == 1),2])
      outSeq1 <- rep("A",thisSeqLeng)
      outSeq1[snpInds] <- thisSeq1
      outSeq2 <- rep("A",thisSeqLeng)
      outSeq2[snpInds] <- thisSeq2
      seqs[ticker] <- sampVec[i]
      ticker <- ticker + 1
      seqs[ticker] <- paste0(outSeq1,collapse="")
      ticker <- ticker + 1
      seqs[ticker] <- sampVec[i]
      ticker <- ticker + 1
      seqs[ticker] <- paste0(outSeq2,collapse="")
      ticker <- ticker + 1
      print(i)
    }
    outSeqs <- NULL
    outSeqs <- cbind(outSeqs,seqs)
    write.table(outSeqs,file=paste("orcaChr",z,"FastEPRR_",j,"_repeatMassked.fas",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)
  }
}
