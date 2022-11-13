# exclude loci in repeat regions from data for recombination rate estimates with FastEPRR
setwd("/Users/martin.kardos/Documents/orca/orca/fastEPRR_estimates/filterForRepeatMasking")
repeats <- read.table("repeatMask_kw_ref.bed")



for(i in 1:22){       # loop through chromosomes
  theseSNPs <- read.table(paste("kw_151.snp.final.minGQ30.minGQ30_maxMiss20Perc_srkw_indvMissMax20Perc_readDepCent95Perc_HWPValMin005_chr",i,"_012.012.pos",sep=""),header=FALSE) 
  theseExclude <- rep(0,nrow(theseSNPs))
  theseRepeats <- repeats[repeats[,1] == paste("chr",i,sep=""),]
  for(j in 1:nrow(theseSNPs)){
    if(sum(theseSNPs[j,2] >= theseRepeats[,2] & theseSNPs[j,2] <= theseRepeats[,3]) > 0) theseExclude[j] <- 1
    #print(j)
  }
  
  if(is.null(theseExclude) == FALSE){
    outExclude <- theseSNPs[theseExclude==1,]
    write.table(outExclude,file=paste("excludeLociRepeats_chr",i,sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
  }
  print(i)
}