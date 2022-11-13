bqs <- NULL
bqs <- cbind(bqs,rep("KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK",242090369))
write.table(bqs,file="ASM322739v1_HiC_pruned_reads1_BQs",quote=FALSE,row.names=FALSE,col.names=FALSE)

rm(bqs)
gc()

pluses <- NULL      
pluses <- cbind(pluses,rep("+",242090369))
write.table(pluses,file="ASM322739v1_HiC_pruned_reads1_pluses",quote=FALSE,row.names=FALSE,col.names=FALSE)

