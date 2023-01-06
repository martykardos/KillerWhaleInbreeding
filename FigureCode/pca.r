rm(list=ls())
library(ggplot2)
test <- read.table('kw_147.snp_QC_ld_pca.eigenvec1')
rownames(test) <- test[,1]
test <- test[,-1]
colnames(test) <- c("sample.id","PC1","PC2","PC3","Group")


################# with background #################
###
pdf('pc1_pc2_backgorund.pdf')
ggplot(test,aes(x=PC1,y=PC2,color=Group))+ geom_point() + labs(title = "Principal Components Analysis")+
  theme(plot.title = element_text(hjust = 0.5))+theme(axis.line = element_line(size=0.5, colour = "black"))
dev.off()
####
pdf('pc1_pc3_backgorund.pdf')
ggplot(test,aes(x=PC1,y=PC3,color=Group))+ geom_point() + labs(title = "Principal Components Analysis")+
  theme(plot.title = element_text(hjust = 0.5))+theme(axis.line = element_line(size=0.5, colour = "black")) 
dev.off()
###
pdf('pc2_pc3_backgorund.pdf')
ggplot(test,aes(x=PC2,y=PC3,color=Group))+ geom_point() + labs(title = "Principal Components Analysis")+
  theme(plot.title = element_text(hjust = 0.5))+theme(axis.line = element_line(size=0.5, colour = "black")) 
dev.off()
