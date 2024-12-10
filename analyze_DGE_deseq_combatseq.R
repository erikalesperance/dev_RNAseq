# Name: analyze DGE deseq Bidens dev stages 
# Author: EY (based off of code written by E. Dittmar)
# Date: 08/30/2024
# Description: Will analyze the output from the DESeq DGE with combatseq for pairwise dev stages
# need Functions.R written by ED in the working directory for this to run properly 


setwd('/Users/erikalesperance/Desktop/DEseq/gene_count_data')

library(dplyr)
library(ggplot2)
library(UpSetR)
library(Glimma)
source("~/Desktop/Computational/dev_RNAseq-main/Functions.R")

# now analyze 
# read in the data, CHANGE PATH ERIKA To where results are on computer, results in a folder on their own, all csv in folder will
##read##
DEData_pairwise_cs<-ImportCSVs('deseq_results',0.05)
# filter out significant results
mydataSig_pairwise_cs<-lapply(DEData_pairwise_cs,SigDEdf,PvaluesCol=7,CritP=0.05)


# see which genes overlap (input into upset plot)
SigOverlap_pairwise_cs<-GeneSets(mydataSig_pairwise_cs$result_I_v_II[1], mydataSig_pairwise_cs$result_II_v_III[1], mydataSig_pairwise_cs$result_III_v_IV[1])
names(SigOverlap_pairwise_cs)
lapply(SigOverlap_pairwise_cs,function(x) {length(x$Gene)})
SigOverlapGraph_pairwise_cs<-lapply(mydataSig_pairwise_cs, function(x) {x$Gene})

# create an upset plot of DE expression by pairwise dev_stage
png("sequential_pairwise_upset.png", res=215, width = 1800, height=1000)
upset(fromList(SigOverlapGraph_pairwise_cs),order.by="freq",nsets=13,nintersects=20, text.scale = 1.5)
dev.off()

#save as csv
write.csv(SigOverlap_pairwise_cs[["mydataSig_pairwise_cs$result_I_v_II[1]Only"]], file= 'pairwise_I_v_II.csv')
write.csv(SigOverlap_pairwise_cs[["mydataSig_pairwise_cs$result_II_v_III[1]Only"]], file= 'pairwise_II_v_III.csv')
write.csv(SigOverlap_pairwise_cs[["mydataSig_pairwise_cs$result_III_v_IV[1]Only"]], file= 'pairwise_III_v_IV.csv')
write.csv(SigOverlap_pairwise_cs[["InCommonAll"]], file= 'InCommonAll.csv')
