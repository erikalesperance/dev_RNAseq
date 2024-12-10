# Name: load GC data and sum reps
# Author: ELL
# Date: 8/30/2024
# Description: will load in the parsed gene count data from easley 

# set wd
#install.packages("DESeq2")
library(DESeq2)
R.version

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("DESeq2")
# Description: will load in the parsed gene count data and collapse the replicates
# This is done using DESeq
# This is a lot of data that needs to be read in, so 100GB allocated 
# set cwd
setwd('/Users/erikalesperance/Desktop/DEseq/gene_count_data')

library(DESeq2)
# put all of the gene_count files into a character vector in order to load them 
gene_count_files<-dir(pattern="*\\.tab$")

# read the data into a dataframe
data_table<-data.frame(sampleName=gene_count_files,fileName=gene_count_files)
# create a column that we want to collapse the replicates on (i.e. just the plant #)
remove_names<-sub("-.*",'',data_table$sampleName)
data_table$replicates<-remove_names

# load the DESeq data set
# we read in the data using the DESeqDataSetFromHTSeqCount bc our dataset is in the same format as HTSeq data
dds_set<-DESeqDataSetFromHTSeqCount(sampleTable=data_table,directory='/Users/erikalesperance/Desktop/DEseq/gene_count_data/',design=~0)

# collapse the replicates 
collapsed_reps<-collapseReplicates(dds_set,dds_set$replicates)
collapsed_reps
# save the data set 
saveRDS(collapsed_reps,file="/Users/erikalesperance/Desktop/DEseq/collapsed_rep_dataset.Rdata")
