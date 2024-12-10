# Name: plot_MDS
# Author: ELL
# Date: 8/30/2024
# Version:4.2.1
# Description: will read the data matrix from load_GC_data_and_sum_reps into
# DESeq2, read in metadata,pre-filter and plotMDS

library(readxl)
library(DESeq2)
library(dplyr)
library(Glimma)
library(sva)

#BiocManager::install("sva")
setwd('/Users/erikalesperance/Desktop/DEseq/gene_count_data')

# read in the data matrix
summed_counts<-readRDS("/Users/erikalesperance/Desktop/DEseq/collapsed_rep_dataset.Rdata")
dim(summed_counts)


summed_counts$replicates


samples=c("Stage_I_2A", "Stage_I_3A", "Stage_I_4A", "Stage_II_2A", "Stage_II_3A", "Stage_II_4A", 
          "Stage_III_2A", "Stage_III_3A", "Stage_III_4A", "Stage_IV_2A", "Stage_IV_3A", "Stage_IV_4A")

dev_stage=c("I","I","I", "II","II","II","III","III","III","IV","IV","IV")
tissue_type=c("inflo","inflo","inflo","inflo","inflo","inflo","inflo","inflo","inflo","inflo","inflo","inflo")
#dev_stage<-sub(".*([0-9]{2,2}D).*", "\\1",samples)


metadata<-data.frame(samples, tissue_type,dev_stage)

# create the factors of interest
metadata$dev_stage<-factor(metadata$dev_stage)
metadata$tissue_type<-factor(metadata$tissue_type)



# create the model
summed_counts<-DESeqDataSetFromMatrix(counts(summed_counts),colData = metadata, design=~0+dev_stage)

# pre-filter for reads where at least 3 samples have a count of 1 or higher
keep<-rowSums(counts(summed_counts)>=1)>=3
length(which(keep==1))
summed_counts_filt<-summed_counts[keep,]

# will load the plot...need to save within the html
glimmaMDS(summed_counts_filt, groups=metadata)
