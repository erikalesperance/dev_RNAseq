# Name: run_DGE_combatseq.R
# Author: ELL, adapted from EY
# Date: 08/30/2024
# Description: differential gene expression on the infloresence developmental stages accounting for variation with combatseq

library(dplyr)
library(DESeq2)
library(Glimma)
library(sva)
library(edgeR)
source("~/Desktop/Computational/dev_RNAseq-main/Functions.R")



# read in and process data
setwd('/Users/erikalesperance/Desktop/DEseq/gene_count_data')

# read in the data matrix
summed_counts<-readRDS("/Users/erikalesperance/Desktop/DEseq/collapsed_rep_dataset.Rdata")

samples=c("Stage_I_2A", "Stage_I_3A", "Stage_I_4A", "Stage_II_2A", "Stage_II_3A", "Stage_II_4A", 
           "Stage_III_2A", "Stage_III_3A", "Stage_III_4A", "Stage_IV_2A", "Stage_IV_3A", "Stage_IV_4A")

dev_stage=c("I","I","I","II","II","II","III","III","III","IV","IV","IV")
dev_stage <- as.factor(dev_stage)
metadata<-data.frame(samples, dev_stage)


# create the factors of interest
metadata$dev_stage<-factor(metadata$dev_stage)

# create the model (wrt dev_stage)
summed_counts<-DESeqDataSetFromMatrix(counts(summed_counts),colData = metadata, design=~0+dev_stage)


# pre-filter for reads where at least 3 samples have a count of 1 or higher
keep<-rowSums(counts(summed_counts)>=1)>=3
length(which(keep==1))
summed_counts_filt<-summed_counts[keep,]

# get a dataframe of the counts
count_matrix <- as.matrix(counts(summed_counts_filt))

# sort by batch (sample group) and group (dev_stage) WITHOUT yl, sam, other 4 stages 
batch <- c(2,3,4,2,3,4,2,3,4,2,3,4)
group <- c(1,1,1,2,2,2,3,3,3,4,4,4)

# adjust counts using combat seq...output is a matrix of adjusted counts 
adjusted_counts <- ComBat_seq(count_matrix, batch=batch, group=group)

# write to a CSV file
write.csv(as.data.frame((adjusted_counts)), file='adjusted_counts_combatseq.csv')


# plot MDS of adjusted counts...can compare this with previous plot pre-combat seq
glimmaMDS(adjusted_counts, group=metadata)

#convert matrix to data frame
adjusted_counts_dataframe <- as.data.frame(adjusted_counts)

# plot PCA of adjusted counts
png("pca_combseq.png", res=215, width = 1200, height=1000)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plotPCA(adjusted_counts_dataframe,labels=TRUE, col=dev_stage)
legend("topright",inset=c(-0.4,0),legend=unique(dev_stage), fill=dev_stage)
dev.off()

library(ggplot2)
# Transpose the data so that rows are samples and columns are variables
pca <- prcomp(t(adjusted_counts_dataframe), scale. = TRUE)
# Create a data frame for plotting
pca_data <- as.data.frame(pca$x)
pca_data$dev_stage <- dev_stage
# Plot the first two principal components
p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = dev_stage)) +
  geom_point() +
  labs(title = "PCA Plot", x = "PC1", y = "PC2")


# create the model with adjusted counts
adjusted_counts_deseq<-DESeqDataSetFromMatrix((adjusted_counts),colData = metadata, design=~0+dev_stage)

# run the DGE analysis
DESeq_dataset_results_combatseq<-DESeq(adjusted_counts_deseq,parallel=TRUE)
# save the output Erika's pathway
saveRDS(DESeq_dataset_results_combatseq, file='deseq_results/deseq_dataset_results_combatseq.RData')

# set up the pairwise contrasts (IvII, IIvIII, IIIvIV, IvIV) play w alpha p values and see how data changes##
result_I_v_II_combatseq<-results(DESeq_dataset_results_combatseq,contrast=c("dev_stage","I","II"),alpha=0.05,parallel=TRUE)
result_II_v_III_combatseq<-results(DESeq_dataset_results_combatseq,contrast=c("dev_stage","II","III"),alpha=0.05,parallel=TRUE)
result_III_v_IV_combatseq<-results(DESeq_dataset_results_combatseq,contrast=c("dev_stage","III","IV"),alpha=0.05,parallel=TRUE)
result_I_v_IV_combatseq<-results(DESeq_dataset_results_combatseq,contrast=c("dev_stage","I","IV"),alpha=0.05,parallel=TRUE)
result_I_v_III_combatseq<-results(DESeq_dataset_results_combatseq,contrast=c("dev_stage","I","III"),alpha=0.05,parallel=TRUE)
result_II_v_IV_combatseq<-results(DESeq_dataset_results_combatseq,contrast=c("dev_stage","II","IV"),alpha=0.05,parallel=TRUE)



# write to CSV file
write.csv(as.data.frame(result_I_v_II_combatseq), file='deseq_results/result_I_v_II.csv')
write.csv(as.data.frame(result_II_v_III_combatseq), file='deseq_results/result_II_v_III.csv')
write.csv(as.data.frame(result_III_v_IV_combatseq), file='deseq_results/result_III_v_IV.csv')
write.csv(as.data.frame(result_I_v_IV_combatseq), file='deseq_results/result_I_v_IV.csv')
write.csv(as.data.frame(result_I_v_III_combatseq), file='deseq_results/result_I_v_III.csv')
write.csv(as.data.frame(result_II_v_IV_combatseq), file='deseq_results/result_II_v_IV.csv')

