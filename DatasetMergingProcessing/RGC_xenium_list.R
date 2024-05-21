library(tidyverse)
library(matrixStats) # This library is good for column/row wise operations like rowMeans()

setwd("/home/sam/scRNAseq/SCP509/expression/")
RGCexpMatrix <- read.csv("RGC_Atlas.csv", head = TRUE) #, sep="\t",skip=1, col.names = c("ID", "UMAP1", "UMAP2", "Cluster"))

setwd("/home/sam/scRNAseq/SCP509/cluster/")
RGCclusters <- read.csv("RGC_Atlas_coordinates.txt", head = TRUE, sep="\t",skip=1, col.names = c("ID", "UMAP1", "UMAP2", "Cluster", "Experiment","SampleID", "BatchID","BRGCkground"))

rownames(RGCexpMatrix) <- RGCexpMatrix$GENE 
RGCexpMatrix <- select(RGCexpMatrix, -GENE)
RGCclusters <- RGCclusters %>%
  select(ID, Cluster) %>%
  mutate(ID = str_replace_all(ID, "-", "."))

RGCexpMatrix_clean <- RGCexpMatrix %>%
  t() %>%
  data.frame() %>%
  rownames_to_column("ID") %>%
  left_join(RGCclusters, by = "ID") %>%
  column_to_rownames("ID")

RGC_clustMean <- RGCexpMatrix_clean %>%
  group_by(Cluster) %>%
  summarise_all(mean)

# Calculate the variance of eRGCh mean cluster column
variances <- apply(select(RGC_clustMean,-Cluster), 2, var)
# Find the mean variance of a given gene
mean_var <- mean(variances)
# Subset those genes with variance 10 times greater than the average
big_vars <- variances[variances > mean_var*10] 
RGC_interesting_genes <- names(big_vars)

RGC_clustMean_variant <- RGC_clustMean %>%
  ungroup() %>%
  column_to_rownames("Cluster") %>%
  select(all_of(RGC_interesting_genes))

# Get single cell matrix containing candidate genes
RGCexpMatrix_candidateGenes <- RGCexpMatrix_clean %>%
  select(Cluster, all_of(RGC_interesting_genes))
  

