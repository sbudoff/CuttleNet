library(tidyverse)
library(matrixStats) # This library is good for column/row wise operations like rowMeans()


setwd("/home/sam/scRNAseq/SCP3/metadata/")
BPclusters <- read.csv("clust_retinal_bipolar.txt", head = TRUE, sep="\t",skip=1, col.names = c("ID", "Cluster", "subcluster"))

print(unique(BPclusters$subcluster))
setwd("/home/sam/scRNAseq/SCP3/expression/")
BPexpMatrix <- read.csv("exp_matrix.txt", head = TRUE, sep="\t")
rownames(BPexpMatrix) <- BPexpMatrix$GENE 
BPexpMatrix <- select(BPexpMatrix, -GENE)
BPclusters <- select(BPclusters, -subcluster)

BPexpMatrix_clean <- BPexpMatrix %>%
  t() %>%
  data.frame() %>%
  rownames_to_column("ID") %>%
  left_join(BPclusters, by = "ID") %>%
  filter(Cluster != "AC (Amacrine cell)",
         Cluster != "Doublets/Contaminants" ) %>%
  column_to_rownames("ID")

BP_clustMean <- BPexpMatrix_clean %>%
  group_by(Cluster) %>%
  summarise_all(mean)

# Calculate the variance of each mean cluster column
variances <- apply(select(BP_clustMean,-Cluster), 2, var)
# Find the mean variance of a given gene
mean_var <- mean(variances)
# Subset those genes with variance 10 times greater than the average
big_vars <- variances[variances > mean_var*10] 
BP_interesting_genes <- names(big_vars)

BP_clustMean_variant <- BP_clustMean %>%
  ungroup() %>%
  column_to_rownames("Cluster") %>%
  select(all_of(BP_interesting_genes))

# Get single cell matrix containing candidate genes
BPexpMatrix_candidateGenes <- BPexpMatrix_clean %>%
  select(Cluster, all_of(BP_interesting_genes))
  

