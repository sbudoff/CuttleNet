library(tidyverse)
library(matrixStats) # This library is good for column/row wise operations like rowMeans()

setwd("/home/sam/scRNAseq/")
ACclusters <- read.csv("MouseAC_clusterfile.txt", head = TRUE, sep="\t",skip=1, col.names = c("ID", "UMAP1", "UMAP2", "Cluster"))
ACexpMatrix <- read.csv("MouseAC_gene_expression_matrix.csv", head = TRUE, sep="\t")

rownames(ACexpMatrix) <- ACexpMatrix$GENE 
ACexpMatrix <- select(ACexpMatrix, -GENE)
ACclusters <- ACclusters %>%
  select(ID, Cluster) %>%
  mutate(ID = str_replace_all(ID, "-", "."))

ACexpMatrix_clean <- ACexpMatrix %>%
  t() %>%
  data.frame() %>%
  rownames_to_column("ID") %>%
  left_join(ACclusters, by = "ID") %>%
  column_to_rownames("ID")

AC_clustMean <- ACexpMatrix_clean %>%
  group_by(Cluster) %>%
  summarise_all(mean)

# Calculate the variance of each mean cluster column
variances <- apply(select(AC_clustMean,-Cluster), 2, var)
# Find the mean variance of a given gene
mean_var <- mean(variances)
# Subset those genes with variance 10 times greater than the average
big_vars <- variances[variances > mean_var*10] 
AC_interesting_genes <- names(big_vars)

AC_clustMean_variant <- AC_clustMean %>%
  ungroup() %>%
  column_to_rownames("Cluster") %>%
  select(all_of(AC_interesting_genes))

# Get single cell matrix containing candidate genes
ACexpMatrix_candidateGenes <- ACexpMatrix_clean %>%
  select(Cluster, all_of(AC_interesting_genes))
  

