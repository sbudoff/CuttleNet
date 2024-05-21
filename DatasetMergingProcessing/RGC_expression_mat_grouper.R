library(tidyverse)
library(matrixStats)

setwd("/home/sam/scRNAseq/SCP509/expression/")
expMatrix <- read.csv("RGC_Atlas.csv", head = TRUE) #, sep="\t",skip=1, col.names = c("ID", "UMAP1", "UMAP2", "Cluster"))

setwd("/home/sam/scRNAseq/SCP509/cluster/")
clusters <- read.csv("RGC_Atlas_coordinates.txt", head = TRUE, sep="\t",skip=1, col.names = c("ID", "UMAP1", "UMAP2", "Cluster", "Experiment","SampleID", "BatchID","Background"))


clustIDs <- unique(clusters$Cluster)
RGC_cluster_expression <- select(expMatrix, GENE)
rownames(RGC_cluster_expression) <- RGC_cluster_expression$GENE
RGC_cluster_median = data.frame(RGC_cluster_expression)
RGC_cluster_SD = data.frame(RGC_cluster_expression)
RGC_cluster_expressing = data.frame(RGC_cluster_expression)

for (clust in clustIDs) {

  RGC_type_cells <- clusters %>%
    select(ID, Cluster) %>%
    filter(Cluster == clust) %>%
    pull(ID)
  
  RGC_type_cells = gsub("-", ".", RGC_type_cells)
  
  RGCmat <- expMatrix %>%
    select(all_of(RGC_type_cells))
  
  RGC_cluster_expression[clust] <- rowMeans(RGCmat)
  RGC_cluster_SD[clust] <- rowSds(as.matrix(RGCmat))
  RGC_cluster_median[clust] <- rowMedians(as.matrix(RGCmat))
  
  RGC_cluster_expressing[clust] <- rowMeans(RGCmat != 0)
  
}

save(RGC_cluster_expression, RGC_cluster_expressing, RGC_cluster_SD, RGC_cluster_median, file="RGC_cluster_ExpressionMats.RData")


