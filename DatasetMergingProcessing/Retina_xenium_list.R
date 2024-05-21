library(tidyverse)
library(matrixStats) # This library is good for column/row wise operations like rowMeans()


################################################################################
################################################################################
##########################     RGC     #########################################
################################################################################
################################################################################
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

################################################################################
################################################################################
##########################     BPs     #########################################
################################################################################
################################################################################

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

################################################################################
################################################################################
##########################     ACs     #########################################
################################################################################
################################################################################

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

################################################################################
################################################################################
##########################     Fusion     ######################################
################################################################################
################################################################################

setwd("/home/sam/scRNAseq/Xenium/")

Retina_interesting_genes <- unique(c(BP_interesting_genes, AC_interesting_genes, RGC_interesting_genes, "Rbpms"))

# Get single cell matrix containing candidate genes
ACexpMatrix_candidateGenes <- ACexpMatrix_clean %>%
  select(Cluster, any_of(Retina_interesting_genes))

# Get single cell matrix containing candidate genes
BPexpMatrix_candidateGenes <- BPexpMatrix_clean %>%
  select(Cluster, any_of(Retina_interesting_genes))

# Get single cell matrix containing candidate genes
RGCexpMatrix_candidateGenes <- RGCexpMatrix_clean %>%
  select(Cluster, any_of(Retina_interesting_genes))

# Get the unique column names from all data frames
all_genes <- unique(c(colnames(ACexpMatrix_candidateGenes), colnames(BPexpMatrix_candidateGenes), colnames(RGCexpMatrix_candidateGenes)))

# Ensure all data frames have the same columns in the same order
ACexpMatrix_candidateGenes <- ACexpMatrix_candidateGenes[, all_genes]
BPexpMatrix_candidateGenes <- BPexpMatrix_candidateGenes[, all_genes]
RGCexpMatrix_candidateGenes <- RGCexpMatrix_candidateGenes[, all_genes]

# Combine data frames using rbind
Retina_expMatrix_candidateGenes <- bind_rows(BPexpMatrix_candidateGenes, ACexpMatrix_candidateGenes, RGCexpMatrix_candidateGenes)

Retina_clustMean <- Retina_expMatrix_candidateGenes %>%
  group_by(Cluster) %>%
  summarise_all(mean)

save(Retina_expMatrix_candidateGenes, Retina_clustMean, file="Retina_expMatrix_candidateGenes.RData")

Retina_expMatrix_candidateGenes <- replace(Retina_expMatrix_candidateGenes, is.na(Retina_expMatrix_candidateGenes), 0)

save(Retina_expMatrix_candidateGenes, file="Retina_expMatrix_clean.RData")
