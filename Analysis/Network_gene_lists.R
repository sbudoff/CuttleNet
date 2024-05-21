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

rm(RGCexpMatrix, RGCclusters)
gc()

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

rm(BPexpMatrix, BPclusters)
gc()

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

rm(ACexpMatrix, ACclusters)
gc()

################################################################################
################################################################################
##########################     Fusion     ######################################
################################################################################
################################################################################

setwd("/home/sam/scRNAseq/Xenium/")

Retina_interesting_genes <- read.csv('Xenium Gene List_final2.csv') %>%
  mutate_all(~str_replace_all(., " ", "")) %>%  # Remove spaces in all columns
  mutate(gene = str_replace_all(gene, "-", "."))

# Get single cell matrix containing candidate genes
ACexpMatrix_candidateGenes <- ACexpMatrix_clean %>%
  select(Cluster, any_of(Retina_interesting_genes$gene)) %>%
  mutate(Dataset = 'Yan2020')

# Get single cell matrix containing candidate genes
BPexpMatrix_candidateGenes <- BPexpMatrix_clean %>%
  select(Cluster, any_of(Retina_interesting_genes$gene)) %>%
  mutate(Dataset = 'Karthik2016')

# Get single cell matrix containing candidate genes
RGCexpMatrix_candidateGenes <- RGCexpMatrix_clean %>%
  select(Cluster, any_of(Retina_interesting_genes$gene)) %>%
  mutate(Dataset = 'Tran2019')

# Identify genes not used for a respective network
RGC_missing <- setdiff(Retina_interesting_genes$gene, colnames(RGCexpMatrix_candidateGenes))
BC_missing <- setdiff(Retina_interesting_genes$gene, colnames(BPexpMatrix_candidateGenes))
AC_missing <- setdiff(Retina_interesting_genes$gene, colnames(ACexpMatrix_candidateGenes))
Classes_missing <- Retina_interesting_genes$gene[Retina_interesting_genes$gene %in% c(RGC_missing, BC_missing, AC_missing)]

# Save a dataframe with the genes to be used in each network
Network_genes_df <- Retina_interesting_genes %>%
  mutate(Classes = ifelse(gene %in% Classes_missing, 0, 1),
         RGC = ifelse(gene %in% RGC_missing, 0, 1),
         BC = ifelse(gene %in% BC_missing, 0, 1),
         AC = ifelse(gene %in% AC_missing, 0, 1))

save(Network_genes_df, file="Network_genes_df.RData")

# Extract the gene names from the columns of 'Retina_expMatrix_candidateGenes'
load('Retina_expMatrix_clean_final2.RData')
genes_order <- colnames(Retina_expMatrix_candidateGenes)
genes_order <- genes_order[genes_order != c('Dataset', 'Cluster')]

# Confirm gene order
print(all(Network_genes_df$gene == genes_order))

Network_genes_df <- Network_genes_df %>%
  mutate(Python_Index = seq(0, nrow(Network_genes_df)-1))

Class_indices <- Network_genes_df %>%
  filter(Classes == 1) %>%
  pull(Python_Index)

RGC_indices <- Network_genes_df %>%
  filter(RGC == 1) %>%
  pull(Python_Index)

AC_indices <- Network_genes_df %>%
  filter(AC == 1) %>%
  pull(Python_Index)

BC_indices <- Network_genes_df %>%
  filter(BC == 1) %>%
  pull(Python_Index)

# Save teh expression matrix with the network specific indices list
save(Retina_expMatrix_candidateGenes, Network_genes_df, Class_indices, RGC_indices, AC_indices, BC_indices, file="Network_genes.RData")

# Save 'Network_genes_df' as a CSV file
write.csv(Network_genes_df, "Network_genes_df.csv", row.names = FALSE)
