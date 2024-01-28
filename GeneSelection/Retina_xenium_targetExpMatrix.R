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

# Get the unique column names from all data frames
all_genes <- unique(c(colnames(ACexpMatrix_candidateGenes), colnames(BPexpMatrix_candidateGenes), colnames(RGCexpMatrix_candidateGenes)))

# Check which genes from retina_genes are not in all_genes
genes_not_in_all <- Retina_interesting_genes$gene[!Retina_interesting_genes$gene %in% all_genes]

# Print the genes that are not in all_genes
if (length(genes_not_in_all) > 0) {
  print(genes_not_in_all)
}

# Add missing genes to each dataframe
# Identify missing columns
missing_columns <- setdiff(all_genes, colnames(ACexpMatrix_candidateGenes))
# Add missing columns and set their values to 0
ACexpMatrix_candidateGenes[, missing_columns] <- 0
# Identify missing columns
missing_columns <- setdiff(all_genes, colnames(BPexpMatrix_candidateGenes))
# Add missing columns and set their values to 0
BPexpMatrix_candidateGenes[, missing_columns] <- 0
# Identify missing columns
missing_columns <- setdiff(all_genes, colnames(RGCexpMatrix_candidateGenes))
# Add missing columns and set their values to 0
RGCexpMatrix_candidateGenes[, missing_columns] <- 0


# Ensure all data frames have the same columns in the same order
ACexpMatrix_candidateGenes <- ACexpMatrix_candidateGenes[, c("Dataset", "Cluster", Retina_interesting_genes$gene)]
BPexpMatrix_candidateGenes <- BPexpMatrix_candidateGenes[,  c("Dataset", "Cluster", Retina_interesting_genes$gene)]
RGCexpMatrix_candidateGenes <- RGCexpMatrix_candidateGenes[,  c("Dataset", "Cluster", Retina_interesting_genes$gene)]

# Combine data frames using rbind
Retina_expMatrix_candidateGenes <- bind_rows(BPexpMatrix_candidateGenes, ACexpMatrix_candidateGenes, RGCexpMatrix_candidateGenes)

Retina_expMatrix_candidateGenes <- replace(Retina_expMatrix_candidateGenes, is.na(Retina_expMatrix_candidateGenes), 0)

save(Retina_expMatrix_candidateGenes, file="Retina_expMatrix_clean_final2.RData")
