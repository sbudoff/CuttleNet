library(tidyverse)

useful_genes <- read.csv('/home/sam/scRNAseq/Xenium/NeurIPS/AUCPAR_MAST_Retina_Genes.csv', header = F) %>%
  pull(V1)

# Load RGCs
expMatrix_RGC <- read.csv("/home/sam/scRNAseq/SCP509/expression/RGC_Atlas.csv", head = TRUE)
rownames(expMatrix_RGC) <- expMatrix_RGC$GENE
expMatrix_RGC <- expMatrix_RGC[useful_genes,]
expMatrix_RGC <- na.omit(expMatrix_RGC)

geneNames <- expMatrix_RGC$GENE  # Store gene names before removal
expMatrix_RGC$GENE <- NULL  # Remove the gene column from dataframe

expMatrix_RGC <- t(expMatrix_RGC)
colnames(expMatrix_RGC) <- geneNames

clusters <- read.csv("/home/sam/scRNAseq/SCP509/cluster/RGC_Atlas_coordinates.txt", head = TRUE, sep="\t",skip=1, col.names = c("ID", "UMAP1", "UMAP2", "Cluster", "Experiment","SampleID", "BatchID","Background"))
clusters$ID <- gsub("-", ".", clusters$ID)

RGC_df <- as.data.frame(expMatrix_RGC)
rownames(RGC_df) <- rownames(expMatrix_RGC)
RGC_df$Label <- clusters$Cluster[match(rownames(RGC_df), clusters$ID)]
RGC_df$Label <- as.factor(RGC_df$Label)


# Load BCs
expMatrix_BC <- read.csv("/home/sam/scRNAseq/SCP3/expression/exp_matrix.txt", head = TRUE, sep="\t")
rownames(expMatrix_BC) <- expMatrix_BC$GENE
expMatrix_BC <- expMatrix_BC[useful_genes,]
expMatrix_BC <- na.omit(expMatrix_BC)

geneNames <- expMatrix_BC$GENE  # Store gene names before removal
expMatrix_BC$GENE <- NULL  # Remove the gene column from dataframe

expMatrix_BC <- t(expMatrix_BC)
colnames(expMatrix_BC) <- geneNames

clusters <- read.csv("/home/sam/scRNAseq/SCP3/metadata/clust_retinal_bipolar.txt", head = TRUE, sep="\t",skip=1, col.names = c("ID", "Cluster", "subcluster"))
clusters$ID <- gsub("-", ".", clusters$ID)

BC_df <- as.data.frame(expMatrix_BC)
rownames(BC_df) <- rownames(expMatrix_BC)
BC_df$Label <- clusters$Cluster[match(rownames(BC_df), clusters$ID)]
BC_df$Label <- as.factor(BC_df$Label)

# Load ACs
expMatrix_AC <- read.csv("/home/sam/scRNAseq/MouseAC_gene_expression_matrix.csv", head = TRUE, sep="\t")
rownames(expMatrix_AC) <- expMatrix_AC$GENE
expMatrix_AC <- expMatrix_AC[useful_genes,]
expMatrix_AC <- na.omit(expMatrix_AC)

geneNames <- expMatrix_AC$GENE  # Store gene names before removal
expMatrix_AC$GENE <- NULL  # Remove the gene column from dataframe

expMatrix_AC <- t(expMatrix_AC)
colnames(expMatrix_AC) <- geneNames

clusters <- read.csv("/home/sam/scRNAseq/MouseAC_clusterfile.txt", head = TRUE, sep="\t",skip=1, col.names = c("ID", "UMAP1", "UMAP2", "Cluster"))
clusters$ID <- gsub("-", ".", clusters$ID)

AC_df <- as.data.frame(expMatrix_AC)
rownames(AC_df) <- rownames(expMatrix_AC)
AC_df$Label <- clusters$Cluster[match(rownames(AC_df), clusters$ID)]
AC_df$Label <- as.factor(AC_df$Label)

# Function to add missing columns from useful_genes to a dataframe and fill with zeros
add_missing_columns <- function(df, columns) {
  missing_cols <- setdiff(columns, names(df))
  df[missing_cols] <- 0  # Assign zero to all missing columns
  df <- df[, columns]  # Ensure consistent column order
  return(df)
}

# Apply this function to each dataframe
useful_genes <- c("Label", useful_genes)
RGC_df <- add_missing_columns(RGC_df, useful_genes)
AC_df <- add_missing_columns(AC_df, useful_genes)
BC_df <- add_missing_columns(BC_df, useful_genes)

# Combine the dataframes
combined_df <- bind_rows(RGC_df, AC_df, BC_df)
# Save the combined_df as an RData file
save(combined_df, file = '/home/sam/scRNAseq/Xenium/NeurIPS/AUCPRExpressionMats.RData')
