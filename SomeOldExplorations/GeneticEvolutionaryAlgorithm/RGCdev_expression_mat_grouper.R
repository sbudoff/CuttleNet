library(tidyverse)
library(matrixStats) # This library is good for column/row wise operations like rowMeans()

n_vals = 50
gene_subset = NA

setwd("/home/sam/scRNAseq/RGC_developement_sc/SCP1706/cluster/")
clusters <- read.csv("rgcP56_metadata.txt", head = TRUE, sep="\t",skip=1, col.names = c("ID", "UMAP.X","UMAP.Y","Cluster", "NGENE","NUMI"))
setwd("/home/sam/scRNAseq/RGC_developement_sc/SCP1706/expression/61dfa0bb771a5b0ddcd0ab78/")
# expMat_gene <- read.csv("features_P56.tsv", head = F, sep="\t")
# expMat_cells <- read.csv("barcodes_P56.tsv", head = F, sep="\t")
# expMatrix <- data.frame(matrix(0,ncol = length(expMat_cells), nrow = length(expMat_gene)))

# rownames(expMatrix) <- expMat_gene$V1
# colnames(expMatrix) <- expMat_cells$V1
# expMatrix_values <- Matrix::readMM("rgcP56_processed.mtx")
# library(adjHelpR)
# expMatrix_values <- as_tibble(expMatrix_values)
# 
# for (i in 1:nrow(expMatrix_values)) {
#   x = expMatrix_values$row[i]
#   y = expMatrix_values$column[i]
#   z = expMatrix_values$value[i]
#   expMatrix[x,y] = z
# }
# del(expMatrix_values, expMat_cells, expMat_cells)
# expMatrix$GENE <- expMat_gene
# save(expMatrix, file="RGC_P56_expMat.RData")
load("RGC_P56_expMat.RData")

clustIDs <- unique(clusters$Cluster)
cluster_expression <- select(expMatrix, GENE)
rownames(cluster_expression) <- cluster_expression$GENEV1
cluster_expression <- cluster_expression %>%
  mutate(GENE = GENE$V1)


ExpressionTableMaker <- function(expMatrix, gene_subset, n_vals = 100, max_val = 5.7) {
  # This function will create a frequency vector out of each genes binned into n_vals segments
  # Each gene's frequency vector will be gathered into the output df
  
  # Init the empty df organized by the genes being evaluated
  expMatrix <- na.omit(expMatrix)
  RGC_cluster_expression <- select(expMatrix, GENE)
  expMatrix <- t(select(expMatrix, -GENE))
  out <- as_tibble(t(select(RGC_cluster_expression,-GENE))) %>%
    add_row(Rbpms=1:n_vals) %>%
    mutate(Rbpms=NA)
  # Init the binning sequence
  min_val = min(min(expMatrix)) - 0.01
  if (is.na(max_val)) { max_val = max(max(expMatrix))}
  sequence <- seq.int(from=min_val,to=max_val,length.out = n_vals+1)
  
  # if gene_subset not given
  if (sum(is.na(gene_subset)) == 1) {
    gene_subset = RGC_cluster_expression$GENE
  } 
  
  # Iterate through each gene and store the identified frequency vector
  for (gene in gene_subset){
    out[,gene] <- as.numeric(table(cut(expMatrix[,gene],breaks=sequence)))
  }
  # Store the binning sequence
  out[,'value'] <- seq.int(from=min_val,to=max_val,length.out = n_vals)
  out
}

# Init the empty list to store frequency dfs
cluster_hist <- list()

# Init and begin populating the Ns vectors
cluster_hist[['N']]['Full'] <- length(expMatrix)-1
cluster_hist[['N ID']]['Full'] <- 'Full'
# Compute frequency dist for unclustered cell population
# RGC_cluster_hist[['Full']] <- ExpressionTableMaker(expMatrix, gene_subset, n_vals = n_vals)
# save(RGC_cluster_hist, file="RGC_cluster_ExpressionFrequencies_allGenes_FullOnly.RData")

# # Init dfs containing summary statistics for use in genetic algorithm
cluster_expression <- select(expMatrix, GENE)
cluster_median = data.frame(cluster_expression)
cluster_SD = data.frame(cluster_expression)
cluster_expressing = data.frame(cluster_expression)

# Iterate through all clusters and compute/store frequency dfs
for (clust in clustIDs) {
  # Get relevant sample Ids
  type_cells <- clusters %>%
    select(ID, Cluster) %>%
    filter(Cluster == clust) %>%
    pull(ID)
  # Convert sample ID to vector
  # type_cells = gsub("-", ".", type_cells)
  # Subset full dataset on samples in the given cluster
  mat <- expMatrix %>%
    select(all_of(type_cells))
  mat$GENE <- cluster_expression$GENE
  # Compute and the store the frequency table
  cluster_hist[[clust]] <- ExpressionTableMaker(mat, gene_subset, n_vals = n_vals)
  cluster_hist[['N']][clust] <- length(mat)-1
  cluster_hist[['N ID']][clust] <- clust
  # # Compute summary statistics 
  mat <- select(mat, -GENE)
  cluster_expression[clust] <- rowMeans(mat)
  cluster_SD[clust] <- rowSds(as.matrix(mat))
  cluster_median[clust] <- rowMedians(as.matrix(mat))
  cluster_expressing[clust] <- rowMeans(mat != 0)
  # Print to tell user all is good
  print(sprintf("%s completed, %s cells present to make frequency table from", clust, length((mat))))
}

cluster_hist[['Full']] <- cluster_hist[[clust]] * 0
for (clust in clustIDs) {
  # Add RGC
  cluster_hist[['Full']] <- cluster_hist[[clust]] + cluster_hist[['Full']]

  print(sprintf("%s added to full set", clust))
}

# Save the list of frequency dataframes for easy loading elsewhere
save(cluster_hist, file="RGC_P56__cluster_ExpressionFrequencies_allGenes.RData")
# Save summary stats
save(cluster_expression, cluster_expressing, cluster_SD, cluster_median, file="RGC_P56__cluster_ExpressionMats_allGenes.RData")

# Save the gene frequency tables in an excel sheet where each sheet is a cluster
library(openxlsx)

# load(file="AC_cluster_ExpressionFrequencies_allGenes.RData")

# create workbook
wb <- createWorkbook()

#Iterate through each list element to create a sheet based on that list (creating an anonymous function inside Map())
Map(function(data, nameofsheet){     
  addWorksheet(wb, nameofsheet)
  writeData(wb, nameofsheet, data)
}, cluster_hist, names(cluster_hist))

## Save workbook to excel file 
saveWorkbook(wb, file = "RGC_P56__cluster_hist_allGenes.xlsx", overwrite = TRUE)
