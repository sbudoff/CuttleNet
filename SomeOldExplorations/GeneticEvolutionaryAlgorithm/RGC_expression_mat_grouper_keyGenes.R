library(tidyverse)
library(matrixStats) # This library is good for column/row wise operations like rowMeans()

n_vals = 10

gene_subset <- c('Rbpms', 'Opn1sw','Grm6',
                 'Spp1','Tbr1','Pcdh20','Opn4','Il1rapl2','Adcyap1','Mmp17','Col25a1','Gpr88',
                 'Igf1', 'Igfbp2', 'Ppp1r1c', 'Parm1', 'Neurod2', 'Meis2', 'S100b', 
                 'Spon1', 'Chl1', 'Prkcq', 'Foxp2', 'Ctxn3', 'Irx4', 'Penk', 'Gal', 'Fxyd6',  
                 'Zeb2', 'Tfap2d', 'Chrna3', 'RP23-407N2.2', 'Dkk3', 'Cck', 'Trp53i11', 'Irx3',
                 'Cdkn1c', 'Pcdh10', 'Khdrbs2', 'Necab2', 'Ndnf', 'Isl2', 'Pcdh11x', 'Crhbp', 
                 'Pcdh9', 'Kcnip4', 'Bex1', 'Rprm', 'Ppp1r17','Evc2', 'Synpr', 'Mt1', 'Pcp4l1',
                 'Ramp3', 'Lpl','Lypd6','Syt2', 'Lmo2', 'Coch','Eomes', 'Irx6', 'Cd83', 'Gabra2', 
                 'Reln', 'Slc6a1', 'Sorl1', 'Kcnb2', 'Kcnab3', 'Marcksl1', 'Gprc5b', 'Six6',
                 'Sh3bgr', 'Atp2b4', 'Fgf1', 'Pou6f2', 'Kitl', 'Clstn2', 'Kcna1', 'Kcnd2', 
                 'Lmo1', 'Sv2c', 'Pvalb', 'Ntng1', 'Lmo4', 'Kcnip1', 'Isl1', 
                 'Dlgap1', 'Cpne4', 'Calb1', 'BC048546', 'Zmat4', 'Stk32a',
                 'Gnai1', 'Dlg2', 'Cplx2', 'Vamp1', 'Slc24a2', 'Satb1', 'Pcsk2', 'Gabrg3',
                 'Fgf13', 'Cacng3', 'Alcam', 'Prph', 'Nrxn3', 'Mgat4c', 'Bdnf', 'Mafb',
                 'Pou4f1', 'Junb', 'Magi1', 'Syndig1l', 'Stk32c', 'Car10', '6330403K07Rik',
                 'Cartpt', 'Zfhx3', 'Vgf', 'Nrgn', 'Lgals1', 'Ebf3', 'Diras2', 'Ano3',
                 '2510009E07Rik', 'Sema5a', 'Myo1b', 'Ly6h', 'Sptssb', 'Ndrg2', 'Lxn',
                 'Chrnb3')

setwd("/home/sam/scRNAseq/SCP509/expression/")
# This code block that is commented out reads in the full dataset=
expMatrix <- read.csv("RGC_Atlas.csv", head = TRUE) #, sep="\t",skip=1, col.names = c("ID", "UMAP1", "UMAP2", "Cluster"))
rownames(expMatrix) <- expMatrix$GENE
gene_subset = NA
# # The below line will subset the ful raw dataset by the key gene list, if this is not desired omit it
# expMatrix <- expMatrix[gene_subset,]
# expMatrix <- na.omit(expMatrix)
# save(expMatrix, file="RGC_cluster_ExpressionMats_keyGenes.RData")

# # If raw data has already been partitioned by key genes run the below
# load("RGC_cluster_ExpressionMats_keyGenes.RData")

# Open and find the cluster specific sample IDs
setwd("/home/sam/scRNAseq/SCP509/cluster/")
clusters <- read.csv("RGC_Atlas_coordinates.txt", head = TRUE, sep="\t",skip=1, col.names = c("ID", "UMAP1", "UMAP2", "Cluster", "Experiment","SampleID", "BatchID","Background"))
clustIDs <- unique(clusters$Cluster)
RGC_cluster_expression <- select(expMatrix, GENE)
rownames(RGC_cluster_expression) <- RGC_cluster_expression$GENE


ExpressionTableMaker <- function(expMatrix, gene_subset, n_vals = 100) {
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
  max_val = max(max(expMatrix))
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
RGC_cluster_hist <- list()

# Init and begin populating the Ns vectors
RGC_cluster_hist[['N']]['Full'] <- length(expMatrix)-1
RGC_cluster_hist[['N ID']]['Full'] <- 'Full'
# Compute frequency dist for unclustered cell population
RGC_cluster_hist[['Full']] <- ExpressionTableMaker(expMatrix, gene_subset, n_vals = n_vals)

# Init dfs containing summary statistics for use in genetic algorithm
RGC_cluster_expression <- select(expMatrix, GENE)
RGC_cluster_median = data.frame(RGC_cluster_expression)
RGC_cluster_SD = data.frame(RGC_cluster_expression)
RGC_cluster_expressing = data.frame(RGC_cluster_expression)

# Iterate through all clusters and compute/store frequency dfs
for (clust in clustIDs) {
  # Get relevant sample Ids
  RGC_type_cells <- clusters %>%
    select(ID, Cluster) %>%
    filter(Cluster == clust) %>%
    pull(ID)
  # Convert sample ID to vector
  RGC_type_cells = gsub("-", ".", RGC_type_cells)
  # Subset full dataset on samples in the given cluster
  RGCmat <- expMatrix %>%
    select(all_of(RGC_type_cells))
  RGCmat$GENE <- RGC_cluster_expression$GENE
  # Compute and the store the frequency table
  RGC_cluster_hist[[clust]] <- ExpressionTableMaker(RGCmat, gene_subset, n_vals = n_vals)
  RGC_cluster_hist[['N']][clust] <- length(RGCmat)-1
  RGC_cluster_hist[['N ID']][clust] <- clust
  # Compute summary statistics 
  RGCmat <- select(RGCmat, -GENE)
  RGC_cluster_expression[clust] <- rowMeans(RGCmat)
  RGC_cluster_SD[clust] <- rowSds(as.matrix(RGCmat))
  RGC_cluster_median[clust] <- rowMedians(as.matrix(RGCmat))
  RGC_cluster_expressing[clust] <- rowMeans(RGCmat != 0)
  # Print to tell user all is good
  print(sprintf("%s completed, %s cells present to make frequency table from", clust, length((RGCmat))))
}

# Save the list of frequency dataframes for easy loading elsewhere
save(RGC_cluster_hist, file="RGC_cluster_ExpressionFrequencies_allGenes.RData")
# Save summary stats
save(RGC_cluster_expression, RGC_cluster_expressing, RGC_cluster_SD, RGC_cluster_median, file="RGC_cluster_ExpressionMats_allGenes.RData")

# Save the gene frequency tables in an excel sheet where each sheet is a cluster
library(openxlsx)

# create workbook
wb <- createWorkbook()

#Iterate through each list element to create a sheet based on that list (creating an anonymous function inside Map())
Map(function(data, nameofsheet){     
  addWorksheet(wb, nameofsheet)
  writeData(wb, nameofsheet, data)
}, RGC_cluster_hist, names(RGC_cluster_hist))

## Save workbook to excel file 
saveWorkbook(wb, file = "RGC_cluster_hist_allGenes.xlsx", overwrite = TRUE)
