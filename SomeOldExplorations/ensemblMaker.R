# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("biomaRt")

library(biomaRt)
library(dplyr)

ensembl=useMart("ENSEMBL_MART_ENSEMBL")
listDatasets(ensembl)[1]

ensembl = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

ensembl = useEnsembl(biomart = "genes",  dataset = "mmusculus_gene_ensembl")


root <- "/home/sam/scRNAseq/Xenium/"
setwd(root)
gene_path <- c('tier_3_genes.csv' ,'extra_57.csv', 'top_225_genes.csv')

gene_path <- c('tier_4_genes.csv')
for (gp in gene_path){
  topNgenes <- read.csv(gp, head=FALSE) %>%
    pull(V1) %>%
    tolower()
  
  
  ensembl_ids <- getBM(attributes = c("external_gene_name","ensembl_gene_id"),
                       filters = 'mgi_symbol',
                       values = topNgenes,
                       mart = ensembl)
  
  # Specify the file path where you want to save the CSV file
  output_file <- paste0(root,"ensembl_",gp)
  
  # Save the data frame as a CSV file
  write.csv(ensembl_ids, file = output_file, row.names = FALSE)
  
  print("Failed to match:")
  print(topNgenes[!(tolower(topNgenes) %in% tolower(ensembl_ids$external_gene_name))])
  
}


 