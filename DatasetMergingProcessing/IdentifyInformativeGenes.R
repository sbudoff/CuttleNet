library(tidyverse)
library(matrixStats) 
setwd('/home/sam/scRNAseq/FreqTables/')
load('Retina_Freq_Table_ExpressingGenes.RData')

# This chunk of code creates a list of gene by cellType matrixes with mean expression value 
gene_mats = list()
for (cell in unique(Retina_freq_table_expressingGenes$CellType)) {
  gene_mats[[cell]] <- Retina_freq_table_expressingGenes %>%
    filter(CellType == cell) %>%
    mutate(ExpressionLevel = case_when(ExpressionLevel == -0.01 ~ 0,
                                       T ~ ExpressionLevel),
           Val = ExpressionLevel*NormalizedFreq) %>%
    group_by(CellType, CellClass, Gene) %>%
    summarize(Mean = mean(Val)) %>%
    ungroup() %>%
    filter (! duplicated(Mean)) %>%
    select(CellClass,Gene,Mean) %>%
    pivot_wider(names_from = CellClass,  id_cols = Gene,values_from = Mean) %>%
    column_to_rownames("Gene")
  print(sprintf('Made the %s matrix', cell))
}

# The below computes The expression Matrices on a gene by cellclass  
#################################################################################
# This chunk of code creates a list of gene by cellType matrixes with mean expressing value 
gene_expression_mats = list()
gene_expression_binary = list()
ubiquitous_genes = list()
threshold = 0.5
for (cell in unique(Retina_freq_table_expressingGenes$CellType)) {
  expr_mat <- Retina_freq_table_expressingGenes %>%
    filter(CellType == cell, ExpressionLevel != -0.01 ) %>%
    group_by(CellType, CellClass, Gene) %>%
    summarize(Expr = sum(Frequency),
              N = N) %>%
    filter (! duplicated(Expr)) %>%
    summarize(Expression = Expr/N)%>%
    ungroup() %>%
    select(CellClass,Gene,Expression) %>%
    pivot_wider(names_from = CellClass,  id_cols = Gene,values_from = Expression) %>%
    column_to_rownames("Gene")
  
  bin_mat <- expr_mat >= threshold
  gene_expression_mats[[cell]] <- expr_mat
  gene_expression_binary[[cell]] <- bin_mat
  ubiquitous_genes[[cell]] <- names(which(rowSums(bin_mat) == ncol(expr_mat)))
  
  print(sprintf('Made the %s expression matrix', cell))
}


houseKeepingGenes <- table(unlist(ubiquitous_genes))
houseKeepingGenes <- names(houseKeepingGenes[houseKeepingGenes>4])
neuralGenes <- ubiquitous_genes
neuralGenes[['Glia']] <- NULL
neuralGenes <- table(unlist(neuralGenes))
neuralGenes <- names(neuralGenes[neuralGenes>3])

expressed_genes <- list()
for (cell in unique(Retina_freq_table_expressingGenes$CellType)) {
  expressed_genes[[cell]] <- Retina_freq_table_expressingGenes %>%
    filter(CellType == cell) %>%
    pull(Gene)
  
  print(sprintf('Pulled the %s expressed genes', cell))
}

expressed_genes <- unique(unlist(expressed_genes))

useful_genes <- expressed_genes[!expressed_genes %in% neuralGenes]

save(useful_genes, houseKeepingGenes, neuralGenes, ubiquitous_genes, gene_expression_binary, gene_expression_mats, gene_mats, file = 'ExpressingExpressionMats.RData')
