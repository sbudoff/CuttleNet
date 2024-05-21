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

# The below computes a df of distance scores by cellclass for each gene. 
# The idea is to sort this list and cut off genes with very low distance 
#################################################################################
dist_table_maker <- function(Retina_freq_table_expressingGenes, verbose = F)
  {
  # THis function of code computes, visualizes and stores the 
  dist_table <-  data.frame(matrix(ncol = 3, nrow = 0))
  colnames(dist_table) <- c('CellType', 'Gene', 'Distance')
  for (celltype in unique(Retina_freq_table_expressingGenes$CellType)) 
    {
    Cell_freqs <- Retina_freq_table_expressingGenes %>%
    filter(CellType == celltype)
    for (gene in unique(Cell_freqs$Gene)) 
      {
      test <- Cell_freqs %>%
               filter(Gene == gene) %>%
        select(CellClass,NormalizedFreq,ExpressionLevel) %>%
        pivot_wider(names_from = ExpressionLevel, id_cols = CellClass, values_from = NormalizedFreq) %>%
        column_to_rownames("CellClass")
      
      
      d_mat = data.frame(matrix(0, length(rownames(test)), length(rownames(test))), 
                         row.names = rownames(test))
      colnames(d_mat) = rownames(test)
      current_row <- c(celltype,gene,sum(d_mat))
      print(current_row)
      dist_table <- rbind(dist_table, current_row)
      for (cell_i in rownames(test)) 
        {
        for (cell_j in rownames(test)) 
        {
          d_mat[cell_i,cell_j] = sum((test[cell_i,] - test[cell_j,])**2)
        }
        if (verbose) {
          barcode <- d_mat %>% 
            as.data.frame() %>% 
            rownames_to_column(var = "class") %>% 
            gather(Class, val, -class) %>% 
            ggplot(aes(Class, class)) + 
            geom_tile(aes(fill = val)) + 
            coord_fixed() + 
            guides(fill = "none") +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
            ggtitle(gene)
          print(barcode)
          
          print(sum(d_mat))
        }
      }
    }
    print(sprintf("%s complete", celltype))
  }
}
dist_table <- dist_table_maker(Retina_freq_table_expressingGenes)
save(dist_table, file='dist_table.RData')
#################################################################################


cell='Glia'
Retina_freq_table_expressingGenes %>%
  filter(CellType == cell) %>%
  mutate(ExpressionLevel = case_when(ExpressionLevel == -0.01 ~ 0,
                                     T ~ ExpressionLevel),
         Val = ExpressionLevel*NormalizedFreq) %>%
  group_by(CellType, CellClass, Gene) %>%
  mutate(Mean = mean(Val)) %>%
  ggplot(aes(x = CellClass, y = Gene, size = Mean)) +
    geom_point()




test <- Retina_freq_table_expressingGenes %>%
  mutate(ExpressionLevel = case_when(ExpressionLevel == -0.01 ~ 0,
                                     T ~ ExpressionLevel),
         Val = ExpressionLevel*NormalizedFreq) %>%
  group_by(CellType, CellClass, Gene) %>%
  summarize(N=N, Mean = mean(Val)) %>%
  ungroup() %>%
  filter (! duplicated(Mean)) %>%

test %>%
  group_by(CellType) %>%
  ggplot(aes(x=CellClass, y= Gene, size = Mean)) +
  geom_point() +
  facet_grid(~CellType)


for (i in 1:100) {
  test <- Retina_freq_table_expressingGenes %>%
    filter(CellType == 'RGC') %>%
    filter(Gene == 'Spp1') %>%
    ggplot(aes(x=ExpressionLevel, y= NormalizedFreq, color = CellClass)) +
    geom_line()
  print(test)
}

test <- Retina_freq_table_expressingGenes %>%
  group_by(CellType, CellClass, ExpressionLevel) %>%
  summarise(max(Frequency))

bin_levels <- unique(Retina_freq_table_expressingGenes$ExpressionLevel)
test2 <- Retina_freq_table_expressingGenes %>%
  group_by(CellType, CellClass, Gene, ExpressionLevel) %>%
  filter(NormalizedFreq > 0.2) %>%
  ungroup() %>%
  group_by(CellType, CellClass, Gene) %>%
  summarise(TopBin = max(ExpressionLevel)) %>%
  ungroup() %>%
  arrange(CellType,Gene, CellClass)

test<- left_join(test2, Retina_freq_table_expressingGenes, by = c("CellType" = "CellType", 
                                                                  "CellClass" = "CellClass",
                                                                  "Gene" = "Gene",
                                                                  "TopBin" = "ExpressionLevel"))

test %>%
  ggplot(aes(x = TopBin, y = NormalizedFreq, color = CellClass)) +
    geom_point() +
    facet_grid(~CellType) + 
    theme(legend.position = "none")


for (bin in bin_levels) {
  plot_test <- test2 %>%
    ggplot(aes(x=ExpressionLevel, y = val, color = CellType)) +
    geom_point() +
     facet_wrap(~CellType)
}