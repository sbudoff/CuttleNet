library(tidyverse)
library(matrixStats)
library(Matrix)

nDELTArow2col2 <- function(expMat){
  # This loss function computes a value between -Inf and Inf where the lower is the better separation
  expMat <- na.omit(expMat)
  n = nrow(expMat)
  classSum2 = sum(rowSums(t(expMat))**2)
  geneSum2 = sum(rowSums(expMat)**2)
  l = n*(geneSum2-classSum2)
}
barcoder <- function(class_mat, verbose = F, binary = F){
  # This Function creates a 'barcode' for each vector of genes defining a class. 
  # This barcode is simply a square matrix of booleans that checks each class against every other class to see if they share a vector of genes.
  # If verbose is true, a visualization of the barcode will be displayed
  class_mat <- na.omit(class_mat)
  n_classes = nrow(t(class_mat))
  n_genes = nrow(class_mat)
  barcode_compare <- data.frame(matrix(ncol = n_classes, nrow = n_classes))
  rownames(barcode_compare) <- colnames(class_mat)
  colnames(barcode_compare) <- colnames(class_mat)
  for (i in 1:n_classes) {
    barcode_i = round(class_mat[,i],1)
    for (j in 1:n_classes) {
      barcode_j = round(class_mat[,j],1)
      if (binary) {
        barcode_compare[i,j] = sum(barcode_i == barcode_j) == nrow(class_mat)
      } else {
        barcode_compare[i,j] = sqrt(sum(((barcode_i-barcode_j)**2)))
      }
    }
  }
  if (verbose) {
    barcode <- barcode_compare %>% 
      as.data.frame() %>% 
      rownames_to_column(var = "class") %>% 
      gather(Class, val, -class) %>% 
      ggplot(aes(Class, class)) + 
      geom_tile(aes(fill = val)) + 
      coord_fixed() + 
      guides(fill = "none") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    print(barcode)
  }
  sum(rowSums(barcode_compare))
}

RedundantRows <- function(df, precision = 1) {
  df <- round(df,1)
  out <- c(1:nrow(df))
  for (i in 1:nrow(df)) {
    out[i] = length(unique(as.numeric(df[i,])))
  }
  which(out < 1)
}

PlotDistance <- function(cluster_expression, cluster_expressing, GOIs, verbose = F) {
  thresh <- t(scale(t(as.matrix(cluster_expression))) * scale(t(as.matrix(cluster_expressing))))
  
  GOIs_ind <- which(rownames(thresh) %in% GOIs)
  
  thresh <- thresh[GOIs_ind,]
  
  barcoder(thresh, verbose = verbose) 
}

PlotBubbles <- function(cluster_expression, cluster_expressing, GOIs) {
  cluster_by_gene <- cluster_expression %>% 
    rownames_to_column('GENE') %>%
    filter(GENE %in% GOIs) %>%
    as.data.frame() %>%
    gather(Class, mean, -GENE)
  cluster_by_gene <- cluster_expressing %>% 
    rownames_to_column('GENE') %>%
    filter(GENE %in% GOIs) %>%
    as.data.frame() %>%
    gather(Class, expressing, -GENE) %>%
    merge(cluster_by_gene, by = c("GENE", "Class"))
  
  # visualize the resulting unique clusters
  visual_cluster <- cluster_by_gene %>% 
    ggplot(aes(Class, GENE)) + 
    geom_point(aes(color = mean, size = expressing)) + 
    coord_fixed() + 
    guides(fill = "none") +
    scale_color_gradient(low = "blue", high = "red") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  print(visual_cluster)
}

GOIextractor <- function(data_path, cellClass = 'RGC', split = 0.8) {
  load(data_path)
  
  cluster_expressing <- gene_mats[[cellClass]]
  cluster_expression <- gene_expression_mats[[cellClass]]
  cluster_expressing[is.na(cluster_expressing)] <- 0
  cluster_expression[is.na(cluster_expression)] <- 0
  
  GOIs <- cluster_expression %>% 
    rownames_to_column('GENE') %>%
    as.data.frame() %>%
    gather(Class, expressing, -GENE) %>%
    group_by(GENE) %>%
    summarize(high = max(expressing),
              low = min(expressing)) %>%
    mutate(split = high-low) %>%
    filter(split > 0.8) %>%
    pull(GENE)
}

################################################################################

setwd('/home/sam/scRNAseq/FreqTables/')
data_path = 'ExpressingExpressionMats.RData'

key_genes <- c('Rbpms', 'Opn1sw','Grm6', 'Prkcq', 'Coch', 'Spp1', 'Meis2', 'Mafb', 'Penk', 'Opn4', 'Etv1', 'Zic1')
GOIs <- key_genes
GOIs <- c(GOIs, GOIextractor(data_path, cellClass = 'RGC', split = 0.8))
GOIs <- c(GOIs, GOIextractor(data_path, cellClass = 'AC', split = 0.8))
GOIs <- c(GOIs, GOIextractor(data_path, cellClass = 'BP', split = 0.8))

GOIs <- unique(GOIs)

save(key_genes, GOIs, file = 'GeneShortList.RData')

################################################################################

load(data_path)
cellClass = 'RGC'

cluster_expressing <- gene_mats[[cellClass]]
cluster_expression <- gene_expression_mats[[cellClass]]
cluster_expressing[is.na(cluster_expressing)] <- 0
cluster_expression[is.na(cluster_expression)] <- 0




key_genes <- c('Rbpms', 'Opn1sw','Grm6', 'Prkcq', 'Coch', 'Spp1', 'Meis2', 'Mafb', 'Penk', 'Opn4', 'Etv1', 'Zic1')
n_probes = 24
probes = length(key_genes):n_probes-2  
for (probe in probes) {
  dist = 0
  for (gen in 1:length(GOIs)) {
    gene_list = c(key_genes, GOIs[gen])
    dist_i = PlotDistance(cluster_expression, cluster_expressing, gene_list, verbose = F)
    if (dist_i > dist) {
      best = gen
      dist = dist_i
    }
  }
 print(sprintf("%s added with new distance score of %s", GOIs[best], dist))
 key_genes = c(key_genes, GOIs[best])
 GOIs = GOIs[-best]
 PlotDistance(cluster_expression, cluster_expressing, key_genes, verbose = T)
}
key_genes
PlotBubbles(cluster_expressing, cluster_expression, GOIs) 

test <- cluster_expressing[key_genes,] %>%
  na.omit() %>%
  rownames_to_column('GENE') %>%
  gather(Class, val, -GENE) %>%
  group_by(Class) %>%
  mutate(mean = mean(val),
         val = val/mean) %>%
  select(-mean) %>%
  ungroup() %>%
  spread(Class, val) %>%
  column_to_rownames('GENE') %>%
  round()

test2 <- cluster_expression[key_genes,] %>%
  na.omit() %>%
  rownames_to_column('GENE') %>%
  gather(Class, val, -GENE) %>%
  group_by(Class) %>%
  mutate(mean = mean(val),
         val = val/mean) %>%
  select(-mean) %>%
  ungroup() %>%
  spread(Class, val) %>%
  column_to_rownames('GENE')

PlotBubbles(test, test2, key_genes) 

test3<- test[ order(max.col(test, "first")), ]

library(umap)


umap_list <- test %>%
  t() %>%
  scale() %>%
  umap() 

umap_df <- umap_list[["layout"]] %>%
  as.data.frame()%>%
  rename(UMAP1="V1",
         UMAP2="V2")

umap_df %>%
  ggplot(aes(x = UMAP1,  y = UMAP2))+
  geom_point()
