library(tidyverse)
library(matrixStats)
library(Matrix)
library(umap)

# # This commented out code was used to create the Retina Expression Matrix and Cluster to ID directory
# load("/home/sam/scRNAseq/Expresion_80percent_Split/BP.RData")
# load("/home/sam/scRNAseq/Expresion_80percent_Split/AC.RData")
# load("/home/sam/scRNAseq/Expresion_80percent_Split/RGC.RData")
# load('/home/sam/scRNAseq/FreqTables/GeneShortList.RData')
# 
# rownames(RGCexpMatrix) <- GOIs
# rownames(BPexpMatrix) <- GOIs
# BPexpMatrix[is.na(BPexpMatrix)] <- 0
# RGCexpMatrix[is.na(RGCexpMatrix)] <- 0
# BPexpMatrix$GENE <- GOIs
# RGCexpMatrix$GENE <- GOIs
# 
# RetinaExpMatrix <- merge(BPexpMatrix,ACexpMatrix,by = "GENE")
# RetinaExpMatrix <- merge(RetinaExpMatrix,RGCexpMatrix,by="GENE")
# rownames(RetinaExpMatrix) <- RetinaExpMatrix$GENE
# RetinaExpMatrix <- select(RetinaExpMatrix, -GENE)
# 
# BPclusters <- select(BPclusters, ID, Cluster)
# ACclusters <- select(ACclusters, ID, Cluster)
# RGCclusters <- select(RGCclusters, ID, Cluster)
# RetinaClusters <- rbind(BPclusters, ACclusters)
# RetinaClusters <- rbind(RetinaClusters, RGCclusters)
# 
# save(RetinaClusters,RetinaExpMatrix,key_genes, file = "/home/sam/scRNAseq/Expresion_80percent_Split/Retina302Genes.RData")
load("/home/sam/scRNAseq/Expresion_80percent_Split/Retina302Genes.RData")

RetinaGathered <- RetinaExpMatrix %>%
  rownames_to_column('GENE') %>%
  gather(ID, val, -GENE)

umap_list <- RetinaExpMatrix %>%
  t() %>%
  scale() %>%
  umap() 

umap_df <- umap_list[["layout"]] %>%
  as.data.frame()%>%
  rename(UMAP1="V1",
         UMAP2="V2") 
  

umap_df %>%
  ggplot(aes(x = UMAP1,  y = UMAP2))+
  geom_point(alpha = 0.1) +
  ylim(-30,-20) +
  xlim(-10,10)


# This conditions the retina on Rbpms and then makes a histogram for the specified gene
RetinaExpMatrix %>%
  t() %>%
  as.data.frame() %>%
  filter(Rbpms > 0) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column('GENE') %>%
  gather(ID, val, -GENE) %>%
  filter(GENE == 'Spp1' ) %>%
  ggplot(aes(val)) +
    geom_histogram()
