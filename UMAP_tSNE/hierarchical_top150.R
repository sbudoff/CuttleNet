library( tidyverse)
library(Rtsne)

withinClustVar <- function(df, clusts){
  col_names <- names(df)
  w <- df %>%
    mutate(Cluster = clusts,
           x = pmap(list(!!!syms(col_names)), c)) %>%
    select(Cluster, x) %>%
    group_by(Cluster) %>%
    mutate(x_bar = map_dbl(x, ~ mean(.)) ,
           delta2 = map_dbl(x, ~ sum((. - x_bar)^2)) ,
           ssq = sum(delta2)) %>%
    ungroup() %>%
    summarize(w = sum(ssq)) %>%
    pull(w)
  return(w)
} 
################################################################################
elbow <- function(ws){
  # this function identifies the elbow of an elbow plot given a vector
  scaled <- abs(1-(ws-min(ws))/(max(ws)-min(ws)))
  difcurve <-  scaled-seq(0,1,length.out=length(scaled))
  difcurve == max(difcurve)
}


load('/home/sam/scRNAseq/Xenium/Retina_expMatrix_clean.RData')
topNgenes <- read.csv('/home/sam/scRNAseq/Xenium/top_147_genes.csv', head=FALSE) %>%
  pull(V1)

Retina_150 <- Retina_expMatrix_candidateGenes %>%
  select(Cluster, all_of(topNgenes)) %>%
  mutate(Class = case_when(substr(Cluster, 1, 2) == "AC" ~ "AC",
                           substr(Cluster, 1, 2) == "BC" ~ "BC",
                           substr(Cluster, 1, 2) == "RB" ~ "BC",
                           substr(Cluster, 1, 2) == "Ro" ~ "Ph",
                           substr(Cluster, 1, 2) == "Co" ~ "Ph",
                           substr(Cluster, 1, 2) == "MG" ~ "Glia",
                           T~ "RGC"))

###################################################################################
###################################################################################
###################################################################################
rgc_150 <- Retina_150  %>%
  filter(Class == "RGC") %>%
  select(-Class, -Cluster)

dist_rgc <- rgc_150 %>%
  as.matrix() %>%
  normalize_input() %>%
  dist()
rgc_hierarchy <- hclust(dist_rgc, method = "centroid")

rgc_clust_list <- list()
ws <- c()
for (k in 2:70) {
  rgc_clust_list[[k]] <- cutree(rgc_hierarchy, k)
  ws <- c(ws, withinClustVar(rgc_150, rgc_clust_list[[k]]))
}

best.w=elbow(ws)

rgc_w_plot <- tibble(k = 2:70, w.index=ws, best.w=elbow(ws)) %>%
  ggplot(aes(x=k,y=w.index)) +
  geom_line() +
  geom_point(aes(color=best.w)) +
  scale_color_manual(values=c('#000000', 'red')) +
  ggtitle("RGC")+
  theme(legend.position='none')+
  theme_minimal()

rgc_w_plot
###################################################################################
ac_150 <- Retina_150  %>%
  filter(Class == "AC") %>%
  select(-Class, -Cluster)

dist_ac <- ac_150 %>%
  as.matrix() %>%
  normalize_input() %>%
  dist()
ac_hierarchy <- hclust(dist_ac)

ac_clust_list <- list()
ws <- c()
for (k in 2:70) {
  ac_clust_list[[k]] <- cutree(ac_hierarchy, k)
  ws <- c(ws, withinClustVar(ac_150, ac_clust_list[[k]]))
}

best.w=elbow(ws)

ac_w_plot <- tibble(k = 2:70, w.index=ws, best.w=elbow(ws)) %>%
  ggplot(aes(x=k,y=w.index)) +
  geom_line() +
  geom_point(aes(color=best.w)) +
  scale_color_manual(values=c('#000000', 'red')) +
  ggtitle("AC")+
  theme(legend.position='none')+
  theme_minimal()

ac_w_plot
###################################################################################
bc_150 <- Retina_150  %>%
  filter(Class == "BC") %>%
  select(-Class, -Cluster)

dist_bc <- bc_150 %>%
  as.matrix() %>%
  normalize_input() %>%
  dist()
bc_hierarchy <- hclust(dist_bc)

bc_clust_list <- list()
ws <- c()
for (k in 2:70) {
  bc_clust_list[[k]] <- cutree(bc_hierarchy, k)
  ws <- c(ws, withinClustVar(bc_150, bc_clust_list[[k]]))
}

best.w=elbow(ws)

bc_w_plot <- tibble(k = 2:70, w.index=ws, best.w=elbow(ws)) %>%
  ggplot(aes(x=k,y=w.index)) +
  geom_line() +
  geom_point(aes(color=best.w)) +
  scale_color_manual(values=c('#000000', 'red')) +
  ggtitle("BC")+
  theme(legend.position='none')+
  theme_minimal()

bc_w_plot
###################################################################################