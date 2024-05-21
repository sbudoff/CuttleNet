library(tidyverse)

load('/home/sam/scRNAseq/Xenium/Retina_expMatrix_clean_final.RData')
Retina_interesting_genes <- read.csv('Xenium Gene List_final.csv') %>%
  mutate_all(~str_replace_all(., " ", "")) %>%  # Remove spaces in all columns
  mutate(gene = str_replace_all(gene, "-", "."))
topNgenes <- Retina_interesting_genes$gene

Retina_300 <- Retina_expMatrix_candidateGenes %>%
  select(Cluster, all_of(topNgenes))

round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  (df)
}

Retina_300 <- round_df(Retina_300, digits=2)

library(Rtsne)

retina_tsne <- Retina_300 %>%
  select(-Cluster) %>%
  Rtsne()

retina_tsne_coordinates <- retina_tsne[['Y']] %>%
  data.frame() %>%
  mutate(Cluster = Retina_300$Cluster,
         X = X1,
         Y = X2,
         Class = case_when(substr(Cluster, 1, 2) == "AC" ~ "AC",
                           substr(Cluster, 1, 2) == "BC" ~ "BC",
                           substr(Cluster, 1, 2) == "RB" ~ "BC",
                           substr(Cluster, 1, 2) == "Ro" ~ "Ph",
                           substr(Cluster, 1, 2) == "Co" ~ "Ph",
                           substr(Cluster, 1, 2) == "MG" ~ "Glia",
                           T~ "RGC")) %>%
  select(-X1,-X2) 

distinct_colors <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
  "#393b79", "#31a354", "#756bb1", "#252525", "#fed9a6", 
  "#fdd0a2", "#525252", "#084081", "#b30000", "#7f0000", 
  "#bdbdbd", "#969696", '#000075', "#f768a1", "#599861", 
  "#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", 
  "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", 
  "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", 
  "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", 
  "#C84248", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", 
  '#3cb44b', '#ffe119', '#911eb4', '#46f0f0', '#f032e6', 
  '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', 
  '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', 
  "#808080"
)

retina_tsne_colors <- retina_tsne_coordinates %>%
  group_by(Class) %>%
  select(Cluster) %>%
  distinct() %>%
  mutate(n = 1:n(),
         Color = distinct_colors[n]) %>%
  select(-n) %>%
  ungroup()

retina_tsne_coordinates <- left_join(retina_tsne_coordinates,
                                     retina_tsne_colors, by=c("Cluster", "Class"))


RGC_tsne <- retina_tsne_coordinates %>%
  filter(Class == "RGC") 
BC_tsne <- retina_tsne_coordinates %>%
  filter(Class == "BC") 
AC_tsne <- retina_tsne_coordinates %>%
  filter(Class == "AC")
Ph_tsne <- retina_tsne_coordinates %>%
  filter(Class == "Ph") 

alpha_bkgrnd = 0.05
alpha_frgrnd = 0.4


RGC_plot <- retina_tsne_coordinates %>% 
  filter(Class != "Glia",
         Class != "RGC") %>%
  ggplot(aes(x=X,y=Y)) +
  geom_point(alpha = alpha_bkgrnd) +
  geom_point(data = RGC_tsne, aes(x=X,y=Y,color = Cluster), alpha = alpha_frgrnd) +
  theme_minimal() 

AC_plot <- retina_tsne_coordinates %>% 
  filter(Class != "Glia",
         Class != "AC") %>%
  ggplot(aes(x=X,y=Y)) +
  geom_point(alpha = alpha_bkgrnd) +
  geom_point(data = AC_tsne, aes(x=X,y=Y,color = Cluster), alpha = alpha_frgrnd) +
  theme_minimal() 

BC_plot <- retina_tsne_coordinates %>% 
  filter(Class != "Glia",
         Class != "BC") %>%
  ggplot(aes(x=X,y=Y)) +
  geom_point(alpha = alpha_bkgrnd) +
  geom_point(data = BC_tsne, aes(x=X,y=Y,color = Cluster), alpha = alpha_frgrnd) +
  theme_minimal() 

Ph_plot <- retina_tsne_coordinates %>% 
  filter(Class != "Glia",
         Class != "Ph") %>%
  ggplot(aes(x=X,y=Y)) +
  geom_point(alpha = alpha_bkgrnd) +
  geom_point(data = Ph_tsne, aes(x=X,y=Y,color = Cluster), alpha = alpha_frgrnd) +
  theme_minimal() 

library(patchwork)

(Ph_plot + BC_plot) / (AC_plot + RGC_plot)

###################################################################################
BC_tsne <- Retina_300 %>%
  mutate(Class = case_when(substr(Cluster, 1, 2) == "AC" ~ "AC",
                           substr(Cluster, 1, 2) == "BC" ~ "BC",
                           substr(Cluster, 1, 2) == "RB" ~ "BC",
                           substr(Cluster, 1, 2) == "Ro" ~ "Ph",
                           substr(Cluster, 1, 2) == "Co" ~ "Ph",
                           substr(Cluster, 1, 2) == "MG" ~ "Glia",
                           T~ "RGC")) %>%
  filter(Class == "BC") 

BC_tsne_labels <- BC_tsne$Cluster

BC_tsne <- BC_tsne %>%
  select(-Cluster, -Class) %>%
  Rtsne()

BC_tsne_coordinates <- BC_tsne[['Y']] %>%
  data.frame() %>%
  mutate(Cluster = BC_tsne_labels,
         X = X1,
         Y = X2) %>%
  select(-X1,-X2) 

BC_tsne_coordinates %>%
  ggplot(aes(x=X,y=Y, color = Cluster)) +
  geom_point(alpha = alpha_frgrnd) +
  theme_minimal() 

###################################################################################
AC_tsne <- Retina_300 %>%
  mutate(Class = case_when(substr(Cluster, 1, 2) == "AC" ~ "AC",
                           substr(Cluster, 1, 2) == "BC" ~ "BC",
                           substr(Cluster, 1, 2) == "RB" ~ "BC",
                           substr(Cluster, 1, 2) == "Ro" ~ "Ph",
                           substr(Cluster, 1, 2) == "Co" ~ "Ph",
                           substr(Cluster, 1, 2) == "MG" ~ "Glia",
                           T~ "RGC")) %>%
  filter(Class == "AC") 

AC_tsne_labels <- AC_tsne$Cluster

AC_tsne <- AC_tsne %>%
  select(-Cluster, -Class) %>%
  Rtsne()

AC_tsne_coordinates <- AC_tsne[['Y']] %>%
  data.frame() %>%
  mutate(Cluster = AC_tsne_labels,
         X = X1,
         Y = X2) %>%
  select(-X1,-X2) 

AC_tsne_coordinates %>%
  ggplot(aes(x=X,y=Y, color = Cluster)) +
  geom_point(alpha = alpha_frgrnd) +
  theme_minimal()

###################################################################################
RGC_tsne <- Retina_300 %>%
  mutate(Class = case_when(substr(Cluster, 1, 2) == "AC" ~ "AC",
                           substr(Cluster, 1, 2) == "BC" ~ "BC",
                           substr(Cluster, 1, 2) == "RB" ~ "BC",
                           substr(Cluster, 1, 2) == "Ro" ~ "Ph",
                           substr(Cluster, 1, 2) == "Co" ~ "Ph",
                           substr(Cluster, 1, 2) == "MG" ~ "Glia",
                           T~ "RGC")) %>%
  filter(Class == "RGC") 

RGC_tsne_labels <- RGC_tsne$Cluster

RGC_tsne <- RGC_tsne %>%
  select(-Cluster, -Class) %>%
  Rtsne()

RGC_tsne_coordinates <- RGC_tsne[['Y']] %>%
  data.frame() %>%
  mutate(Cluster = RGC_tsne_labels,
         X = X1,
         Y = X2) %>%
  select(-X1,-X2) 

RGC_tsne_coordinates %>%
  ggplot(aes(x=X,y=Y, color = Cluster)) +
  geom_point(alpha = alpha_frgrnd) +
  theme_minimal() 

###################################################################################
save(RGC_tsne_coordinates, AC_tsne_coordinates, BC_tsne_coordinates, 
     retina_tsne_coordinates, file = "/home/sam/scRNAseq/Xenium/tsne_coordinate_dfs_final.RData")

###################################################################################
###################################################################################
###################################################################################

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



###################################################################################
###################################################################################
###################################################################################
dist_rgc <- RGC_tsne_coordinates %>%
  select(X,Y) %>%
  as.matrix() %>%
  normalize_input() %>%
  dist()
rgc_hierarchy <- hclust(dist_rgc)

rgc_clust_list <- list()
ws <- c()
for (k in 2:70) {
  rgc_clust_list[[k]] <- cutree(rgc_hierarchy, k)
  ws <- c(ws, withinClustVar(select(RGC_tsne_coordinates,-Cluster), rgc_clust_list[[k]]))
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

dist_bc <- BC_tsne_coordinates %>%
  select(X,Y) %>%
  as.matrix() %>%
  normalize_input() %>%
  dist()
bc_hierarchy <- hclust(dist_bc)

bc_clust_list <- list()
ws <- c()
for (k in 2:70) {
  bc_clust_list[[k]] <- cutree(bc_hierarchy, k)
  ws <- c(ws, withinClustVar(BC_tsne_coordinates, bc_clust_list[[k]]))
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

dist_ac <- AC_tsne_coordinates %>%
  select(X,Y) %>%
  as.matrix() %>%
  normalize_input() %>%
  dist()
ac_hierarchy <- hclust(dist_ac)

ac_clust_list <- list()
ws <- c()
for (k in 2:70) {
  ac_clust_list[[k]] <- cutree(ac_hierarchy, k)
  ws <- c(ws, withinClustVar(AC_tsne_coordinates, ac_clust_list[[k]]))
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

bc_w_plot + ac_w_plot + rgc_w_plot
