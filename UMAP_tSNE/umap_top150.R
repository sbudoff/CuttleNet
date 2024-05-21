library(tidyverse)

load('/home/sam/scRNAseq/Xenium/Retina_expMatrix_clean.RData')
topNgenes <- read.csv('/home/sam/scRNAseq/Xenium/top_225_genes.csv', head=FALSE) %>%
  pull(V1)

Retina_150 <- Retina_expMatrix_candidateGenes %>%
  select(Cluster, all_of(topNgenes))

round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  (df)
}

Retina_150 <- round_df(Retina_150, digits=2)

library(umap)

retina_umap <- Retina_150 %>%
  select(-Cluster) %>%
  umap()

retina_umap_coordinates <- retina_umap[['layout']] %>%
  data.frame() %>%
  mutate(Cluster = Retina_150$Cluster,
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

retina_umap_colors <- retina_umap_coordinates %>%
  group_by(Class) %>%
  select(Cluster) %>%
  distinct() %>%
  mutate(n = 1:n(),
         Color = distinct_colors[n]) %>%
  select(-n) %>%
  ungroup()

retina_umap_coordinates <- left_join(retina_umap_coordinates,
                                     retina_umap_colors, by=c("Cluster", "Class"))


RGC_umap <- retina_umap_coordinates %>%
  filter(Class == "RGC") 
BC_umap <- retina_umap_coordinates %>%
  filter(Class == "BC") 
AC_umap <- retina_umap_coordinates %>%
  filter(Class == "AC")
Ph_umap <- retina_umap_coordinates %>%
  filter(Class == "Ph") 

alpha_bkgrnd = 0.05
alpha_frgrnd = 0.4


RGC_plot <- retina_umap_coordinates %>% 
  filter(Class != "Glia",
         Class != "RGC") %>%
  ggplot(aes(x=X,y=Y)) +
  geom_point(alpha = alpha_bkgrnd) +
  geom_point(data = RGC_umap, aes(x=X,y=Y,color = Cluster), alpha = alpha_frgrnd) +
  xlim(-15,15) +
  theme_minimal() 

AC_plot <- retina_umap_coordinates %>% 
  filter(Class != "Glia",
         Class != "AC") %>%
  ggplot(aes(x=X,y=Y)) +
  geom_point(alpha = alpha_bkgrnd) +
  geom_point(data = AC_umap, aes(x=X,y=Y,color = Cluster), alpha = alpha_frgrnd) +
  xlim(-15,15) +
  theme_minimal() 

BC_plot <- retina_umap_coordinates %>% 
  filter(Class != "Glia",
         Class != "BC",) %>%
  ggplot(aes(x=X,y=Y)) +
  geom_point(alpha = alpha_bkgrnd) +
  geom_point(data = BC_umap, aes(x=X,y=Y,color = Cluster), alpha = alpha_frgrnd) +
  xlim(-15,15) +
  theme_minimal() 

Ph_plot <- retina_umap_coordinates %>% 
  filter(Class != "Glia",
         Class != "Ph") %>%
  ggplot(aes(x=X,y=Y)) +
  geom_point(alpha = alpha_bkgrnd) +
  geom_point(data = Ph_umap, aes(x=X,y=Y,color = Cluster), alpha = alpha_frgrnd) +
  xlim(-15,15) +
  theme_minimal() 

library(patchwork)

(Ph_plot + BC_plot) / (AC_plot + RGC_plot)

###################################################################################
BC_umap <- Retina_150 %>%
  mutate(Class = case_when(substr(Cluster, 1, 2) == "AC" ~ "AC",
                           substr(Cluster, 1, 2) == "BC" ~ "BC",
                           substr(Cluster, 1, 2) == "RB" ~ "BC",
                           substr(Cluster, 1, 2) == "Ro" ~ "Ph",
                           substr(Cluster, 1, 2) == "Co" ~ "Ph",
                           substr(Cluster, 1, 2) == "MG" ~ "Glia",
                           T~ "RGC")) %>%
  filter(Class == "BC") 

BC_umap_labels <- BC_umap$Cluster

BC_umap <- BC_umap %>%
  select(-Cluster, -Class) %>%
  umap()

BC_umap_coordinates <- BC_umap[['layout']] %>%
  data.frame() %>%
  mutate(Cluster = BC_umap_labels,
         X = X1,
         Y = X2) %>%
  select(-X1,-X2) 

BC_umap_coordinates %>%
  ggplot(aes(x=X,y=Y, color = Cluster)) +
  geom_point(alpha = alpha_frgrnd) +
  theme_minimal() 

###################################################################################
AC_umap <- Retina_150 %>%
  mutate(Class = case_when(substr(Cluster, 1, 2) == "AC" ~ "AC",
                           substr(Cluster, 1, 2) == "BC" ~ "BC",
                           substr(Cluster, 1, 2) == "RB" ~ "BC",
                           substr(Cluster, 1, 2) == "Ro" ~ "Ph",
                           substr(Cluster, 1, 2) == "Co" ~ "Ph",
                           substr(Cluster, 1, 2) == "MG" ~ "Glia",
                           T~ "RGC")) %>%
  filter(Class == "AC") 

AC_umap_labels <- AC_umap$Cluster

AC_umap <- AC_umap %>%
  select(-Cluster, -Class) %>%
  umap()

AC_umap_coordinates <- AC_umap[['layout']] %>%
  data.frame() %>%
  mutate(Cluster = AC_umap_labels,
         X = X1,
         Y = X2) %>%
  select(-X1,-X2) 

AC_umap_coordinates %>%
  ggplot(aes(x=X,y=Y, color = Cluster)) +
  geom_point(alpha = alpha_frgrnd) +
  theme_minimal()

###################################################################################
RGC_umap <- Retina_150 %>%
  mutate(Class = case_when(substr(Cluster, 1, 2) == "AC" ~ "AC",
                           substr(Cluster, 1, 2) == "BC" ~ "BC",
                           substr(Cluster, 1, 2) == "RB" ~ "BC",
                           substr(Cluster, 1, 2) == "Ro" ~ "Ph",
                           substr(Cluster, 1, 2) == "Co" ~ "Ph",
                           substr(Cluster, 1, 2) == "MG" ~ "Glia",
                           T~ "RGC")) %>%
  filter(Class == "RGC") 

RGC_umap_labels <- RGC_umap$Cluster

RGC_umap <- RGC_umap %>%
  select(-Cluster, -Class) %>%
  umap()

RGC_umap_coordinates <- RGC_umap[['layout']] %>%
  data.frame() %>%
  mutate(Cluster = RGC_umap_labels,
         X = X1,
         Y = X2) %>%
  select(-X1,-X2) 

RGC_umap_coordinates %>%
  ggplot(aes(x=X,y=Y, color = Cluster)) +
  geom_point(alpha = alpha_frgrnd) +
  theme_minimal() 
###################################################################################

save(RGC_umap_coordinates, AC_umap_coordinates, BC_umap_coordinates, 
     retina_umap_coordinates, file = "/home/sam/scRNAseq/Xenium/umap_coordinate_dfs.RData")

