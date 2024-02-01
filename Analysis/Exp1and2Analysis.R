library(tidyverse)
library(patchwork)

in_path <- '/home/sam/scRNAseq/Xenium/ClassVsSubclass/'

network_results = read_csv(paste0(in_path, 'Exp1and2_Network_Results_withClass.csv')) %>%
  mutate(Experiment = as.factor(Experiment),
         Class = as.factor(Class),
         Cluster = as.factor(Cluster))

str(network_results)

ClassOnly <- function(network_results) {
  network_results %>%
    select(-c(Cluster, cells, mean_TP, sd_TP, mean_TN, sd_TN, mean_FP, sd_FP, mean_FN,
              sd_FN, mean_prec, sd_prec, mean_tpr, sd_tpr, mean_F1, sd_F1, cell_fraction)) %>%
    unique()
}

# Histogram of Mean True Positive Rate (TPR)
ggplot(network_results, aes(x = mean_tpr)) + 
  geom_histogram(bins = 30, fill = "blue", color = "black") +
  ggtitle("Distribution of Mean True Positive Rate (TPR)") +
  facet_grid(~Experiment)

# Histogram of Mean F1 Score
ggplot(network_results, aes(x = mean_F1)) + 
  geom_histogram(bins = 30, fill = "green", color = "black") +
  ggtitle("Distribution of Mean F1 Score")+
  facet_grid(~Experiment)

# Histogram of Mean Precision
ggplot(network_results, aes(x = mean_prec)) + 
  geom_histogram(bins = 30, fill = "red", color = "black") +
  ggtitle("Distribution of Mean Prec Score")+
  facet_grid(~Experiment)


# Histogram of Mean True Positive Rate (TPR)
ggplot(ClassOnly(network_results), aes(x = mean_tpr_class)) + 
  geom_histogram(bins = 30, fill = "blue", color = "black") +
  ggtitle("Distribution of Mean True Positive Rate (TPR) By Class")+
  facet_grid(~Experiment)

# Histogram of Mean F1 Score
ggplot(ClassOnly(network_results), aes(x = mean_F1_class)) + 
  geom_histogram(bins = 30, fill = "green", color = "black") +
  ggtitle("Distribution of Mean F1 Score By Class")+
  facet_grid(~Experiment)

# Histogram of Mean Precision
ggplot(ClassOnly(network_results), aes(x = mean_prec_class)) + 
  geom_histogram(bins = 30, fill = "red", color = "black") +
  ggtitle("Distribution of Mean Prec Score By Class")+
  facet_grid(~Experiment)

# Comparing mean F1 scores across experiments
ggplot(network_results, aes(x = Experiment, y = mean_F1, fill = Experiment)) + 
  geom_boxplot() +
  ggtitle("Comparison of Mean F1 Scores Across Experiments")

# Plotting ClassWeight/DeltaRate against Mean F1 Score for class weight data
p1 <- network_results %>% 
  filter(Experiment == "ClassWeight") %>%
  ggplot(aes(x = ClassWeight, y = mean_F1, color = Class)) + 
  geom_point(position = position_dodge(width = 0.1), alpha = 0.8) +  # Add an offset
  scale_x_log10() +
  ggtitle("ClassWeight vs Mean F1 Score (Class Weight Data)")

p2 <- network_results %>% 
  ClassOnly() %>%
  filter(Experiment == "ClassWeight") %>%
  ggplot(aes(x = ClassWeight, y = mean_F1_class, color = Class)) + 
  geom_point(position = position_dodge(width = 0.1), alpha = 0.8) +  # Add an offset
  scale_x_log10() +
  ggtitle("ClassWeight vs Mean F1 Score (Class Weight Data) Class Only")

p1+p2

# Plotting ClassWeight/DeltaRate against Mean Prec Score for class weight data
p1 <- network_results %>% 
  filter(Experiment == "ClassWeight") %>%
  ggplot(aes(x = ClassWeight, y = mean_prec, color = Class)) + 
  geom_point(position = position_dodge(width = 0.1), alpha = 0.8) +  # Add an offset
  scale_x_log10() +
  ggtitle("ClassWeight vs Mean Prec Score (Class Weight Data)")

p2 <- network_results %>% 
  ClassOnly() %>%
  filter(Experiment == "ClassWeight") %>%
  ggplot(aes(x = ClassWeight, y = mean_prec_class, color = Class)) + 
  geom_point(position = position_dodge(width = 0.1), alpha = 0.8) +  # Add an offset
  scale_x_log10() +
  ggtitle("ClassWeight vs Mean Prec Score (Class Weight Data) Class Only")

p1+p2

# Plotting ClassWeight/DeltaRate against Mean TPR Score for class weight data
p1 <- network_results %>% 
  filter(Experiment == "ClassWeight") %>%
  ggplot(aes(x = ClassWeight, y = mean_tpr, color = Class)) + 
  geom_point(position = position_dodge(width = 0.1), alpha = 0.8) +  # Add an offset
  scale_x_log10() +
  ggtitle("ClassWeight vs Mean TPR Score (Class Weight Data)")

p2 <- network_results %>% 
  ClassOnly() %>%
  filter(Experiment == "ClassWeight") %>%
  ggplot(aes(x = ClassWeight, y = mean_tpr_class, color = Class)) + 
  geom_point(position = position_dodge(width = 0.1), alpha = 0.8) +  # Add an offset
  scale_x_log10() +
  ggtitle("ClassWeight vs Mean TPR Score (Class Weight Data) Class Only")

p1+p2
