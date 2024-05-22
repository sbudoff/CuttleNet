library(tidyverse)
library(ggsignif)
library(patchwork)
library(grid)
library(png)

# Load the PNG image
img <- readPNG("/home/sam/scRNAseq/Xenium/NeurIPS/Carl_NeurIPS.png")
img_grob <- rasterGrob(img, interpolate = TRUE)

# Convert the grob into a ggplot object
fig_a <- ggplot() +
  annotation_custom(img_grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  theme_void()


data <- read_csv('/home/sam/scRNAseq/Xenium/ClassVsSubclass/Merged_Experiments_FINAL.csv') 

data[] <- lapply(data, function(x) ifelse(is.na(x), "none", x))

data <- data %>%
  mutate(Class = as.factor(Class),
         Cluster = as.factor(Cluster),
         NumEpochs = as.factor(NumEpochs),
         L1Lambda = as.factor(L1Lambda),
         EarlyStopping = as.factor(EarlyStopping),
         deltaRate = as.factor(deltaRate),
         ClassWeight = as.factor(ClassWeight),
         Layer1 = as.factor(Layer1),
         Layer2 = as.factor(Layer2),
         Layer3 = as.factor(Layer3),
         Layers = as.factor(Layers),
         Experiment = ifelse(ImputedGenes==F & 
                               Experiment=='CellTICS' & 
                               Layers==1, paste(Experiment, '1 Layer No 0s'), Experiment),
         Experiment = ifelse(ImputedGenes==F & 
                               Experiment=='CellTICS' & 
                               Layers==5, paste(Experiment, '5 Layer No 0s'), Experiment),
         Experiment = ifelse(ImputedGenes==T & 
                               Experiment=='CellTICS' & 
                               Layers==1, paste(Experiment, '1 Layer'), Experiment),
         Experiment = ifelse(ImputedGenes==T & 
                               Experiment=='CellTICS' & 
                               Layers==5, paste(Experiment, '5 Layer'), Experiment),
         Experiment = ifelse(ImputedGenes==F & 
                               Experiment=='Ablation' & 
                               deltaRate==0.75, paste(Experiment, 'deltaRate'), Experiment),
         Experiment = as.factor(Experiment),)


# Subset the data to the list of best conditions
best_conds <- data %>%
  group_by(NumEpochs, L1Lambda, EarlyStopping, deltaRate, ClassWeight, Experiment, 
           Layer1,Layer2, Layer3, Layers, ImputedGenes) %>%
  reframe(F1 = median(mean_F1),
          F1_sd = mean(sd_F1),
          F1_class = median(mean_F1_class),
          F1_class_sd = mean(sd_F1_class),
          
          prec = median(mean_prec),
          prec_sd = mean(sd_prec),
          prec_class = median(mean_prec_class),
          prec_class_sd = mean(sd_prec_class),
          
          tpr = median(mean_tpr),
          tpr_sd = mean(sd_tpr),
          tpr_class = median(mean_tpr_class),
          tpr_class_sd = mean(sd_tpr_class)
          
  ) %>%
  mutate(score = (F1-F1_sd) + (F1_class-F1_class_sd) +
           (prec-prec_sd) + (prec_class-prec_class_sd) +
           (tpr-tpr_sd) + (tpr_class-tpr_class_sd)) %>%
  group_by(Experiment) %>%
  filter(score == max(score))


OptimalExtractor <- function(data, X = 'deltaRate') {
  best_conds <- data %>%
    group_by(NumEpochs, L1Lambda, EarlyStopping, deltaRate, ClassWeight, Experiment, 
             Layer1,Layer2, Layer3, Layers, ImputedGenes) %>%
    reframe(F1 = median(mean_F1),
            F1_sd = mean(sd_F1),
            F1_class = median(mean_F1_class),
            F1_class_sd = mean(sd_F1_class),
            
            prec = median(mean_prec),
            prec_sd = mean(sd_prec),
            prec_class = median(mean_prec_class),
            prec_class_sd = mean(sd_prec_class),
            
            tpr = median(mean_tpr),
            tpr_sd = mean(sd_tpr),
            tpr_class = median(mean_tpr_class),
            tpr_class_sd = mean(sd_tpr_class)  ) %>%
    mutate(score = (F1-F1_sd) + (F1_class-F1_class_sd) +
             (prec-prec_sd) + (prec_class-prec_class_sd) +
             (tpr-tpr_sd) + (tpr_class-tpr_class_sd)) %>%
    group_by(Experiment) %>%
    filter(score == max(score))
  
  out_df <- data %>%
    filter(NumEpochs == pull(filter(best_conds, Experiment==X),   NumEpochs),
           L1Lambda == pull(filter(best_conds, Experiment==X),   L1Lambda),
           EarlyStopping == pull(filter(best_conds, Experiment==X),   EarlyStopping),
           deltaRate == pull(filter(best_conds, Experiment==X),   deltaRate),
           ClassWeight == pull(filter(best_conds, Experiment==X),   ClassWeight),
           Layer1 == pull(filter(best_conds, Experiment==X),   Layer1),
           Layer2 == pull(filter(best_conds, Experiment==X),   Layer2),
           Layer3 == pull(filter(best_conds, Experiment==X),   Layer3),
           Layers == pull(filter(best_conds, Experiment==X),   Layers),
           ImputedGenes == pull(filter(best_conds, Experiment==X), ImputedGenes),
           Experiment == X,
    )  
  
  return(out_df)
}

Exp0_df <- OptimalExtractor(data, X = 'FF') 
Exp1_df <- OptimalExtractor(data, X = 'ClassWeight') 
Exp2_df <- OptimalExtractor(data, X = 'deltaRate') 
Exp2a_df <- OptimalExtractor(data, X = 'Ablation deltaRate') 
Exp3_df <- OptimalExtractor(data, X = 'CuttleNet') 
Con_df <- OptimalExtractor(data, X = 'Control') 
Ct1_df <- OptimalExtractor(data, X = 'CellTICS 1 Layer') 
Ct5_df <- OptimalExtractor(data, X = 'CellTICS 5 Layer') 
Ct1c_df <- OptimalExtractor(data, X = 'CellTICS 1 Layer No 0s') 
Ct5c_df <- OptimalExtractor(data, X = 'CellTICS 5 Layer No 0s') 


# Ct_df <- rbind(Ct1c_df, Ct5c_df)
# Ct_df <- rbind(Ct1_df, Ct_df)
# Ct_df <- rbind(Ct5_df, Ct_df)
# Ct_df <- rbind(Exp3_df, Ct_df) 
# Ct_df <- rbind(Exp2_df, Ct_df) 


Ct_df <- rbind(Exp3_df, Ct5_df) 
Ct_df <- rbind(Ct1_df, Ct_df)


Ct_dfe <- Ct_df %>%
  mutate(Experiment = as.factor(Experiment)) %>%
  dplyr::select(Experiment, mean_F1, mean_F1_class) %>%
  rename(Subclass = mean_F1, 
         Class = mean_F1_class,
         Network = Experiment) %>%
  gather("Metric", "F1", -Network) 

fig_b <- ggplot(Ct_dfe, aes(x = Network, y = F1, fill = Metric)) +
  geom_boxplot() +
  theme_classic() +
  theme(legend.position = 'bottom')

# Perform ANOVA
anova_model <- aov(F1 ~ Network, data = filter(Ct_dfe, Metric == "Subclass"))
anova_summary <- summary(anova_model)

# Perform Tukey's HSD test
tukey_result <- TukeyHSD(anova_model)

# Extract the pairwise comparisons and p-values
tukey_df <- as.data.frame(tukey_result$Network)
tukey_df$comparison <- rownames(tukey_df)
tukey_df <- tukey_df %>%
  mutate(significance = case_when(
    `p adj` < 0.001 ~ "***",
    `p adj` < 0.01 ~ "**",
    `p adj` < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Create a list of comparisons and corresponding significance levels
comparison_list <- strsplit(as.character(tukey_df$comparison), "-")
annotations_list <- tukey_df$significance

# Determine y-position for the annotations
y_positions <- seq(from = max(Ct_dfe$F1) + 0.1, length.out = nrow(tukey_df), by = 0.05)

# Add statistical significance annotations to the plot
fig_b <- fig_b +
  geom_signif(
    comparisons = comparison_list,
    map_signif_level = TRUE,
    y_position = y_positions,
    annotations = annotations_list
  )


combined_plot <- fig_a + fig_b + plot_annotation(tag_levels = 'A')
combined_plot

combined_plot


ggsave('/home/sam/scRNAseq/Xenium/NeurIPS/NeurIPS_Fig3.pdf', combined_plot, device = "pdf", width = 8, height = 4)


anova_model <- aov(F1 ~ Network, data = filter(Ct_dfe, Metric == "Class"))
anova_summary <- summary(anova_model)

# Perform Tukey's HSD test
tukey_result <- TukeyHSD(anova_model)

# Extract the pairwise comparisons and p-values
tukey_df <- as.data.frame(tukey_result$Network)

##########################################################################################################

library(multcomp)

# Function to perform ANOVA and posthoc analysis
perform_analysis <- function(data, metric_name) {
  # Filter the dataframe for the specific metric
  filtered_data <- filter(data, Metric == metric_name)
  
  # Perform ANOVA
  anova_result <- aov(F1 ~ Network, data = filtered_data)
  print(paste("ANOVA for", metric_name))
  print(summary(anova_result))
  
  # Perform posthoc test if ANOVA is significant
  if (summary(anova_result)[[1]][["Pr(>F)"]][1] < 0.05) {
    posthoc_result <- TukeyHSD(anova_result)
    print(paste("Posthoc analysis for", metric_name))
    print(posthoc_result)
  } else {
    print(paste("No significant differences found for", metric_name))
  }
}

# Perform analysis for "Subclass"
perform_analysis(Ct_dfe, "Subclass")

# Perform analysis for "Class"
perform_analysis(Ct_dfe, "Class")

subclass_data <- Ct_dfe %>%
  filter(Metric == "Subclass", Network %in% c("deltaRate", "CuttleNet"))

# Compute the effect size (Cohen's d)
# This requires mean and standard deviation of each group
mean_deltaRate <- mean(subclass_data$F1[subclass_data$Network == "deltaRate"])
mean_CuttleNet <- mean(subclass_data$F1[subclass_data$Network == "CuttleNet"])
sd_deltaRate <- sd(subclass_data$F1[subclass_data$Network == "deltaRate"])
sd_CuttleNet <- sd(subclass_data$F1[subclass_data$Network == "CuttleNet"])

pooled_sd <- sqrt(((length(subclass_data$F1[subclass_data$Network == "deltaRate"]) - 1) * sd_deltaRate^2 + 
                     (length(subclass_data$F1[subclass_data$Network == "CuttleNet"]) - 1) * sd_CuttleNet^2) / 
                    (length(subclass_data$F1[subclass_data$Network == "deltaRate"]) + 
                       length(subclass_data$F1[subclass_data$Network == "CuttleNet"]) - 2))

effect_size <- (mean_deltaRate - mean_CuttleNet) / pooled_sd



###########################################################################################################



opt_Ct <- data %>%
  group_by(ImputedGenes) %>%
  OptimalExtractor(X = 'CellTICS')  %>%
  group_by(Experiment, Layers, ImputedGenes) %>%
  summarize(F1_sc = median(mean_F1),
            sd_sc = mean(sd_F1),
            F1_c = median(mean_F1_class),
            sd_c = mean(sd_F1_class),)

opt_Cut <- OptimalExtractor(data, X = 'CuttleNet')  %>%
  group_by(Experiment) %>%
  summarize(F1_sc = median(mean_F1),
            sd_sc = mean(sd_F1),
            F1_c = median(mean_F1_class),
            sd_c = mean(sd_F1_class),)

opt_Con <- OptimalExtractor(data, X = 'Control')  %>%
  group_by(Experiment) %>%
  summarize(F1_sc = median(mean_F1),
            sd_sc = mean(sd_F1),
            F1_c = median(mean_F1_class),
            sd_c = mean(sd_F1_class),)

opt_Cl <- OptimalExtractor(data, X = 'ClassWeight')  %>%
  group_by(Experiment) %>%
  summarize(F1_sc = median(mean_F1),
            sd_sc = mean(sd_F1),
            F1_c = median(mean_F1_class),
            sd_c = mean(sd_F1_class),)

opt_dr <- OptimalExtractor(data, X = 'deltaRate')  %>%
  group_by(Experiment) %>%
  summarize(F1_sc = median(mean_F1),
            sd_sc = mean(sd_F1),
            F1_c = median(mean_F1_class),
            sd_c = mean(sd_F1_class),)
