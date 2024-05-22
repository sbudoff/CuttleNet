library(tidyverse)
library(ggplot2)
library(ggsignif)
library(multcomp)

threshold <- 0.9

# Load the gradient_df and df_f1s dataframes
gradient_df <- read_csv('/home/sam/scRNAseq/Xenium/NeurIPS/All_Variants_GradientAlgo_thresh_summary.csv')
df_f1s <- read_csv('/home/sam/scRNAseq/Xenium/NeurIPS/All_Variants_AUCPR_thresh_summary.csv')

# Integrate gradient_df and df_f1s with summary_df
summary_df <- gradient_df %>%
  filter(Thresh == threshold, N_genes == 300) %>%
  dplyr::select(Seed, Counts) %>%
  rename("GraSP" = "Counts") 

# Integrate df_f1s into summary_df
df_f1s_long <- df_f1s %>%
  pivot_wider(names_from = Key, values_from = Count_Above_Threshold)

summary_df <- summary_df %>%
  left_join(df_f1s_long, by = "Seed") %>%
  mutate(GraSP = GraSP/130,
         Retina = Retina/130,
         BC = BC/17,
         AC = AC/63,
         RGC= RGC/45)

# Reshape the summary_df for plotting and statistical analysis
long_df <- summary_df %>%
  gather(key = "Method", value = "Subclasses", -Seed)

# Update the graph, ANOVA, and Tukey HSD test
method_order <- long_df %>%
  group_by(Method) %>%
  summarize(mean_subclasses = mean(Subclasses, na.rm = TRUE)) %>%
  arrange(desc(mean_subclasses)) %>%
  pull(Method)

long_df$Method <- factor(long_df$Method, levels = method_order)

plot <- ggplot(long_df, aes(x = Method, y = Subclasses)) +
  geom_boxplot() +
  theme_classic()

anova_model <- aov(Subclasses ~ Method, data = long_df)
tukey_result <- TukeyHSD(anova_model)

tukey_df <- as.data.frame(tukey_result$Method)
tukey_df$comparison <- rownames(tukey_df)
tukey_df <- tukey_df %>%
  mutate(significance = case_when(
    `p adj` < 0.001 ~ "***",
    `p adj` < 0.01 ~ "**",
    `p adj` < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Filter to only keep comparisons with "GraSP"
tukey_df <- tukey_df %>%
  filter(grepl("GraSP", comparison))

y_positions <- seq(from = max(long_df$Subclasses) + 0.05, length.out = nrow(tukey_df), by = 0.05)
comparison_list <- strsplit(as.character(tukey_df$comparison), "-")
annotations_list <- tukey_df$significance

plot <- plot +
  ylab('% of Possible Subclasses Classified above 0.9 TPR and Prec')+
  geom_signif(
    comparisons = comparison_list,
    map_signif_level = TRUE,
    y_position = y_positions,
    annotations = annotations_list
  )

print(plot)


# Save the plot as a PDF
ggsave("/home/sam/scRNAseq/Xenium/NeurIPS/supfig_MASTAUCPR.pdf", plot, width = 10, height = 6)
