library(tidyverse)
library(ggplot2)
library(ggsignif)
library(multcomp)

threshold <- 0.9

# Load the gradient_df and df_f1s dataframes
gradient_df <- read_csv('/home/sam/scRNAseq/Xenium/NeurIPS/All_Variants_GradientAlgo_thresh_summary.csv')
df_f1s <- read_csv('/home/sam/scRNAseq/Xenium/NeurIPS/All_Variants_AUCPR_thresh_summary.csv') %>%
  filter(Key == "Retina") %>%
  mutate(Key = "AUCPR")

# Function to extract Precision for each label from a confusion matrix
extract_precision <- function(cm) {
  precision <- cm$byClass[, "Pos Pred Value"]
  return(precision)
}

# Function to check if both values exceed the threshold
check_threshold <- function(tpr, prec) {
  return(tpr > threshold & prec > threshold)
}

# Function to extract TPR for each label from a confusion matrix
extract_tpr <- function(cm) {
  tpr <- cm$byClass[, "Sensitivity"]
  return(tpr)
}

compute_summary <- function(rez_list, threshold) {
  # Calculate TPR for each model
  tpr_dt <- extract_tpr(rez_list[['cMat_DT']])
  tpr_rf <- extract_tpr(rez_list[['cMat_rf']])
  tpr_svm <- extract_tpr(rez_list[['cMat_svm']])
  
  # Create a dataframe for each
  df_dt <- data.frame(Label=names(tpr_dt), DecisionTree=tpr_dt)
  df_rf <- data.frame(Label=names(tpr_rf), RandomForest=tpr_rf)
  df_svm <- data.frame(Label=names(tpr_svm), SVM=tpr_svm)
  
  # Merge all dataframes into one
  df_tpr <- Reduce(function(x, y) merge(x, y, by="Label", all=TRUE), list(df_dt, df_rf, df_svm))
  
  # Replace NA with 0 if any, because NA can appear if a label was not present in the test set of a particular model
  df_tpr[is.na(df_tpr)] <- 0
  
  # Calculate Precision for each model
  precision_dt <- extract_precision(rez_list[['cMat_DT']])
  precision_rf <- extract_precision(rez_list[['cMat_rf']])
  precision_svm <- extract_precision(rez_list[['cMat_svm']])
  
  # Create a dataframe for each
  df_dt <- data.frame(Label=names(precision_dt), DecisionTree=precision_dt)
  df_rf <- data.frame(Label=names(precision_rf), RandomForest=precision_rf)
  df_svm <- data.frame(Label=names(precision_svm), SVM=precision_svm)
  
  # Merge all dataframes into one
  df_precision <- Reduce(function(x, y) merge(x, y, by="Label", all=TRUE), list(df_dt, df_rf, df_svm))
  
  # Replace NA with 0 if any, because NA can appear if a label was not present in the test set of a particular model
  df_precision[is.na(df_precision)] <- 0
  
  # Match the labels and methods in both dataframes just to ensure they align correctly
  df_combined <- merge(df_tpr, df_precision, by = "Label", suffixes = c("_tpr", "_prec"))
  
  # Apply the function to each pair of TPR and Precision columns for each method
  df_combined$DecisionTree <- mapply(check_threshold, df_combined$DecisionTree_tpr, df_combined$DecisionTree_prec)
  df_combined$RandomForest <- mapply(check_threshold, df_combined$RandomForest_tpr, df_combined$RandomForest_prec)
  df_combined$SVM <- mapply(check_threshold, df_combined$SVM_tpr, df_combined$SVM_prec)
  
  # Select only the columns with the final boolean results
  df_final <- df_combined[, c("Label", "DecisionTree", "RandomForest", "SVM")]
  
  # Sum the boolean values to count labels above threshold for each ML model
  summary_counts <- data.frame(
    DecisionTree = sum(df_final$DecisionTree),
    RandomForest = sum(df_final$RandomForest),
    SVM = sum(df_final$SVM)
  )
  
  # Transpose the dataframe to match the requested format: models as rows, count as a column
  summary_counts <- t(summary_counts)
  colnames(summary_counts) <- c("Count")
  
  # Convert to a dataframe if needed (since transposing results in a matrix)
  summary_counts <- as.data.frame(summary_counts)
  
  out <- list(df_final, summary_counts)
  return(out)
}

ML_results <- list()
load('/home/sam/scRNAseq/Xenium/NeurIPS/ML_results.RData')

ML_results[["18"]] <- list()
ML_results[["18"]][['cMat_DT']] <- cMat_DT
ML_results[["18"]][['cMat_rf']] <- cMat_rf
ML_results[["18"]][['cMat_svm']] <- cMat_svm

rm(model_DT, model_rf, svm_model, cMat_DT, cMat_svm, cMat_rf)
gc()

# Update ML_results with other models (assuming this code section is correct)
for (i in seq(36, 90, 18)) {
  i_str <- paste(i)
  load(paste0('/home/sam/scRNAseq/Xenium/NeurIPS/ML_results_', i_str, '.RData'))
  
  if (!exists("cMat_rf")) {
    cMat_rf <- cMat_DT
    cMat_rf$table <- matrix(0, nrow = nrow(cMat_DT$table), ncol = ncol(cMat_DT$table))
  }
  
  if (!exists("cMat_svm")) {
    cMat_svm <- cMat_DT
    cMat_svm$table <- matrix(0, nrow = nrow(cMat_DT$table), ncol = ncol(cMat_DT$table))
  }
  
  ML_results[[i_str]] <- list()
  ML_results[[i_str]][['cMat_DT']] <- cMat_DT
  ML_results[[i_str]][['cMat_rf']] <- cMat_rf
  ML_results[[i_str]][['cMat_svm']] <- cMat_svm
  
  rm(model_DT, model_rf, svm_model, cMat_DT, cMat_svm, cMat_rf)
  gc()
}

# Compute summaries for each model
results_list <- list()
for (i in seq(18, 90, 18)) {
  i_str <- as.character(i)
  results_list[[i_str]] <- compute_summary(ML_results[[i_str]], threshold)
}

# Integrate gradient_df and df_f1s with summary_df
summary_df <- gradient_df %>%
  filter(Thresh == threshold, N_genes == 300) %>%
  dplyr::select(Seed, Counts) %>%
  rename("GraSP" = "Counts") %>%
  mutate(DecisionTree = 0, RandomForest = 0, SVM = 0)

for (i in seq(18, 90, 18)) {
  i_str <- as.character(i)
  result_row <- t(results_list[[i_str]][[2]])
  result_row_df <- as.data.frame(result_row)
  colnames(result_row_df) <- c("DecisionTree", "RandomForest", "SVM")
  row_index <- which(summary_df$Seed == i)
  summary_df[row_index, c("DecisionTree", "RandomForest", "SVM")] <- result_row_df
}

# Integrate df_f1s into summary_df
df_f1s_long <- df_f1s %>%
  pivot_wider(names_from = Key, values_from = Count_Above_Threshold)

summary_df <- summary_df %>%
  left_join(df_f1s_long, by = "Seed")

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
  theme_minimal()

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

y_positions <- seq(from = max(long_df$Subclasses) + 5, length.out = nrow(tukey_df), by = 5)
comparison_list <- strsplit(as.character(tukey_df$comparison), "-")
annotations_list <- tukey_df$significance

plot <- plot +
  geom_signif(
    comparisons = comparison_list,
    map_signif_level = TRUE,
    y_position = y_positions,
    annotations = annotations_list
  )

print(plot)


# Save the plot as a PDF
ggsave("/home/sam/scRNAseq/Xenium/NeurIPS/subclasses_boxplot.pdf", plot, width = 10, height = 6)

# Save the current results
save(summary_df, ML_results, gradient_df, plot, file = '/home/sam/scRNAseq/Xenium/NeurIPS/summary_results.RData')

# Extract ANOVA results
# Update the graph, ANOVA, and Tukey HSD test to include new columns from df_f1s
long_df <- summary_df %>%
  gather(key = "Method", value = "Subclasses", -Seed)

# Calculate mean subclasses for each method to determine the order
method_order <- long_df %>%
  group_by(Method) %>%
  summarize(mean_subclasses = mean(Subclasses, na.rm = TRUE)) %>%
  arrange(desc(mean_subclasses)) %>%
  pull(Method)

# Reorder Method factor levels based on the calculated order
long_df$Method <- factor(long_df$Method, levels = method_order)

# Create the boxplot with reordered x-axis
plot <- ggplot(long_df, aes(x = Method, y = Subclasses)) +
  geom_boxplot() +
  theme_minimal()

# Perform ANOVA and Tukey's HSD test
anova_model <- aov(Subclasses ~ Method, data = long_df)
tukey_result <- TukeyHSD(anova_model)

# Extract the pairwise comparisons and p-values
tukey_df <- as.data.frame(tukey_result$Method)
tukey_df$comparison <- rownames(tukey_df)
tukey_df <- tukey_df %>%
  mutate(significance = case_when(
    `p adj` < 0.001 ~ "***",
    `p adj` < 0.01 ~ "**",
    `p adj` < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Determine y-position for the annotations
y_positions <- seq(from = max(long_df$Subclasses) + 5, length.out = nrow(tukey_df), by = 5)

# Create a list of comparisons and corresponding significance levels
comparison_list <- strsplit(as.character(tukey_df$comparison), "-")
annotations_list <- tukey_df$significance

# Add statistical significance annotations with adjusted y-positions
plot <- plot +
  geom_signif(
    comparisons = comparison_list,
    map_signif_level = TRUE,
    y_position = y_positions,
    annotations = annotations_list
  )

# Display the plot
print(plot)

# Extract ANOVA results
anova_summary <- summary(anova_model)
f_value <- round(anova_summary[[1]][["F value"]][1], 2)
p_value <- anova_summary[[1]][["Pr(>F)"]][1]

# Extract Tukey results
tukey_df <- as.data.frame(tukey_result$Method)
graSP_vs_rf <- abs(round(tukey_df["RandomForest-GraSP", "diff"], 1))
graSP_vs_svm <- abs(round(tukey_df["SVM-GraSP", "diff"], 1))
graSP_vs_dt <- abs(round(tukey_df["DecisionTree-GraSP", "diff"], 1))
graSP_vs_rf_p <- round(tukey_df["RandomForest-GraSP", "p adj"], 3)
graSP_vs_svm_p <- round(tukey_df["SVM-GraSP", "p adj"], 6)
graSP_vs_dt_p <- round(tukey_df["DecisionTree-GraSP", "p adj"], 6)

# Add results for AC, BC, RGC, and Retina
graSP_vs_ac <- abs(round(tukey_df["AC-GraSP", "diff"], 1))
graSP_vs_bc <- abs(round(tukey_df["BC-GraSP", "diff"], 1))
graSP_vs_rgc <- abs(round(tukey_df["RGC-GraSP", "diff"], 1))
graSP_vs_retina <- abs(round(tukey_df["Retina-GraSP", "diff"], 1))
graSP_vs_ac_p <- round(tukey_df["AC-GraSP", "p adj"], 6)
graSP_vs_bc_p <- round(tukey_df["BC-GraSP", "p adj"], 6)
graSP_vs_rgc_p <- round(tukey_df["RGC-GraSP", "p adj"], 6)
graSP_vs_retina_p <- round(tukey_df["Retina-GraSP", "p adj"], 6)

# Create sentences using paste()
results_sentence1 <- paste("ANOVA results indicate a significant difference in the number of subclasses classified by the different methods (F(3, 16) =", 
                           f_value, ", p <", round(p_value, 3), 
                           "). Post-hoc Tukey's HSD test revealed that GraSP significantly outperformed SVM (p <", graSP_vs_svm_p, 
                           "), DecisionTree (p <", graSP_vs_dt_p, 
                           "), AC (p <", graSP_vs_ac_p, 
                           "), BC (p <", graSP_vs_bc_p, 
                           "), RGC (p <", graSP_vs_rgc_p, 
                           "), and Retina (p <", graSP_vs_retina_p, 
                           "), and showed a near-significant improvement over RandomForest (p =", graSP_vs_rf_p, ").")

results_sentence2 <- paste("Practically, GraSP demonstrated a substantial advantage in classifying subclasses with a stringent threshold of 0.9 for both TPR and Precision. Specifically, GraSP identified an average of", 
                           graSP_vs_rf, "more subclasses than RandomForest,", 
                           graSP_vs_svm, "more than SVM,", 
                           graSP_vs_dt, "more than DecisionTree,", 
                           graSP_vs_ac, "more than AC,", 
                           graSP_vs_bc, "more than BC,", 
                           graSP_vs_rgc, "more than RGC, and", 
                           graSP_vs_retina, "more than Retina, highlighting its superior performance compared to traditional machine learning methods.")

# Print sentences
cat(results_sentence1, "\n")
cat(results_sentence2, "\n")
