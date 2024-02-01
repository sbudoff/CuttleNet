library(tidyverse)

# Function to read and process each file
process_file <- function(file) {
  read_csv(file) %>%
    mutate(experiment_id = basename(file))
}

path <- '/home/sam/scRNAseq/Xenium/ClassVsSubclass/Experiment2/'

control_data <- read_csv('/home/sam/scRNAseq/Xenium/ClassVsSubclass/Control_Network_Results.csv')

# List all CSV files in the directory
files <- list.files(path, pattern = "*.csv", full.names = TRUE)

# Read and combine all datasets
combined_data <- map_dfr(files, process_file)

# Extract classweight and replicateID from experiment_id
combined_data <- combined_data %>%
  extract(col = experiment_id,
          into = c("deltaRate", "replicateID"),
          regex = ".*_deltaRate_(\\d+\\.?\\d*)_result_replicate(\\d+).*",
          remove = TRUE,  # Set to TRUE if you want to remove the original column
          convert = TRUE)  %>%# Convert to appropriate types (numeric, integer, etc.)
  mutate(Prec = replace_na(Prec, 0),
       F1 = replace_na(F1, 0))


# View the structure of the combined data
glimpse(combined_data)


# Function to analyze trends and create scatterplots
analyze_trends <- function(data, score_column) {
  # Summarize data
  summary_data <- data %>%
    select(Cluster, deltaRate, replicateID, !!sym(score_column)) %>%
    group_by(Cluster, deltaRate) %>%
    summarize(mean_score = mean(!!sym(score_column), na.rm = TRUE),
              sd_score = sd(!!sym(score_column), na.rm = TRUE),
              .groups = 'drop')
  return(summary_data)
}

# Use the function to create summary
f1_score_df <- analyze_trends(combined_data, "F1")

plot_trends <- function(score_df, metric = "F1") {
  plot_list <- list()
  for (cluster in unique(score_df$Cluster)) {
    # Filter data for the current cluster
    cluster_data <- score_df %>% filter(Cluster == cluster)
    
    # Data for deltaRate = 0
    zero_data <- cluster_data %>% filter(deltaRate == 0)
    
    # Create the plot
    plot <- cluster_data %>%
      filter(deltaRate != 0) %>%
      ggplot(aes(x = deltaRate, y = mean_score)) +
      annotate("rect", xmin = 0.0008, xmax = 1.01, ymin = mean(zero_data$mean_score) - mean(zero_data$sd_score),
               ymax = mean(zero_data$mean_score) + mean(zero_data$sd_score), alpha = 0.5, fill = "grey") +
      geom_point() +
      geom_errorbar(aes(ymin = mean_score - sd_score, ymax = mean_score + sd_score), width = 0.1) +
      labs(title = paste(metric, " Score Trends for Cluster", cluster),
           x = "Class Weight",
           y = "Mean", metric, "Score") +
      scale_x_log10() +
      theme_minimal() +
      # Add horizontal line and shaded region for classweight = 0
      geom_hline(yintercept = mean(zero_data$mean_score), linetype = "dotted") 
      
    
    plot_list[[cluster]] <- plot
  }
  return(plot_list)
}

f1_score_plots <- plot_trends(f1_score_df)

# Display plots (example for one cluster, replace 'cluster_name' with actual cluster name)
for(clust in names(f1_score_plots)) {
  print(f1_score_plots[[clust]])
}

f1_interesting_results <- c("AC_56", "AC_38", "AC_33", "30_Novel")

# Use the function to create summary
TPR_score_df <- analyze_trends(combined_data, "TPR")
TPR_score_plots <- plot_trends(TPR_score_df, metric = "TPR")

# Display plots (example for one cluster, replace 'cluster_name' with actual cluster name)
for(clust in names(TPR_score_plots)) {
  print(TPR_score_plots[[clust]])
}

tpr_interesting_results <- c("AC_56", "AC_51", "AC_25", "AC_24", "AC_22", "AC_21", 
                             "AC_2", "AC_19", "AC_18", "AC_17", "AC_16", "AC_25",  
                             "AC_14", "AC_13", "AC_11", "AC_10", "AC_1")


# Use the function to create summary
Prec_score_df <- analyze_trends(combined_data, "Prec")
Prec_score_plots <- plot_trends(Prec_score_df, metric = "Prec")

# Display plots (example for one cluster, replace 'cluster_name' with actual cluster name)
for(clust in names(Prec_score_plots)) {
  print(Prec_score_plots[[clust]])
}

prec_interesting_results <- c("AC_9", "AC_8", "AC_7", "AC_61", "AC_60", "AC_6",
                              "AC_58", "AC_56", "AC_54", "AC_53", "AC_5", "AC_49",
                              "AC_47", "AC_42", "AC_41", "AC_40", "AC_38", "AC_34",
                              "AC_33", "AC_32", "AC_29", "AC_26", "43_AlphaONS",
                              "40_M1dup", "38_FmidiON", "37_Novel", "36_Novel",
                              "35_Novel", "34_Novel", "31_M2", "30_Novel", 
                              "27_Novel", "26_Novel", "24_Novel", "23_W3D2",
                              "22_M5", "21_Tbr1_S2", "20_Novel", "15_Novel", 
                              "14_ooDS_Cck", "13_Novel", "12_ooDS_NT", "11_Novel")
                              
                    
# Function to perform ANOVA for each cluster and return the results
metric_anova <- function(data, metric) {
  
  data_f <- data %>% filter(classweight!=1)
  anova_results <- list()
  
  for (cluster in unique(data$Cluster)) {
    data_cluster <- data_f %>% 
      filter(Cluster == cluster) %>%
      mutate(score = !!sym(metric),
             classweight = as.factor(classweight)) 

    anova_test <- aov(score ~ classweight, data = data_cluster)
    anova_summary <- broom::tidy(anova_test)
    anova_summary$Cluster <- cluster
    anova_results[[cluster]] <- anova_summary
  }
  
  # Combine results into a dataframe
  results_df <- bind_rows(anova_results) %>%
    filter(term == "classweight") %>%
    dplyr::select(-term) %>%
    mutate(Metric = metric)
  return(results_df)
}


post_hoc_analysis <- function(anova_data, original_data, metric, alpha=0.05) {
  # Filter significant ANOVA results
  significant_clusters <- anova_data %>% 
    filter(p.value < alpha) %>%
    pull(Cluster)
  
  tukey_results <- list()
  
  for (cluster in significant_clusters) {
    data_cluster <- original_data %>%
      filter(Cluster == cluster, classweight != 1) %>%
      mutate(score = !!sym(metric),
             classweight = as.factor(classweight))
    
    # Perform Tukey HSD test
    tukey_test <- TukeyHSD(aov(score ~ classweight, data = data_cluster))
    
    # Extract and format the Tukey HSD results
    tukey_summary <- data.frame(tukey_test$`classweight`)
    tukey_summary <- tukey_summary %>%
      rownames_to_column(var = "Comparison") %>%
      separate(Comparison, into = c("classweight", "comparison"), sep = "-") %>%
      filter(comparison == 0) %>%
      dplyr::select(classweight, p.adj) %>%
      filter(p.adj <= alpha)
    if (nrow(tukey_summary > 0)) {
      tukey_summary$Cluster <- cluster
      tukey_summary$metric <- metric
      tukey_results[[cluster]] <- tukey_summary
    }
  }
  
  # Combine results into a dataframe
  results_df <- bind_rows(tukey_results)
  return(results_df)
}

append_mean_differences <- function(tukey_df, original_data, metric) {
  # Add a new column for mean differences
  tukey_df$mean_diff <- NA
  tukey_df$mean_score0 <- NA
  tukey_df$mean_score1 <- NA
  
  # Iterate through each row in tukey_df
  for(i in 1:nrow(tukey_df)) {
    # Extract relevant details from the current row
    current_cluster <- tukey_df$Cluster[i]
    current_classweight <- as.numeric(tukey_df$classweight[i])
    
    # Filter original_data for the current cluster and classweights
    data_cluster <- original_data %>%
      filter(Cluster == current_cluster, 
             classweight %in% c(0, current_classweight)) %>%
      mutate(score = !!sym(metric))
    
    # Summarize mean for each classweight
    summary_data <- data_cluster %>%
      group_by(classweight) %>%
      summarize(mean_score = mean(score, na.rm = TRUE), .groups = 'drop')
    
    # Calculate and store the mean difference
    mean_diff <- summary_data$mean_score[1] - summary_data$mean_score[2]
    tukey_df$mean_diff[i] <- mean_diff
    tukey_df$mean_score0[i] <- summary_data$mean_score[1]
    tukey_df$mean_score1[i] <- summary_data$mean_score[2]
  }
  
  return(tukey_df)
}

# Statistical testing of F1
f1_anova <- metric_anova(combined_data, "F1")
f1_tukey_results <- post_hoc_analysis(f1_anova, combined_data, "F1")
tukey_with_means_f1 <- append_mean_differences(tukey_df = f1_tukey_results, 
                                            original_data = combined_data, 
                                            metric = "F1")
f1_interesting <- tukey_with_means_f1 %>%
  filter(mean_diff < 0)
  

# Statistical testing of TPR
TPR_anova <- metric_anova(combined_data, "TPR")
TPR_tukey_results <- post_hoc_analysis(TPR_anova, combined_data, "TPR")
tukey_with_means_TPR <- append_mean_differences(tukey_df = TPR_tukey_results, 
                                               original_data = combined_data, 
                                               metric = "TPR")
TPR_interesting <- tukey_with_means_TPR %>%
  filter(mean_diff < 0)

# Statistical testing of Prec
prec_anova <- metric_anova(combined_data, "Prec")
Prec_tukey_results <- post_hoc_analysis(prec_anova, combined_data, "Prec")
tukey_with_means_Prec <- append_mean_differences(tukey_df = Prec_tukey_results, 
                                               original_data = combined_data, 
                                               metric = "Prec")
Prec_interesting <- tukey_with_means_Prec %>%
  filter(mean_diff < 0)


cell_numbers <- combined_data %>%
  select(Cluster, cells) %>%
  unique()

combined_tukey_df <- rbind(tukey_with_means_Prec, tukey_with_means_TPR, tukey_with_means_f1) %>%
  left_join(cell_numbers, by = "Cluster") %>%
  mutate(Interpretation = ifelse(mean_diff < 0, "Improvement", "Degradation"),
         Interpretation = as.factor(Interpretation))

ggplot(combined_tukey_df, aes(x=cells, y = mean_diff, shape = Interpretation, color = classweight)) +
  geom_point() +
  theme_minimal() +
  scale_x_log10() +
  facet_wrap(~metric)
