library(tidyverse)

# Function to read, process, and summarize Experiment 4 subclass data
process_experiment_subclass <- function(file) {
  read_csv(file) %>%
    mutate(Class = case_when(
      grepl("^\\d{2}_", Cluster) ~ "RGC",
      startsWith(Cluster, "AC_") ~ "AC",
      grepl("Photoreceptors$", Cluster) ~ "Ph",
      Cluster == "0MG (Mueller Glia)" ~ "MG",
      startsWith(Cluster, "0BC") ~ "BC",
      startsWith(Cluster, "0RBC") ~ "BC",
      TRUE ~ "Other")) %>%
    group_by(Class, Cluster, num_epochs, early_stopping, l1_lambda, arm_length,skip) %>%
    reframe(TPR = replace_na(TPR, 0),  # Replace NaN values with 0
            Prec = replace_na(Prec, 0),
            F1 = replace_na(F1, 0),
            cells = mean(cells),
            mean_TP = mean(TP),
            sd_TP = sd(TP),
            mean_TN = mean(TN),
            sd_TN = sd(TN),
            mean_FP = mean(FP),
            sd_FP = sd(FP),
            mean_FN = mean(FN),
            sd_FN = sd(FN),
            mean_prec = mean(Prec, na.rm = TRUE),
            sd_prec = sd(Prec, na.rm = TRUE),
            mean_tpr = mean(TPR, na.rm = TRUE),
            sd_tpr = sd(TPR, na.rm = TRUE),
            mean_F1 = mean(F1, na.rm = TRUE),
            sd_F1 = sd(F1, na.rm = TRUE)) %>%
    select(-TPR, -Prec, -F1) %>%
    ungroup() %>%
    unique() %>%
    group_by(num_epochs, early_stopping, l1_lambda, arm_length, skip) %>%
    mutate(cell_fraction = cells/sum(cells)) %>%
    ungroup()
}

# Function to read, process, and summarize Experiment 4 class data
process_experiment_class <- function(file) {
  read_csv(file) %>%
    select(-Accuracy, -TNR) %>%
    group_by(Class, num_epochs, early_stopping, l1_lambda, arm_length, skip) %>%
    reframe(cells_class = mean(cells),
              TPR = replace_na(TPR, 0),  # Replace NaN values with 0
              Prec = replace_na(Prec, 0),
              F1 = replace_na(F1, 0),
              cells_class = mean(cells),
              mean_TP_class = mean(TP),
              sd_TP_class = sd(TP),
              mean_TN_class = mean(TN),
              sd_TN_class = sd(TN),
              mean_FP_class = mean(FP),
              sd_FP_class = sd(FP),
              mean_FN_class = mean(FN),
              sd_FN_class = sd(FN),
              mean_prec_class = mean(Prec, na.rm = TRUE),
              sd_prec_class = sd(Prec, na.rm = TRUE),
              mean_tpr_class = mean(TPR, na.rm = TRUE),
              sd_tpr_class = sd(TPR, na.rm = TRUE),
              mean_F1_class = mean(F1, na.rm = TRUE),
              sd_F1_class = sd(F1, na.rm = TRUE)) %>%
    select(-TPR, -Prec, -F1) %>%
    ungroup() %>%
    unique() %>%
    group_by(num_epochs, early_stopping, l1_lambda, arm_length, skip) %>%
    mutate(cell_fraction_class = cells_class/sum(cells_class)) %>%
    ungroup()
}

# Path to Experiment  results
exp_path <- '/home/sam/scRNAseq/Xenium/AlonNN/NoiseInj/'

# Process Experiment  subclass and class data
exp_subclass_file <- paste0(exp_path, "experiment_results_perSubtype.csv")
exp_class_file <- paste0(exp_path, "experiment_results_perClass.csv")

exp_subclass_data <- process_experiment_subclass(exp_subclass_file)
exp_class_data <- process_experiment_class(exp_class_file)

# Merge subclass and class data from Experiment 
exp_combined_data <- left_join(exp_subclass_data, exp_class_data, 
                                by = c("Class", "num_epochs", 
                                       "early_stopping", "l1_lambda", 
                                       "arm_length", "skip")) %>%
  mutate(Experiment = "Deep", ClassWeight = NA, deltaRate = NA) %>%
  group_by(Class, num_epochs, early_stopping, l1_lambda, arm_length) %>%
  mutate(class_members = n())

# Export the merged dataset
write_csv(exp_combined_data, '/home/sam/scRNAseq/Xenium/AlonNN/NoiseInj/Merged_Experiments_DeepSkip.csv')

# Print success message
print("Merged dataset saved successfully.")
