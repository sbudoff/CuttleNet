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
    group_by(Class, Cluster, NumEpochs, EarlyStopping, L1Lambda) %>%
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
    group_by(NumEpochs, EarlyStopping, L1Lambda) %>%
    mutate(cell_fraction = cells/sum(cells)) %>%
    ungroup()
}

# Function to read, process, and summarize Experiment 4 class data
process_experiment_class <- function(file) {
  read_csv(file) %>%
    select(-Accuracy, -TNR) %>%
    group_by(Class, NumEpochs, EarlyStopping, L1Lambda) %>%
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
    group_by(NumEpochs, EarlyStopping, L1Lambda) %>%
    mutate(cell_fraction_class = cells_class/sum(cells_class)) %>%
    ungroup()
}

# Path to Experiment 4 results
exp4_path <- '/home/sam/scRNAseq/Xenium/ClassVsSubclass/Experiment4/'

# Process Experiment 4 subclass and class data
exp4_subclass_file <- paste0(exp4_path, "experiment_results_perSubtype.csv")
exp4_class_file <- paste0(exp4_path, "experiment_results_perClass.csv")

exp4_subclass_data <- process_experiment_subclass(exp4_subclass_file)
exp4_class_data <- process_experiment_class(exp4_class_file)

# Merge subclass and class data from Experiment 4
exp4_combined_data <- left_join(exp4_subclass_data, exp4_class_data, 
                                by = c("Class", "NumEpochs", 
                                       "EarlyStopping", "L1Lambda")) %>%
  mutate(Experiment = "CuttleNet", ClassWeight = NA, deltaRate = NA) %>%
  group_by(Class, NumEpochs, EarlyStopping, L1Lambda) %>%
  mutate(class_members = n())

# Load merged data from Experiments 1 and 2
exp12_data_path <- '/home/sam/scRNAseq/Xenium/ClassVsSubclass/Exp1and2_Network_Results_withClass.csv'
exp12_data <- read_csv(exp12_data_path) %>%
  mutate(NumEpochs = NA, EarlyStopping = NA, L1Lambda = NA)

# Merge the datasets
merged_data <- bind_rows(exp12_data, exp4_combined_data)

# Export the merged dataset
write_csv(merged_data, '/home/sam/scRNAseq/Xenium/ClassVsSubclass/Merged_Experiments_1_2_3.csv')

# Print success message
print("Merged dataset saved successfully.")
