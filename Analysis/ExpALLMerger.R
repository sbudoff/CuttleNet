library(tidyverse)

# Function to read, process, and summarize Experiment 4 subclass data
process_experiment_subclass <- function(file, group_cols) {
  # group_cols = c("Cluster", "Class", "NumEpochs", "EarlyStopping", 
  #                "L1Lambda", "Layer1", "Layer2", "Layer3")
  
  # Convert group_cols to symbols
  group_cols_syms <- syms(group_cols)
  
  read_csv(file) %>%
    mutate(Class = case_when(
      grepl("^\\d{2}_", Cluster) ~ "RGC",
      startsWith(Cluster, "AC_") ~ "AC",
      grepl("Photoreceptors$", Cluster) ~ "Ph",
      Cluster == "0MG (Mueller Glia)" ~ "MG",
      startsWith(Cluster, "0BC") ~ "BC",
      startsWith(Cluster, "0RBC") ~ "BC",
      TRUE ~ "Other")) %>%
    group_by(!!!group_cols_syms) %>%
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
    dplyr::select(-TPR, -Prec, -F1) %>%
    ungroup() %>%
    unique() 
}





# Function to read, process, and summarize Experiment 4 class data
process_experiment_class <- function(file, group_cols) {
  # 
  # group_cols = c("Class", "NumEpochs", "EarlyStopping", 
  #                "L1Lambda", "Layer1", "Layer2", "Layer3")
  group_cols_syms <- syms(group_cols)
  
  read_csv(file) %>%
    dplyr::select(-Accuracy, -TNR) %>%
    group_by(!!!group_cols_syms)  %>%
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
    dplyr::select(-TPR, -Prec, -F1) %>%
    ungroup() %>%
    unique() 
}

# Path to Experiment 4 results
exp0_path <- '/home/sam/scRNAseq/Xenium/ClassVsSubclass/Experiment0/'

# Process Experiment 4 subclass and class data
exp0_subclass_file <- paste0(exp0_path, "experiment_results_perSubtype.csv")
exp0_class_file <- paste0(exp0_path, "experiment_results_perClass.csv")

exp0_subclass_data <- process_experiment_subclass(exp0_subclass_file, group_cols = c("Class", "Cluster", "num_epochs", "early_stopping", 
                                                                                     "l1_lambda", "layer1", "layer2", "layer3"))%>%
  rename(NumEpochs = num_epochs, EarlyStopping=early_stopping,
         L1Lambda = l1_lambda, Layer1 = layer1,
         Layer2 = layer2, Layer3 = layer3
  )%>%
  group_by(NumEpochs, EarlyStopping, L1Lambda, Layer1, Layer2, Layer3) %>%
  mutate(cell_fraction = cells/sum(cells)) %>%
  ungroup()
exp0_class_data <- process_experiment_class(exp0_class_file, group_cols = c("Class", "num_epochs", "early_stopping", 
                                                                            "l1_lambda", "layer1", "layer2", "layer3")) %>%
  rename(NumEpochs = num_epochs, EarlyStopping=early_stopping,
         L1Lambda = l1_lambda, Layer1 = layer1,
         Layer2 = layer2, Layer3 = layer3
         )%>%
  group_by(NumEpochs, EarlyStopping, L1Lambda, Layer1, Layer2, Layer3) %>%
  mutate(cell_fraction_class = cells_class/sum(cells_class)) %>%
  ungroup()

# Merge subclass and class data from Experiment 4
exp0_combined_data <- left_join(exp0_subclass_data, exp0_class_data, 
                                  by = c("Class", "NumEpochs", 
                                         "EarlyStopping", "L1Lambda",
                                         "Layer1","Layer2","Layer3")) %>%
  mutate(Experiment = "FF", ClassWeight = NA, deltaRate = NA) %>%
  group_by(Class, NumEpochs, EarlyStopping, L1Lambda) %>%
  mutate(class_members = n())

# Load merged data from Experiments 1 and 2
exp123_data_path <- '/home/sam/scRNAseq/Xenium/ClassVsSubclass/Merged_Experiments_1_2_3.csv'
exp123_data <- read_csv(exp123_data_path) %>%
  mutate(Layer1 = NA, Layer2 = NA, Layer3 = NA)

# Load CellTICS
cellTICS_stats_path <- '/home/sam/scRNAseq/Xenium/ClassVsSubclass/Merged_CellTICS_stats.csv'
celltics_data <- read_csv(cellTICS_stats_path) 
  
# Merge the datasets
merged_data <- bind_rows(exp123_data, exp0_combined_data) %>%
  mutate(Layers=NA,
         ImputedGenes = ifelse(Experiment == "CuttleNet", F, T))

merged_data <- bind_rows(merged_data, celltics_data)

# Export the merged dataset
write_csv(merged_data, '/home/sam/scRNAseq/Xenium/ClassVsSubclass/Merged_Experiments_FINAL.csv')

# Print success message
print("Merged dataset saved successfully.")

names(merged_data)


# Path to Experiment 6 results
exp6_path <- '/home/sam/scRNAseq/Xenium/ClassVsSubclass/Experiment6/'

# Process Experiment 6 subclass and class data
exp6_subclass_file <- paste0(exp6_path, "experiment_results_perSubtype.csv")
exp6_class_file <- paste0(exp6_path, "experiment_results_perClass.csv")

exp6_subclass_data <- process_experiment_subclass(exp6_subclass_file, group_cols = c("Class", "Cluster", "num_epochs", "early_stopping", 
                                                                                     "l1_lambda", "deltaRate"))%>%
  rename(NumEpochs = num_epochs, EarlyStopping=early_stopping,
         L1Lambda = l1_lambda
  )%>%
  group_by(NumEpochs, EarlyStopping, L1Lambda, deltaRate) %>%
  mutate(cell_fraction = cells/sum(cells)) %>%
  ungroup()
exp6_class_data <- process_experiment_class(exp6_class_file, group_cols = c("Class", "num_epochs", "early_stopping", 
                                                                            "l1_lambda", "deltaRate")) %>%
  rename(NumEpochs = num_epochs, EarlyStopping=early_stopping,
         L1Lambda = l1_lambda
  )%>%
  group_by(NumEpochs, EarlyStopping, L1Lambda, deltaRate) %>%
  mutate(cell_fraction_class = cells_class/sum(cells_class)) %>%
  ungroup()

# Merge subclass and class data from Experiment 6
exp6_combined_data <- left_join(exp6_subclass_data, exp6_class_data, 
                                by = c("Class", "NumEpochs", 
                                       "EarlyStopping", "L1Lambda")) %>%
  mutate(Experiment = "CuttleNet deltaRate", ClassWeight = NA, 
         Layer1 = NA, Layer2 = NA, Layer3 = NA, ImputedGenes=F) %>%
  group_by(Class, NumEpochs, EarlyStopping, L1Lambda) %>%
  mutate(class_members = n())

merged_data <- bind_rows(merged_data, exp6_combined_data)




# Path to Experiment 6 results
exp7_path <- '/home/sam/scRNAseq/Xenium/ClassVsSubclass/Experiment7/'

# Process Experiment 6 subclass and class data
exp7_subclass_file <- paste0(exp7_path, "experiment_results_perSubtype.csv")
exp7_class_file <- paste0(exp7_path, "experiment_results_perClass.csv")

names(read_csv(exp7_class_file))

exp7_subclass_data <- process_experiment_subclass(exp7_subclass_file, group_cols = c("Class", "Cluster", "num_epochs", "early_stopping", 
                                                                                     "l1_lambda", "deltaRate", "Layer1"))%>%
  rename(NumEpochs = num_epochs, EarlyStopping=early_stopping,
         L1Lambda = l1_lambda
  )%>%
  group_by(NumEpochs, EarlyStopping, L1Lambda, deltaRate, Layer1) %>%
  mutate(cell_fraction = cells/sum(cells)) %>%
  ungroup()
exp7_class_data <- process_experiment_class(exp7_class_file, group_cols = c("Class", "num_epochs", "early_stopping", 
                                                                            "l1_lambda", "deltaRate", "Layer1")) %>%
  rename(NumEpochs = num_epochs, EarlyStopping=early_stopping,
         L1Lambda = l1_lambda
  )%>%
  group_by(NumEpochs, EarlyStopping, L1Lambda, deltaRate, Layer1) %>%
  mutate(cell_fraction_class = cells_class/sum(cells_class)) %>%
  ungroup()

# Merge subclass and class data from Experiment 6
exp7_combined_data <- left_join(exp7_subclass_data, exp7_class_data, 
                                by = c("Class", "NumEpochs", 
                                       "EarlyStopping", "L1Lambda", 'Layer1', 'deltaRate')) %>%
  mutate(Experiment = "Ablation", ClassWeight = NA, 
         Layer2 = NA, Layer3 = NA, ImputedGenes=F) %>%
  group_by(Class, NumEpochs, EarlyStopping, L1Lambda, deltaRate, Layer1) %>%
  mutate(class_members = n())

merged_data <- bind_rows(merged_data, exp7_combined_data)

# Export the merged dataset
write_csv(merged_data, '/home/sam/scRNAseq/Xenium/ClassVsSubclass/Merged_Experiments_FINAL.csv')

# Print success message
print("Merged dataset saved successfully.")

names(merged_data)

