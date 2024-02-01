library(tidyverse)

# Load the data for each experiment
all_data_path <- '/home/sam/scRNAseq/Xenium/ClassVsSubclass/Exp1and2_Network_Results.csv'
experiment1_path <- '/home/sam/scRNAseq/Xenium/ClassVsSubclass/Experiment1_results/Experiment1_AllClass_Network_results.csv'
experiment2_path <- '/home/sam/scRNAseq/Xenium/ClassVsSubclass/Experiment2_results/Experiment2_AllClass_Network_results.csv'

experiment1_data <- read_csv(experiment1_path)
experiment2_data <- read_csv(experiment2_path)
all_data <- read_csv(all_data_path)

df <- all_data %>%
  select(-matches(".*_class"), -class_cells, -class_cell_frac)

class_summary1 <- experiment1_data %>%
  select(-Accuracy, -TNR) %>%
  group_by(Class, ClassWeight) %>%
  reframe(cells_class = mean(cells),
         TPR = replace_na(TPR, 0),  # Replace NaN values with 0
         Prec = replace_na(Prec, 0),
         F1 = replace_na(F1, 0),
         mean_TP_class = mean(TP),
         sd_TP_class = sd(TP),
         mean_TN_class = mean(TN),
         sd_TN_class = sd(TN),
         mean_FP_class = mean(FP),
         sd_FP_class = sd(FP),
         mean_FN_class = mean(FN),
         sd_FN_class = sd(FN),
         mean_prec_class = mean(Prec),
         sd_prec_class = sd(Prec),
         mean_tpr_class = mean(TPR),
         sd_tpr_class = sd(TPR),
         mean_F1_class = mean(F1),
         sd_F1_class = sd(F1),
         deltaRate = 0) %>%
  select(-TPR, -Prec, -F1) %>%
  ungroup() %>%
  unique() %>%
  group_by(ClassWeight) %>%
  mutate(cell_fraction_class = cells_class/sum(cells_class)) %>%
  ungroup()

class_summary2 <- experiment2_data %>%
  select(-Accuracy, -TNR) %>%
  group_by(Class, deltaRate) %>%
  reframe(cells_class = mean(cells),
          TPR = replace_na(TPR, 0),  # Replace NaN values with 0
          Prec = replace_na(Prec, 0),
          F1 = replace_na(F1, 0),
          mean_TP_class = mean(TP),
          sd_TP_class = sd(TP),
          mean_TN_class = mean(TN),
          sd_TN_class = sd(TN),
          mean_FP_class = mean(FP),
          sd_FP_class = sd(FP),
          mean_FN_class = mean(FN),
          sd_FN_class = sd(FN),
          mean_prec_class = mean(Prec),
          sd_prec_class = sd(Prec),
          mean_tpr_class = mean(TPR),
          sd_tpr_class = sd(TPR),
          mean_F1_class = mean(F1),
          sd_F1_class = sd(F1),
          ClassWeight = NA) %>%
  select(-TPR, -Prec, -F1) %>%
  ungroup() %>%
  unique() %>%
  group_by(deltaRate) %>%
  mutate(cell_fraction_class = cells_class/sum(cells_class)) %>%
  ungroup()

experiments_df <- rbind(class_summary1, class_summary2)

final_df <- left_join(df, experiments_df, by = c("Class" , "ClassWeight" , "deltaRate" ))

write_csv(final_df, '/home/sam/scRNAseq/Xenium/ClassVsSubclass/Exp1and2_Network_Results_withClass.csv')
