library(tidyverse)
library(patchwork)


# ABLATION DATA, COMPARISON OF deltaRate Protocol to Control
# Path to Experiment 6 results
exp7_path <- '/home/sam/scRNAseq/Xenium/ClassVsSubclass/Experiment7/'

# Process Experiment 6 subclass and class data
exp7_subclass_file <- paste0(exp7_path, "experiment_results_perSubtype.csv")
exp7_class_file <- paste0(exp7_path, "experiment_results_perClass.csv")

df_c <- read_csv(exp7_class_file) %>%
  filter(l1_lambda == 0,
         Layer1 == 600) %>%
  dplyr::select(deltaRate, F1)



# Perform the t-test
t_test_result <- t.test(F1 ~ deltaRate, data = df_c)

# View the results
print(t_test_result)



# Perform the t-test
t_test_result <- t.test(F1 ~ deltaRate, data = df_c)

# View the results
print(t_test_result)

ablation_class_df <- read_csv(exp7_class_file) %>%
  filter(l1_lambda == 0,
         deltaRate==0,
         Layer1 == 600) %>%
  mutate(condition = "Ablation") %>%
  dplyr::select(condition, Class, F1, replicate)

ablation_sc_df <- read_csv(exp7_subclass_file) %>%
  filter(l1_lambda == 0,
         deltaRate==0,
         Layer1 == 600) %>%
  mutate(condition = "Ablation",
         Class = case_when(
           grepl("^\\d{2}_", Cluster) ~ "RGC",
           startsWith(Cluster, "AC_") ~ "AC",
           grepl("Photoreceptors$", Cluster) ~ "Ph",
           Cluster == "0MG (Mueller Glia)" ~ "MG",
           startsWith(Cluster, "0BC") ~ "BC",
           startsWith(Cluster, "0RBC") ~ "BC",
           TRUE ~ "Other")) %>%
  dplyr::select(condition, Class, Cluster, F1, replicate)
#################################################################################


# Full DATA, COMPARISON OF deltaRate Protocol to Control
# Path to Experiment 6 results
exp0_path <- '/home/sam/scRNAseq/Xenium/ClassVsSubclass/Experiment0/'

# Process Experiment 6 subclass and class data
exp0_subclass_file <- paste0(exp0_path, "experiment_results_perSubtype.csv")
exp0_class_file <- paste0(exp0_path, "experiment_results_perClass.csv")

imputed_class_df <- read_csv(exp0_class_file) %>%
  filter(l1_lambda == 0,
         layer1 == 250,
         layer2==0,
         layer3==0,
         early_stopping==0,
         num_epochs==50)  %>%
  mutate(condition = "Imputed") %>%
  dplyr::select(condition, Class, F1, replicate)


imputed_sc_df <- read_csv(exp0_subclass_file) %>%
  filter(l1_lambda == 0,
         layer1 == 250,
         layer2==0,
         layer3==0,
         early_stopping==0,
         num_epochs==50) %>%
  mutate(condition = "Imputed",
         Class = case_when(
             grepl("^\\d{2}_", Cluster) ~ "RGC",
             startsWith(Cluster, "AC_") ~ "AC",
             grepl("Photoreceptors$", Cluster) ~ "Ph",
             Cluster == "0MG (Mueller Glia)" ~ "MG",
             startsWith(Cluster, "0BC") ~ "BC",
             startsWith(Cluster, "0RBC") ~ "BC",
             TRUE ~ "Other")) %>%
  dplyr::select(condition, Class, Cluster, F1, replicate)

################################################################################
################## CLASS T-TEST ################################################

merged_df <- merge(imputed_class_df, ablation_class_df, by = "Class")

# Assuming 'F1_x' and 'F1_y' are the F1 scores from imputed_sc_df and ablation_sc_df respectively
# Perform the paired t-test
t_test_result <- t.test(merged_df$F1.x, merged_df$F1.y, paired = TRUE)


# View the results
print(t_test_result)

################################################################################
################ SUBCLASS T-TEST ###############################################

# Merge dataframes by cluster
merged_df <- merge(imputed_sc_df, ablation_sc_df, by = "Class")

# Assuming 'F1_x' and 'F1_y' are the F1 scores from imputed_sc_df and ablation_sc_df respectively
# Perform the paired t-test
t_test_result <- t.test(merged_df$F1.x, merged_df$F1.y, paired = TRUE)


# View the results
print(t_test_result)

