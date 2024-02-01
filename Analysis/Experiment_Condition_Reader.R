library(tidyverse)

# Function to read and process each file
process_file <- function(file) {
  read_csv(file) %>%
    mutate(experiment_id = basename(file))
}

Experiment_Reader <- function(experiment, condition,
                              out_path, in_path,
                              visualize = T) {

  files <- list.files(in_path, pattern = paste0(experiment, "_result_replicate.*\\.csv$"), full.names = TRUE)
  data <- map_dfr(files, process_file) %>%
    mutate(Class = case_when(
      grepl("^\\d{2}_", Cluster) ~ "RGC",
      startsWith(Cluster, "AC_") ~ "AC",
      grepl("Photoreceptors$", Cluster) ~ "Ph",
      Cluster == "0MG (Mueller Glia)" ~ "MG",
      startsWith(Cluster, "0BC") ~ "BC",
      startsWith(Cluster, "0RBC") ~ "BC",
      TRUE ~ "Other"
    ),
    Class = as.factor(Class))
  
  summary <- data %>%
    mutate(Prec = replace_na(Prec, 0),
           F1 = replace_na(F1, 0)) %>%
    group_by(Class, Cluster) %>%
    summarize(cells = mean(cells),
              mean_TP = mean(TP),
              sd_TP = sd(TP),
              mean_TN = mean(TN),
              sd_TN = sd(TN),
              mean_FP = mean(FP),
              sd_FP = sd(FP),
              mean_FN = mean(FN),
              sd_FN = sd(FN),
              mean_prec = mean(Prec),
              sd_prec = sd(Prec),
              mean_tpr = mean(TPR),
              sd_tpr = sd(TPR),
              mean_F1 = mean(F1),
              sd_F1 = sd(F1)) %>%
    ungroup()  %>%
    mutate(cell_fraction = cells/sum(cells)) %>%
    group_by(Class) %>%
    mutate(class_members = n(),
           class_cells = sum(cells)) %>%
    ungroup() %>%
    mutate(class_cell_frac = class_cells/sum(cells))
  
  class_summary <- data %>%
    select(-Cluster, -Accuracy, -TNR) %>%
    group_by(Class, experiment_id) %>%
    mutate(cells = sum(cells),
           TP = sum(TP),
           TN = sum(TN),
           FP = sum(FP),
           FN = sum(FN),
           TPR = TP / (TP + FN),  # True Positive Rate
           Prec = TP / (TP + FP),  # Precision
           F1 = 2 * ((Prec * TPR) / (Prec + TPR)),  # F1 Score
           TPR = replace_na(TPR, 0),  # Replace NaN values with 0
           Prec = replace_na(Prec, 0),
           F1 = replace_na(F1, 0)
    ) %>%
    ungroup() %>%
    unique() %>%
    group_by(Class) %>%
    summarize(class_cells = mean(cells),
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
              sd_F1_class = sd(F1)) %>%
    ungroup() 
  
  data <- left_join(summary, class_summary, by=c("Class", "class_cells"))
  
  if (visualize) {
    ggplot(data, aes(x=cell_fraction, y = mean_F1, color = Class))+
      geom_point() +
      geom_errorbar(aes(ymin = mean_F1 - sd_F1, ymax = mean_F1 + sd_F1), width = 0.01) +
      ggtitle(paste(condition, "Network")) +
      scale_x_log10()+
      theme_minimal()
  }
  
  write_csv(data, paste0(out_path, condition,'_Network_Results.csv'))
  
  return(data)
}

control_data <- Experiment_Reader(experiment = "Experiment1_classVSsubclass_classweight_0",
                                  condition = "control",
                                  out_path = '/home/sam/scRNAseq/Xenium/ClassVsSubclass/',
                                  in_path = '/home/sam/scRNAseq/Xenium/ClassVsSubclass/Experiment1/',
                                  visualize = T) %>%
  mutate(Experiment = "Control",
         deltaRate = 0,
         ClassWeight = 0
         )

write_csv(control_data, paste0('/home/sam/scRNAseq/Xenium/ClassVsSubclass/', 'Control_Network_Results.csv'))

################################################################################
# Process Experiment 1
# Define the path to the directory
in_path <- '/home/sam/scRNAseq/Xenium/ClassVsSubclass/Experiment1/'
out_path <- '/home/sam/scRNAseq/Xenium/ClassVsSubclass/Experiment1_results/'

# List all files in the directory and extract unique experiment conditions
files <- list.files(in_path, pattern = "Experiment1_classVSsubclass_classweight_.*_result_replicate.*\\.csv$")
unique_conditions <- unique(sub("\\_result_replicate.*", "", sub(".*classweight_", "", files)))
unique_conditions <- unique_conditions[-1]
print(unique_conditions)

# Iterate over each unique condition and process data
for (condition in unique_conditions) {
  experiment <- paste0("Experiment1_classVSsubclass_classweight_", condition)
  
  # Call Experiment_Reader function
  condition_data <- Experiment_Reader(experiment = experiment,
                                      condition = paste0("ClassWeight", condition),
                                      out_path = out_path,
                                      in_path = in_path,
                                      visualize = TRUE) %>%
    mutate(Experiment = "ClassWeight",
           deltaRate = 0,
           ClassWeight = condition)
  write_csv(condition_data, paste0(out_path, paste0("ClassWeight", condition),'_Network_Results.csv'))
  
  print(paste("Processed condition:", condition))
}

################################################################################
# Process Experiment 2
# Define the path to the directory
in_path <- '/home/sam/scRNAseq/Xenium/ClassVsSubclass/Experiment2/'
out_path <- '/home/sam/scRNAseq/Xenium/ClassVsSubclass/Experiment2_results/'

# List all files in the directory and extract unique experiment conditions
files <- list.files(in_path, pattern = "Experiment2_classVSsubclass_deltaRate_.*_result_replicate.*\\.csv$")
unique_conditions <- unique(sub("\\_result_replicate.*", "", sub(".*deltaRate_", "", files)))
print(unique_conditions)

# Iterate over each unique condition and process data
for (condition in unique_conditions) {
  experiment <- paste0("Experiment2_classVSsubclass_deltaRate_", condition)
  
  # Call Experiment_Reader function
  condition_data <- Experiment_Reader(experiment = experiment,
                                      condition = paste0("deltaRate", condition),
                                      out_path = out_path,
                                      in_path = in_path,
                                      visualize = TRUE) %>%
    mutate(Experiment = "deltaRate",
           deltaRate = condition,
           ClassWeight = NA)
  
  write_csv(condition_data, paste0(out_path, paste0("deltaRate", condition),'_Network_Results.csv'))
  
  print(paste("Processed condition:", condition))
}
