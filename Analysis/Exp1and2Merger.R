library(tidyverse)

# Function to read CSV files
read_experiment_data <- function(file_path) {
  read_csv(file_path)
}

# Define paths to the directories
control_path <- '/home/sam/scRNAseq/Xenium/ClassVsSubclass/'
exp1_path <- '/home/sam/scRNAseq/Xenium/ClassVsSubclass/Experiment1_results/'
exp2_path <- '/home/sam/scRNAseq/Xenium/ClassVsSubclass/Experiment2_results/'

# Get list of CSV files from control and experiment directories
control_files <- list.files(control_path, pattern = "Control_Network_Results\\.csv$", full.names = TRUE)
exp1_files <- list.files(exp1_path, pattern = "\\.csv$", full.names = TRUE)
exp2_files <- list.files(exp2_path, pattern = "\\.csv$", full.names = TRUE)

# Read and combine data from Control, Experiment 1, and Experiment 2
control_data <- map_dfr(control_files, read_experiment_data)
exp1_data <- map_dfr(exp1_files, read_experiment_data)
exp2_data <- map_dfr(exp2_files, read_experiment_data)

# Combine data from Control, Experiment 1 and Experiment 2 into one dataframe
all_experiment_data <- bind_rows(control_data, exp1_data, exp2_data)

# Optional: Check the structure of the combined dataframe
glimpse(all_experiment_data)

write_csv(all_experiment_data, paste0(control_path, 'Exp1and2_Network_Results.csv'))
