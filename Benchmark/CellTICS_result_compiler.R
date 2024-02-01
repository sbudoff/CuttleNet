library(tidyverse)

CellTICS_Result_Compiler <- function(true_path, pred_path) {
  # Load the true labels
  true_labels <- read_csv(true_path)
  
  # Load the predictions
  predictions <- read_csv(pred_path)


  combined <- true_labels %>%
    rename(true_celltype = celltype, true_subcelltype = subcelltype) %>%
    bind_cols(predictions %>%
                rename(pred_celltype = celltype, pred_subcelltype = subcelltype))  
  
  class_results<- combined %>%
    group_by(true_celltype) %>%
    summarise(
      cells_class = n(),
      TP_class = sum(true_celltype == pred_celltype & true_subcelltype == pred_subcelltype),
      FP_class = sum(true_celltype != pred_celltype & true_subcelltype == pred_subcelltype),
      FN_class = sum(true_celltype == pred_celltype & true_subcelltype != pred_subcelltype),
      TN_class = sum(true_celltype != pred_celltype & true_subcelltype != pred_subcelltype)
     ) %>%
    mutate(
      TPR_class = TP_class / (TP_class + FN_class),
      Prec_class = TP_class / (TP_class + FP_class),
      F1_class = 2 * (Prec_class * TPR_class) / (Prec_class + TPR_class)
    )
  
  subclass_results<- combined %>%
    group_by(true_subcelltype, true_celltype) %>%
    summarise(
      cells = n(),
      TP = sum(true_subcelltype == pred_subcelltype & true_celltype == pred_celltype),
      FP = sum(true_subcelltype != pred_subcelltype & true_celltype == pred_celltype),
      FN = sum(true_subcelltype == pred_subcelltype & true_celltype != pred_celltype),
      TN = sum(true_subcelltype != pred_subcelltype & true_celltype != pred_celltype)) %>%
    mutate(
      TPR = TP / (TP + FN),
      Prec = TP / (TP + FP),
      F1 = 2 * (Prec * TPR) / (Prec + TPR)
    )
  
  results <- left_join(subclass_results, class_results, by = 'true_celltype') %>%
    rename(Class = true_celltype, 
           Cluster = true_subcelltype)
  
  return(results)
}


results_list = list()
i = 1
for (seed in seq(18,90,18)) {
  true_path = paste0("/home/sam/Poleg-Polsky/ICML/CellTICS/Retina/Retina_shuffle_",seed,"/retina_qlabel.csv")
  pred_path = paste0("/home/sam/Poleg-Polsky/ICML/CellTICS/Retina",seed,"_results/pred_y.csv")
  results_list[[i]] = CellTICS_Result_Compiler(true_path, pred_path) %>%
    mutate(Replicate = seed,
           Layers = 5,
           ImputedGenes = T)
  i = i + 1
}

results_5 <- bind_rows(results_list)


results_list = list()
i = 1
for (seed in seq(18,90,18)) {
  true_path = paste0("/home/sam/Poleg-Polsky/ICML/CellTICS/Retina/Retina_shuffle_",seed,"/retina_qlabel.csv")
  pred_path = paste0("/home/sam/Poleg-Polsky/ICML/CellTICS/Retina",seed,"_1Layer_results/pred_y.csv")
  results_list[[i]] = CellTICS_Result_Compiler(true_path, pred_path) %>%
    mutate(Replicate = seed,
           Layers = 1,
           ImputedGenes = T)
  i = i + 1
}

results_1 <- bind_rows(results_list)

results <- rbind(results_5, results_1)

results_list = list()
i = 1
for (seed in seq(18,90,18)) {
  true_path = paste0("/home/sam/Poleg-Polsky/ICML/CellTICS/Retina/Retina_shuffle_ablated_",seed,"/retina_qlabel.csv")
  pred_path = paste0("/home/sam/Poleg-Polsky/ICML/CellTICS/Retina",seed,"_Ablated_results/pred_y.csv")
  results_list[[i]] = CellTICS_Result_Compiler(true_path, pred_path) %>%
    mutate(Replicate = seed,
           Layers = 5,
           ImputedGenes = F)
  i = i + 1
}

results_5a <- bind_rows(results_list)
results <- rbind(results_5a, results)


results_list = list()
i = 1
for (seed in seq(18,90,18)) {
  true_path = paste0("/home/sam/Poleg-Polsky/ICML/CellTICS/Retina/Retina_shuffle_ablated_",seed,"/retina_qlabel.csv")
  pred_path = paste0("/home/sam/Poleg-Polsky/ICML/CellTICS/Retina",seed,"_Ablated_1layer_results/pred_y.csv")
  results_list[[i]] = CellTICS_Result_Compiler(true_path, pred_path) %>%
    mutate(Replicate = seed,
           Layers = 1,
           ImputedGenes = F)
  i = i + 1
}

results_1a <- bind_rows(results_list)
results <- rbind(results_1a, results)

# Function to read, process, and summarize Experiment 4 subclass data
results_stats <- results %>%
    group_by(Class, Cluster, Layers, ImputedGenes) %>%
    reframe(cells = mean(cells),
            TPR = replace_na(TPR, 0),  # Replace NaN values with 0
            Prec = replace_na(Prec, 0),
            F1 = replace_na(F1, 0),
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
    group_by(Layers) %>%
    mutate(cell_fraction = cells/sum(cells)) %>%
    ungroup()


# Function to read, process, and summarize Experiment 4 class data
results_stats_class <- results %>%
    group_by(Class, Layers, ImputedGenes) %>%
    reframe(cells_class = mean(cells_class),
            TPR_class = replace_na(TPR_class, 0),  # Replace NaN values with 0
            Prec_class = replace_na(Prec_class, 0),
            F1_class = replace_na(F1_class, 0),
            mean_TP_class = mean(TP_class),
            sd_TP_class = sd(TP_class),
            mean_TN_class = mean(TN_class),
            sd_TN_class = sd(TN_class),
            mean_FP_class = mean(FP_class),
            sd_FP_class = sd(FP_class),
            mean_FN_class = mean(FN_class),
            sd_FN_class = sd(FN_class),
            mean_prec_class = mean(Prec_class, na.rm = TRUE),
            sd_prec_class = sd(Prec_class, na.rm = TRUE),
            mean_tpr_class = mean(TPR_class, na.rm = TRUE),
            sd_tpr_class = sd(TPR_class, na.rm = TRUE),
            mean_F1_class = mean(F1_class, na.rm = TRUE),
            sd_F1_class = sd(F1_class, na.rm = TRUE)) %>%
    select(-TPR_class, -Prec_class, -F1_class) %>%
    ungroup() %>%
    unique() %>%
    group_by(Layers) %>%
    mutate(cell_fraction_class = cells_class/sum(cells_class)) %>%
    ungroup()

results_stats_final <- left_join(results_stats, results_stats_class, by =c('Class', "Layers", "ImputedGenes")) %>%
  mutate(Experiment = "CellTICS", 
         ClassWeight=NA, deltaRate=NA, NumEpochs=NA, EarlyStopping=NA, 
         L1Lambda=NA, Layer1=NA, Layer2=NA, Layer3=NA)


# exp123_data_path <- '/home/sam/scRNAseq/Xenium/ClassVsSubclass/Merged_Experiments_FINAL.csv'
# exp123_data <- read_csv(exp123_data_path)  %>%
#   mutate(Layers=0,
#          ImputedGenes = ifelse(Experiment == "CuttleNet", F, T))
# 
# 
# meta <- exp123_data %>%
#   select(Cluster, Class, class_members) %>%
#   unique()
# results_stats_final <- left_join(results_stats_final, meta, by = c('Class', 'Cluster'))

write_csv(results_stats_final, '/home/sam/scRNAseq/Xenium/ClassVsSubclass/Merged_CellTICS_stats.csv')


# # Merge the datasets
# merged_data <- bind_rows(exp123_data, results_stats_final)
# 
# # Export the merged dataset
# write_csv(merged_data, '/home/sam/scRNAseq/Xenium/ClassVsSubclass/Merged_Experiments_FINAL2.csv')



