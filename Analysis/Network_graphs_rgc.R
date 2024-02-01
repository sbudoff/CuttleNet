library(tidyverse)

path = '/home/sam/scRNAseq/Xenium/curriculum_models_results_tran.csv'

df_tran <- read_csv(path) %>%
  filter(Genes == "Tran") %>%
  mutate(Genes = 57,
         Source = "Tran")

path = '/home/sam/scRNAseq/Xenium/curriculum_models_results.csv'

df <- read_csv(path)

# df <- df[grepl("^[0-9]{2}", df$Cluster), ]

df <- df %>%
  mutate(Source = "Net") %>%
  rbind(df_tran)

thresh <- 0.90

# Calculate the median and standard error of TPR and Prec for each Genes value
summary_df <- df %>%
  mutate(Prec = ifelse(is.na(Prec), 0, Prec))%>%
  group_by(Genes, Source) %>%
  summarise(
    median_TPR = median(TPR),
    SE_TPR = sd(TPR) / sqrt(n()),
    median_Prec = median(Prec, na.rm = T),
    SE_Prec = sd(Prec) / sqrt(n())
  )


# Filter the data for "Tran"
tran_data <- filter(summary_df, Source == "Tran")

# Create a plot with medians and error bars for both TPR and Prec
summary_df %>%
  filter(Source == "Net") %>%
  ggplot(aes(x = Genes)) +
  geom_point(aes(y = median_TPR, color = "TPR"), size = 3) + # Plot the median of TPR as blue points
  geom_errorbar(aes(ymin = median_TPR - SE_TPR, ymax = median_TPR + SE_TPR, color = "TPR"), width = 0.2) + # Error bars for TPR
  geom_point(aes(y = median_Prec, color = "Prec"), size = 3) + # Plot the median of Prec as red points
  geom_errorbar(aes(ymin = median_Prec - SE_Prec, ymax = median_Prec + SE_Prec, color = "Prec"), width = 0.2) + # Error bars for Prec
  labs(
    x = "Genes",
    y = "Median Value",
    title = "Median True Positive Rate and Precision For All Classes Per Number of Predictor Genes"
  ) +
  ylim(c(0,1.0)) +
  scale_color_manual(
    values = c("TPR" = "blue", "Prec" = "red"),
    labels = c("TPR", "Prec")
  ) + # Color mapping and legend labels
  theme_minimal() +
  theme(legend.title = element_blank(),
    text = element_text(size = 20) # Adjust the size to your preferred font size
  ) +
  geom_hline(data = tran_data, aes(yintercept= median_TPR), color = "red", alpha = 0.5) +
  geom_hline(data = tran_data, aes(yintercept= median_TPR-SE_TPR), color = "red", alpha = 0.5, linetype="dashed") +
  geom_hline(data = tran_data, aes(yintercept= median_TPR+SE_TPR), color = "red", alpha = 0.5, linetype="dashed") +
  geom_hline(data = tran_data, aes(yintercept= median_Prec), color = "blue", alpha = 0.5) +
  geom_hline(data = tran_data, aes(yintercept= median_Prec-SE_Prec), color = "blue", alpha = 0.5, linetype="dashed") +
  geom_hline(data = tran_data, aes(yintercept= median_Prec+SE_Prec), color = "blue", alpha = 0.5, linetype="dashed") 



df_final <- df %>%
  mutate(TPR = TPR > thresh,
         Prec = Prec > thresh,
         Identifiable = 1-((TPR+Prec) == 0),
         Drop = ifelse(Source=="Net" & Genes < 150, 1, 0),
         Source = ifelse(Source=="Net", "Gradient", "Tran 2019")) %>%
  filter(Drop == 0) %>%
  group_by(Genes, Source) %>%
  summarize(Subtypes = sum(Identifiable)) %>%
  ungroup() %>%
  rbind(data.frame(Subtypes = c(16,48,125), Genes = c(10,39,18222), Source = c("Manual", "Gradient", "scRNAseq")))

final_plot <- ggplot(df_final, aes(x= Genes, y = Subtypes, color = Source)) +
  # ggtitle("Retinal Subtypes Classifiable With Proposed Gene Sets") +
  ylim(1,125) +
  geom_rect(
    aes(xmin = 3, xmax = 40, ymin = -Inf, ymax = Inf),
    fill = "lightblue", alpha = 0.05,
  color = NA  
  )+
  geom_rect(
    aes(xmin = 150, xmax = 300, ymin = -Inf, ymax = Inf),
    fill = "lightgreen", alpha = 0.05,
    color = NA  
  ) +
  geom_rect(
    aes(xmin = 300, xmax = 500, ymin = -Inf, ymax = Inf),
    fill = "yellow", alpha = 0.05,
    color = NA  
  ) +
  geom_rect(
    aes(xmin = 500, xmax = 1000, ymin = -Inf, ymax = Inf),
    fill = "orange", alpha = 0.05,
    color = NA  
  ) +
  geom_rect(
    aes(xmin = 5000, xmax = 50000, ymin = -Inf, ymax = Inf),
    fill = "red", alpha = 0.05,
    color = NA  
  ) +
  geom_text(
    aes(x = 10, y = 125, label = "RNAscope"),
    size = 4, color = "black", hjust = 0.5, family = "Times"
  ) +
  geom_text(
    aes(x = 225, y = 122, label = "10x Xenium"),
    size = 4, color = "black", hjust = 0.5, family = "Times"
  ) +
  geom_text(
    aes(x = 400, y = 116, label = "Vizgen MERFISH"),
    size = 4, color = "black", hjust = 0.5, family = "Times"
  ) +
  geom_text(
    aes(x = 700, y = 110, label = "Nanostring CosMx"),
    size = 4, color = "black", hjust = 0.5, family = "Times"
  ) + 
  geom_text(
    aes(x = 15000, y = 110, label = "scRNAseq"),
    size = 4, color = "black", hjust = 0.5, , family = "Times"
  ) +
  geom_point(size=5) +
  scale_x_log10() +
  theme_classic() +
  theme(
    text = element_text(size = 20, family = "Times"),
    axis.title = element_text(size = 20, family = "Times"),
    axis.text = element_text(size = 20, family = "Times"),
    legend.title = element_text(size = 20, family = "Times"),
    legend.text = element_text(size = 20, family = "Times"),
    plot.title = element_text(size = 20, family = "Times")
  )

# Specify the filename and file type (SVG)
filename <- "/home/sam/Poleg-Polsky/ERM/Presentation Images/visium_plot.svg"
icml_path <- '/home/sam/Poleg-Polsky/ICML/Figures/Fig2.pdf'

# Export the ggplot to an SVG file
ggsave(filename, plot = final_plot, device = "svg")
ggsave(icml_path, plot = final_plot, device = "pdf", width = 10, height = 6, dpi = 300)


classified_df <- df %>%
  mutate(TPR = TPR > thresh,
         Prec = Prec > thresh,
         Ambiguous = (TPR+Prec) == 0) %>%
  group_by(Genes, Source) %>%
  summarize(
    `Ambiguous Subtypes` = sum(Ambiguous),
    `Ambiguous Cells` = sum(cells*Ambiguous),
    `Percent Ambiguous Subtypes` = sum(Ambiguous)/n(),
    `Percent Ambiguous Cells` = sum(cells*Ambiguous)/sum(cells)) 

# Filter the data for "Tran"
tran_data <- filter(classified_df, Source == "Tran")

# Create a dual-axis plot
classified_df%>%
  filter(Source == "Net") %>%
  ggplot(aes(x = Genes)) +
    geom_line(aes(y = 1-`Percent Ambiguous Subtypes`, group = 1, color = "Subtypes")) +
    geom_line(aes(y = 1-`Percent Ambiguous Cells`, group = 1, color = "Cells")) +
    labs(
      x = "Genes",
      y = "Classified Cells",
      color = "Legend",
      fill = "Legend",
      title = paste0("Portion of Cells and Subtypes Classified with Greater than ", 100*thresh, "% TPR or Precision")
    ) +
    scale_color_manual(values = c("Subtypes" = "blue", "Cells" = "red")) +
    theme_minimal() +
    theme(
      axis.text.y.right = element_text(color = "blue"),
      axis.title.y.right = element_text(color = "blue"),
      legend.title = element_blank(),
      text = element_text(size = 20),
      legend.position = 'bottom'
    ) +
    ylim(c(0,1.0))  +
  geom_hline(data = tran_data, aes(yintercept= `Percent Ambiguous Cells`), color = "red", alpha = 0.5) +
  geom_hline(data = tran_data, aes(yintercept=  `Percent Ambiguous Subtypes`), color = "blue", alpha = 0.5) 

