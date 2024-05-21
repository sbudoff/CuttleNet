library(tidyverse)

data_path <- '/home/sam/scRNAseq/Xenium/NeurIPS/ML_df.RData'
load(data_path)

cluster <- unique(combined_df$Label)

combined_df <- combined_df %>%
  mutate(dataset = case_when(
    str_starts(Label, "AC_") ~ "dataset_1",
    str_detect(Label, "^[0-9]+_") ~ "dataset_2",
    TRUE ~ "dataset_3"
  ))

count_zeros <- combined_df %>%
  group_by(dataset) %>%
  summarize(across(everything(), ~ all(. == 0), .names = "all_zeros_{.col}")) %>%
  rowwise() %>%
  mutate(count_zero_columns = sum(c_across(starts_with("all_zeros_")))) %>%
  select(dataset, count_zero_columns)


summary_df <- combined_df %>%
  rename("Subclass" = "Label") %>%
  group_by(Subclass) %>%
  summarize(
    count = n(),
    proportion = round(n() / nrow(combined_df),4)
  ) %>%
  arrange(desc(count))


# Step 2: Create the histogram with ggplot2
plot <- ggplot(summary_df, aes(y = reorder(Subclass, count), x = count)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(trans = 'log10', sec.axis = sec_axis(~ ., name = "Proportion", labels = function(x) round(x / sum(summary_df$count), 4))) +
  labs(x = "Count", y = "Subclass") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 7)) +
  geom_text(aes(label = count), hjust = -0.1, size = 3)


# Step 3: Save the plot as a PDF using ggsave
ggsave("/home/sam/scRNAseq/Xenium/NeurIPS/cell_demographics.pdf", plot = plot, height = 13, width = 10)


write.csv(summary_df, "/home/sam/scRNAseq/Xenium/NeurIPS/cell_demographics.csv", row.names = FALSE)



  