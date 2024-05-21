library(tidyverse)

setwd('/home/sam/scRNAseq/FreqTables/')
load(file='Full_Retina_FreqTable50.RData')

# # Find all RIKEN IDs and fix them
# test <- grep("Rik", genes, value=T)



Retina_freq_table_expressingGenes <- Retina_freq_table[0:0,]
bins <- unique(Retina_freq_table$ExpressionLevel)
low_cut_oof = 3
records_0 = nrow(Retina_freq_table)
for (cell_type in unique(Retina_freq_table$CellType)) 
{
  drop_genes = c()
  n_classes = length(unique(filter(Retina_freq_table, CellType==cell_type)$CellClass))
  for (gene in unique(Retina_freq_table$Gene))
  {
    test <- Retina_freq_table %>%
      filter(CellType == cell_type) %>%
      filter(Gene == gene, ExpressionLevel <= bins[low_cut_oof]) %>%
      group_by(CellClass) %>%
      summarise(inBin = round(sum(NormalizedFreq),2)) %>%
      pull(inBin)
    if (sum(round(test,0)) == n_classes) 
    {
      drop_genes = c(drop_genes, gene)
      print(sprintf("Dropping %s", gene))
    }
  }
  print(sprintf("Completed all %s cells", cell_type))
  temp <- Retina_freq_table %>%
    filter(CellType == cell_type) %>%
    filter(!Gene %in% drop_genes)
  Retina_freq_table_expressingGenes <- rbind(Retina_freq_table_expressingGenes, temp)
  remove(temp)
  print(sprintf("%s genes removed from the %s set",length(drop_genes), cell_type))
}
records_1 = nrow(Retina_freq_table_expressingGenes)
drop_perc = round((records_0-records_1)/records_0,3)*100
print(sprintf('%s%% of the genes were expressing below bin level %s and were removed.', drop_perc,low_cut_oof))
save(Retina_freq_table_expressingGenes, file='Retina_Freq_Table_ExpressingGenes.RData')


test <-  %>%
  filter(CellType == cell_type) %>%
  filter(!Gene %in% drop_genes)

for (i in 1:100) {
  test <- Retina_freq_table_expressingGenes %>%
    filter(CellType == 'RGC') %>%
    filter(Gene == 'Opn4') %>%
    ggplot(aes(x=ExpressionLevel, y= NormalizedFreq, color = CellClass)) +
    geom_line()
  print(test)
}
