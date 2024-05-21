library(tidyverse)

# Load RGC
setwd("/home/sam/scRNAseq/FreqTables/")
data_path = "RGC_cluster_ExpressionFrequencies_allGenes.RData"
load(data_path)
RGC_freq_table = RGC_cluster_hist
remove(RGC_cluster_hist)


# data_path = "RGC_P56__cluster_ExpressionFrequencies_allGenes.RData"
# load(data_path)



# RGC_freq_table2 = cluster_hist
# 
# RGC_freq_table = RGC_freq_table + RGC_freq_table2
# del(RGC_freq_table2)

# Load ACs
data_path = "AC_cluster_ExpressionFrequencies_allGenes.RData"
load(data_path)

AC_freq_table = cluster_hist

# Load BPs
data_path = "BP_cluster_ExpressionFrequencies_allGenes.RData"
load(data_path)

BP_freq_table = cluster_hist
remove(cluster_hist)

# Combine All records into single DF

FindMax <- function(list) {
  ids <- labels(list)
  ids <- ids[!ids %in% c('N', 'N ID', 'Full')]
  max_vec <- c()
  for (id in ids) 
  {
    max_vec <- c(max_vec, max(list[[id]]['value']))
  }
  max_vec
}

PullRecords <- function(list, cellType = 'RGC', skip = c('N', 'N ID', 'Full') ) {
  ids <- labels(list)
  ids <- ids[!ids %in% skip]
  out <- list()
  for (clust in ids) 
  {
    out[[clust]] <- list[[clust]] %>%
      mutate(ExpressionLevel = value) %>%
      select(-value) %>%
      gather(Gene, Frequency, -ExpressionLevel) %>%
      mutate(CellType = cellType,
             CellClass = clust,
             N = list[['N']][clust],
             NormalizedFreq = Frequency/N) %>%
      arrange(Gene, ExpressionLevel) %>%
      select(CellType, CellClass, Gene, ExpressionLevel, Frequency, N, NormalizedFreq)
    print(sprintf('%s gathered and stored', clust))
  }
  out <- data.table::rbindlist(out)
  print(sprintf('%s All merged', cellType))
  out
}

setwd("/home/sam/scRNAseq/FreqTables/")
Retina_freq_table <- list()
Retina_freq_table[['BP']] <- PullRecords(BP_freq_table, cellType = 'BP',
                                         skip = c('N', 'N ID', 'Full', 
                                                  "AC (Amacrine cell)",
                                                  "Doublets/Contaminants",
                                                  "Cone Photoreceptors" ,
                                                  "Rod Photoreceptors",
                                                  "MG (Mueller Glia)") )
Retina_freq_table[['Glia']] <- PullRecords(BP_freq_table, cellType = 'Glia',
                                         skip = labels(BP_freq_table)[!labels(BP_freq_table) %in% "MG (Mueller Glia)" ] )
Retina_freq_table[['Photoreceptors']] <- PullRecords(BP_freq_table, cellType = 'Photoreceptors',
                                         skip = labels(BP_freq_table)[!labels(BP_freq_table) %in% c("Cone Photoreceptors" ,
                                                                                                    "Rod Photoreceptors") ] )
save(Retina_freq_table, file='Full_Retina_FreqTable50.RData')

Retina_freq_table[['AC']] <- PullRecords(AC_freq_table, cellType = 'AC')
save(Retina_freq_table, file='Full_Retina_FreqTable50.RData')

Retina_freq_table[['RGC']] <- PullRecords(RGC_freq_table)
save(Retina_freq_table, file='Full_Retina_FreqTable50.RData')


setwd('/home/sam/scRNAseq/FreqTables/')
load(file='Full_Retina_FreqTable50.RData')
# Save the gene frequency tables in an excel sheet where each sheet is a cluster
# library(openxlsx)
# 
# # create workbook
# wb <- createWorkbook()
# 
# #Iterate through each list element to create a sheet based on that list (creating an anonymous function inside Map())
# Map(function(data, nameofsheet){
#   addWorksheet(wb, nameofsheet)
#   writeData(wb, nameofsheet, data)
# }, Retina_freq_table, names(Retina_freq_table))
# 
# ## Save workbook to excel file
# saveWorkbook(wb, file = "Full_Retina_FreqTable50_allGenes.xlsx", overwrite = TRUE)
library(tidyverse)
Retina_freq_table <- data.table::rbindlist(Retina_freq_table) %>%
  arrange(CellType, CellClass, Gene, ExpressionLevel)
save(Retina_freq_table, file='Full_Retina_FreqTable50.RData')