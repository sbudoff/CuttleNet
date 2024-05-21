## TITLE: How Credible Are Published Cell Classifications?
## ABSTRACT: The NIH funds many important cutting edge basic research projects that open the doors to further, more sophisticated experiments. A particularly important advance, recently, has been the advent of single cell sequencing databases, in which ever more specific cellular identities are revealed on the basis of the many thousands of genes they can express. As this type of analysis is cost prohibitive, and unnecessary for many more specific projects, particularly those related to specific diseases, it is common for investigators to propose small subsets of these genes to categorize the subtypes they discover. Unfortunately, these recommended subsets are not commonly associated with uncertainty measures, even though the genes that are used have variable expressions. Thus the aim of this project is to assess the uncertainty, using Bayesian modeling, of such a published categorization scheme.

library(tidyverse)

df_root = '/home/sam/Classes/Stats/Bayes/assignments/' # Workstation
# df_path = 'keyGeneExpressionBySubtype.RData'
raw_seq_path='/home/sam/scRNAseq/SCP509/expression/RGC_Atlas.csv'
raw_clust_path="/home/sam/scRNAseq/SCP509/cluster/RGC_Atlas_coordinates.txt"
# 
# if (file.exists(paste0(df_root,df_path))) {
#   load(paste0(df_root,df_path))
# } else {
#   # This code block reads in the full dataset from tran 2019
  # expMatrix <- read.csv(raw_seq_path, head = TRUE) #, sep="\t",skip=1, col.names = c("ID", "UMAP1", "UMAP2", "Cluster"))
  # rownames(expMatrix) <- expMatrix$GENE
  # 
  # # Open and find the cluster specific sample IDs
  # clusters <- read.csv(raw_clust_path, head = TRUE, sep="\t",skip=1, col.names = c("ID", "UMAP1", "UMAP2", "Cluster", "Experiment","SampleID", "BatchID","Background"))
  # RGC_cluster_expression <- select(expMatrix, GENE)
  # rownames(RGC_cluster_expression) <- RGC_cluster_expression$GENE
#   
#   # Manually enter in the Tran proposed classification scheme
#   # Proposed marker sets from Tran (2019) Figure 1F
  h = list()
  h[1] <- data.frame(c(Cluster="1_W3D1.1", Gene1='Serpine2', Gene1_pos=T,Gene2='Amigo2',Gene2_pos=T,Gene3=NA,Gene3_pos=NA))
  h[2] <- data.frame(c(Cluster="2_W3D1.2",Gene1='Lypd1',Gene1_pos=T,Gene2='Ntrk1', Gene2_pos=F, Gene3=NA,Gene3_pos=NA))
  h[3] <- data.frame(c(Cluster="3_FminiON",Gene1='Foxp2',Gene1_pos=T,Gene2='Irx4',Gene2_pos=T,Gene3=NA, Gene3_pos=NA))
  h[4] <- data.frame(c(Cluster="4_FminiOFF",Gene1='Pde1a',Gene1_pos=T,Gene2=NA,Gene2_pos=NA,Gene3=NA,Gene3_pos=NA))
  h[5] <- data.frame(c(Cluster="5_J-RGC",Gene1='Tbr1',Gene1_pos=T,Gene2='Pcdh20',Gene2_pos=T, Gene3=NA,Gene3_pos=NA))
  h[6] <- data.frame(c(Cluster="6_W3B",Gene1='Zic1',Gene1_pos=T,Gene2=NA,Gene2_pos=NA, Gene3=NA,Gene3_pos=NA))
  h[7] <- data.frame(c(Cluster="7_Novel",Gene1='Tbx20',Gene1_pos=T,Gene2='Tagln2',Gene2_pos=T,Gene3=NA,Gene3_pos=NA))
  h[8] <- data.frame(c(Cluster="8_Novel",Gene1='Prkcq',Gene1_pos=T,Gene2='Tac1',Gene2_pos=T, Gene3='Spp1',Gene3_pos=F))
  h[9] <- data.frame(c(Cluster="9_Tbr1_Novel",Gene1='Slc7a11',Gene1_pos=T,Gene2='Plpp4',Gene2_pos=T,Gene3=NA,Gene3_pos=NA))
  h[10] <- data.frame(c(Cluster="10_Novel",Gene1='Gpr88',Gene1_pos=T,Gene2=NA,Gene2_pos=NA, Gene3=NA,Gene3_pos=NA))
  h[11] <- data.frame(c(Cluster="11_Novel", Gene1='Serpinb1b',Gene1_pos=T,Gene2='Gm17750',Gene2_pos=T, Gene3=NA,Gene3_pos=NA))
  h[12] <- data.frame(c(Cluster="12_ooDS_NT",Gene1='Mmp17',Gene1_pos=T, Gene2=NA, Gene2_pos=NA, Gene3=NA,Gene3_pos=NA))
  h[13] <- data.frame(c(Cluster="13_Novel", Gene1='Lypd1',Gene1_pos=T,Gene2='Ntrk1',Gene2_pos=T,Gene3=NA,Gene3_pos=NA))
  h[14] <- data.frame(c(Cluster="14_ooDS_Cck",Gene1='Cartpt',Gene1_pos=T,Gene2='Vit',Gene2_pos=T,Gene3=NA,Gene3_pos=NA))
  h[15] <- data.frame(c(Cluster="15_Novel",Gene1='Apela',Gene1_pos=T,Gene2=NA,Gene2_pos=NA,Gene3=NA,Gene3_pos=NA))
  h[16] <- data.frame(c(Cluster="16_ooDS_DV",Gene1='Cartpt', Gene1_pos=T,Gene2='Col25a1',Gene2_pos=T,Gene3=NA,Gene3_pos=NA))
  h[17] <- data.frame(c(Cluster="17_Tbr1_S1",Gene1='Tbr1',Gene1_pos=T,Gene2='Irx4',Gene2_pos=T, Gene3=NA,Gene3_pos=NA))
  h[18] <- data.frame(c(Cluster="18_Novel",Gene1='Pcdh20',Gene1_pos=T,Gene2='4833423E24Rik',Gene2_pos=T,Gene3=NA,Gene3_pos=NA))
  h[19] <- data.frame(c(Cluster="19_Novel",Gene1='Penk',Gene1_pos=T,Gene2='Prdm8', Gene2_pos=T,Gene3='Slc24a2', Gene3_pos=T))
  h[20] <- data.frame(c(Cluster="20_Novel",Gene1='Penk',Gene1_pos=T, Gene2='Gal',Gene2_pos=T, Gene3=NA,Gene3_pos=NA))
  h[21] <- data.frame(c(Cluster="21_Tbr1_S2",Gene1='Tbr1',Gene1_pos=T, Gene2='Calca', Gene2_pos=T,Gene3=NA,Gene3_pos=NA))
  h[22] <- data.frame(c(Cluster="22_M5",Gene1='Serpine2',Gene1_pos=T, Gene2='Cdhr1', Gene2_pos=T, Gene3=NA,Gene3_pos=NA))
  h[23] <- data.frame(c(Cluster="23_W3D2", Gene1='Prokr1',Gene1_pos=T, Gene2=NA,Gene2_pos=NA,Gene3=NA,Gene3_pos=NA))
  h[24] <- data.frame(c(Cluster="24_Novel",Gene1='Fam19a4',Gene1_pos=T,Gene2=NA,Gene2_pos=NA,Gene3=NA, Gene3_pos=NA))
  h[25] <- data.frame(c(Cluster="25_Novel",Gene1='Slc17a7',Gene1_pos=T, Gene2=NA, Gene2_pos=NA,Gene3=NA,Gene3_pos=NA))
  h[26] <- data.frame(c(Cluster="26_Novel",Gene1='Penk',Gene1_pos=T,Gene2='Igfbp5',Gene2_pos=T, Gene3=NA,Gene3_pos=NA))
  h[27] <- data.frame(c(Cluster="27_Novel",Gene1='Prkcq', Gene1_pos=T,Gene2=NA,Gene2_pos=NA,Gene3=NA,Gene3_pos=NA))
  h[28] <- data.frame(c(Cluster="28_FmidiOFF",Gene1='Foxp2',Gene1_pos=T, Gene2='Cdk15',Gene2_pos=T,Gene3=NA,Gene3_pos=NA))
  h[29] <- data.frame(c(Cluster="29_Novel",Gene1='Stxbp6', Gene1_pos=T, Gene2='Prlr', Gene2_pos=T,   Gene3=NA, Gene3_pos=NA))
  h[30] <- data.frame(c(Cluster="30_Novel",Gene1='Postn',Gene1_pos=T,Gene2=NA, Gene2_pos=NA,Gene3=NA,Gene3_pos=NA))
  h[31] <- data.frame(c(Cluster="31_M2",Gene1='Tbx20', Gene1_pos=T,Gene2='Spp1',Gene2_pos=T,Gene3=NA,Gene3_pos=NA))
  h[32] <- data.frame(c(Cluster="32_F_Novel",Gene1='Rhox5',Gene1_pos=T, Gene2=NA, Gene2_pos=NA, Gene3=NA,Gene3_pos=NA))
  h[33] <- data.frame(c(Cluster="33_M1",Gene1='Adcyap1', Gene1_pos=T,Gene2='Opn4',Gene2_pos=T, Gene3='Nmb',  Gene3_pos=F))
  h[34] <- data.frame(c(Cluster="34_Novel",Gene1='Tpbg',  Gene1_pos=T,Gene2='Spp1', Gene2_pos=F,Gene3=NA, Gene3_pos=NA))
  h[35] <- data.frame(c(Cluster="35_Novel",Gene1='Igfbp4',Gene1_pos=T, Gene2='Chrm2', Gene2_pos=T,Gene3=NA, Gene3_pos=NA))
  h[36] <- data.frame(c(Cluster="36_Novel",Gene1='Stxbp6',Gene1_pos=T, Gene2='Coch', Gene2_pos=T, Gene3=NA,Gene3_pos=NA))
  h[37] <- data.frame(c(Cluster="37_Novel",Gene1='Ceacam10', Gene1_pos=T, Gene2=NA,Gene2_pos=NA, Gene3=NA,Gene3_pos=NA))
  h[38] <- data.frame(c(Cluster="38_FmidiON", Gene1='Foxp2',Gene1_pos=T,Gene2='Anxa3', Gene2_pos=T,Gene3=NA, Gene3_pos=NA))
  h[39] <- data.frame(c(Cluster="39_Novel", Gene1='Neurod2',Gene1_pos=T, Gene2='S100b',Gene2_pos=T,Gene3=NA,Gene3_pos=NA))
  h[40] <- data.frame(c(Cluster="40_M1dup",Gene1='Nmb',Gene1_pos=T, Gene2=NA,Gene2_pos=NA, Gene3=NA, Gene3_pos=NA))
  h[41] <- data.frame(c(Cluster="41_AlphaONT", Gene1='Spp1',Gene1_pos=T, Gene2='Kit',Gene2_pos=T, Gene3=NA, Gene3_pos=NA))
  h[42] <- data.frame(c(Cluster="42_AlphaOFFS",Gene1='Spp1',Gene1_pos=T,Gene2='Fes',Gene2_pos=T,Gene3=NA,Gene3_pos=NA))
  h[43] <- data.frame(c(Cluster="43_AlphaONS", Gene1='Spp1', Gene1_pos=T, Gene2='Il1rapl2', Gene2_pos=T,Gene3=NA,Gene3_pos=NA))
  h[44] <- data.frame(c(Cluster="44_Novel",Gene1='Bhlhe22', Gene1_pos=T,Gene2='Fxyd6',Gene2_pos=T,Gene3=NA, Gene3_pos=NA))
  h[45] <- data.frame(c(Cluster="45_AlphaOFFT",Gene1='Spp1', Gene1_pos=T,Gene2='Tpbg',Gene2_pos=T,Gene3=NA, Gene3_pos=NA))

  # Reorganize the manually entered genes into a dataframe
  hypothesis <- t(as.data.frame(h))
  rownames(hypothesis) <- NULL
  hypothesis <- as.data.frame(hypothesis)
  colnames(hypothesis) <- c("Cluster", 'Gene1', 'Gene1_pos', 'Gene2', 'Gene2_pos', 'Gene3', 'Gene3_pos')

  # Extrac the unique genes fro subsetting expMatrix
unique_genes <- hypothesis %>%
  select(Gene1, Gene2, Gene3) %>%
  gather(key, value) %>%
  drop_na() %>%
  distinct(value)%>%
  pull(value)
#   
#   # Add in the extra key genes I care about considering
#   key_genes = c('Rbpms',# 'Opn1sw', 'Grm6',
#                 'Zic1', 'Mafb','Etv1',
#                 'Prkcq', 'Penk', 'Coch',
#                 'Opn4', 'Meis2', 'Spp1',
#                 'Neurod2', 'Kctd4',
#                 'Lmo2',  'Cd24a')
#   key_genes <- unique(c(unique_genes, key_genes))
#   
#   # Clean the cluster by ID pairs for merging
#   clusters <- clusters %>%
#     select(ID, Cluster) %>%
#     mutate(ID = gsub("-", ".", ID))
#   
#   # Subset the data such that I have observations from each cell with genes as regressors and cluster assignment as a column
#   df <- expMatrix %>%
#     filter(GENE %in% key_genes) %>%
#     select(-GENE) %>%
#     t() %>%
#     as.data.frame() %>%
#     rownames_to_column('ID') %>%
#     left_join(clusters, by='ID') %>%
#     select(-ID) %>%
#     mutate(Cluster = as.factor(Cluster))
#   
#   if (sum("4833423E24Rik" == key_genes) == 1) {
#     key_genes = sub("4833423E24Rik" , "four833423E24Rik", key_genes)
#     df <- df %>%
#       mutate(four833423E24Rik = `4833423E24Rik`) %>%
#       select(-`4833423E24Rik`)
#   }
#   
#   # Save the list of frequency dataframes for easy loading elsewhere
#   save(df, hypothesis, key_genes, file=paste0(df_root,df_path))
# }

library(ROCR)
# load('Downloads/keyGeneExpressionBySubtype.RData') #laptop
load('/home/sam/Classes/Stats/Bayes/assignments/keyGeneExpressionBySubtype.RData') 


# Condition data on RBPMS expression to build realistic models
df <- df %>%
  filter(Rbpms != 0) 

for (i in 1:nrow(hypothesis)) {
  h_i <- hypothesis[i,]
  
  df_i <- df %>%
    mutate(y = 0+(Cluster == h_i$Cluster)) %>%
    select(-Cluster)
  
  n_80 = round(nrow(df_i)*0.7)
  train <- na.omit(df_i[1:n_80,])
  test <- na.omit(df_i[n_80+1:nrow(df_i),])
  
  # make a model
  model_i <- glm(y ~.,family=binomial(link='logit'),data=test)
  
  # Get results
  fitted.results <- predict(model_i,test,type='response')
  fitted.results <- ifelse(fitted.results > 0.5,1,0)
  misClasificError <- mean(fitted.results != test$y)
  print(paste0('Accuracy for ', h_i$Cluster, ' using logistic regression is: ', 1-misClasificError))
  

  pr <- prediction(fitted.results, test$y)
  # prf <- performance(pr, measure = "tpr", x.measure = "fpr")
  # plot(prf)
  auc <- performance(pr, measure = "auc")
  auc <- auc@y.values[[1]]
  print(paste0('AUC for ', h_i$Cluster, ' is: ', auc))
  response <- tibble(pred = fitted.results, y=test$y) %>%
    mutate(TP = (y==1 & pred==1),
           TN = (y==0 & pred==0),
           FP = (y==0 & pred==1),
           FN = (y==1 & pred==0))
  resp_sum <- response %>%
    summarize(cells = sum(y),
              TP = sum(TP),
              TN = sum(TN),
              FP = sum(FP),
              FN = sum(FN)) %>%
    mutate(TPR = TP / (TP+FN),
           TNR = TN / (TN+FP),
           Prec = TP/(TP+FP),
           Accurancy = (TP+TN)/(TP+TN+FP+FN))
  print(resp_sum)
}

for (i in 1:nrow(hypothesis)) {
  h_i <- hypothesis[i,]
  h_genes <- h_i %>%
    select(Gene1, Gene2, Gene3) %>%
    slice(1) %>% 
    unlist(., use.names=FALSE)
  
  h_genes <- na.exclude(h_genes)
  
  df_i <- df %>%
    mutate(y = 0+(Cluster == h_i$Cluster)) %>%
    select(y, all_of(h_genes))
  
  n_80 = round(nrow(df_i)*0.7)
  train <- na.omit(df_i[1:n_80,])
  test <- na.omit(df_i[n_80+1:nrow(df_i),])
  
  # make a model
  model_i <- glm(y ~.,family=binomial(link='logit'),data=test)
  
  # Get results
  fitted.results <- predict(model_i,test,type='response')
  fitted.results <- ifelse(fitted.results > 0.5,1,0)
  misClasificError <- mean(fitted.results != test$y)
  print(paste0('Accuracy for ', h_i$Cluster, ' using logistic regression is: ', 1-misClasificError))
  
  
  pr <- prediction(fitted.results, test$y)
  # prf <- performance(pr, measure = "tpr", x.measure = "fpr")
  # plot(prf)
  auc <- performance(pr, measure = "auc")
  auc <- auc@y.values[[1]]
  print(paste0('AUC for ', h_i$Cluster, ' is: ', auc))
  response <- tibble(pred = fitted.results, y=test$y) %>%
    mutate(TP = (y==1 & pred==1),
           TN = (y==0 & pred==0),
           FP = (y==0 & pred==1),
           FN = (y==1 & pred==0))
  resp_sum <- response %>%
    summarize(cells = sum(y),
              TP = sum(TP),
              TN = sum(TN),
              FP = sum(FP),
              FN = sum(FN)) %>%
    mutate(TPR = TP / (TP+FN),
           TNR = TN / (TN+FP),
           Prec = TP/(TP+FP),
           Accurancy = (TP+TN)/(TP+TN+FP+FN))
  print(resp_sum)
}
