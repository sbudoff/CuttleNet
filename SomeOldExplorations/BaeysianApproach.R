## TITLE: How Credible Are Published Cell Classifications?
## ABSTRACT: The NIH funds many important cutting edge basic research projects that open the doors to further, more sophisticated experiments. A particularly important advance, recently, has been the advent of single cell sequencing databases, in which ever more specific cellular identities are revealed on the basis of the many thousands of genes they can express. As this type of analysis is cost prohibitive, and unnecessary for many more specific projects, particularly those related to specific diseases, it is common for investigators to propose small subsets of these genes to categorize the subtypes they discover. Unfortunately, these recommended subsets are not commonly associated with uncertainty measures, even though the genes that are used have variable expressions. Thus the aim of this project is to assess the uncertainty, using Bayesian modeling, of such a published categorization scheme.

library(tidyverse)
library(foreach)
library(brms)
library(loo)

Workstation = T

if (Workstation) {
  df_root = '/home/sam/Classes/Stats/Bayes/assignments/' # Workstation
  # save_path = '/home/sam/Classes/Stats/Bayes/ALLclusts1hypo_50000samples_multiFit.Rdata' # Workstation
  save_path = '/home/sam/Classes/Stats/Bayes/fit1_b10clusts1hypo_50000_multi.Rdata' # Workstation
  df_path = '/home/sam/Classes/Stats/Bayes/dfPLUShypothesis.Rdata'
  load(df_path)
} else {
  df_root = '/home/sam/Classes/bayesAssignments/' # Laptop
  save_path = '/home/sam/Classes/bayesAssignments/fit1_b5clusts1hypo_50000_multi.Rdata' # laptop
}

# load(save_path)

############################# FUNCTIONS ####################################################
one_hot_matrix <-function(y) {
  n_cols = length(levels(y))
  n_rows = length(y)
  out <- matrix(0, nrow=n_rows, ncol=n_cols)
  for (i in 1:n_rows) {
    out[i,as.numeric(y[i])] <- 1
  }
  return(out)
}

train.test.split <- function(df_sub, frac = 0.80){
  # #Shuffle the data
  # dt = sort(sample(nrow(df_sub), nrow(df_sub)*frac))
  # # extract training set
  # df_train <- df_sub[dt,]
  # y_train <- y[dt,]
  # # extract training set
  # df_test <- df_sub[-dt,] 
  # y_test <- y[-dt,]
  # 
  # df_test <- df_test %>%
  #   select(-Cluster)
  # df_train <- df_train %>%
  #   select(-Cluster)
  # Shuffle the data and extract training set
  dt <- sort(sample(nrow(df_sub), nrow(df_sub)*frac))
  y <- one_hot_matrix(df_sub$Cluster)
  X <- as.matrix(df_sub[,!(names(df_sub) == "Cluster")])
  df_train <- data.frame(X[dt,])
  y_train <- y[dt,]
  
  # Extract test set
  df_test <- data.frame(X[-dt,]) 
  y_test <- y[-dt,]
  
  
  return(list(df_train, df_test, y_test, y_train))
}
################################################################################
clusterHypothesis <- function(clust_ID, df_sub, iter=5000, chains=4) {
  # hypothesis_i <-hypothesis %>%
  #   filter(Cluster == clustID)
  # 
  # genes_hi <- hypothesis_i %>%
  #   select(Gene1,Gene2,Gene3) %>%
  #   gather("id", "Gene") %>%
  #   drop_na() %>%
  #   pull(Gene)
  
  # Subset hypothesis data frame
  hypothesis_i <- subset(hypothesis, Cluster == clust_ID)
  
  # Select Gene1, Gene2, and Gene3 columns and reshape to long format
  genes_hi <- reshape(hypothesis_i[, c("Gene1", "Gene2", "Gene3")], 
                      direction = "long", 
                      varying = list(c("Gene1", "Gene2", "Gene3")),
                      idvar = "id", 
                      times = c("Gene1", "Gene2", "Gene3"), 
                      v.names = "Gene")
  
  # Remove rows with missing values and extract Gene column as a vector
  genes_hi <- na.omit(genes_hi)$Gene
  
  df_hi <- df_sub[, c("Cluster", genes_hi)]
  
  # df_hi <- df_sub %>%
  #   select(Cluster, all_of(genes_hi)) 
  
  tts_list <- train.test.split(df_hi)
  df_train <- tts_list[[1]]
  df_test  <- tts_list[[2]]
  y_test <- tts_list[[3]]
  y_train  <- tts_list[[4]]
  df_train$size = 1
  df_train$Cluster <- with(df_train, cbind(y_train))
  
  fit <- brms::brm(Cluster | trials(size) ~ .,
                   data=df_train,
                   family = brms::multinomial(), #prior = prior,
                   chains = chains, iter = iter, cores = chains)
  
  list(clust_ID, fit, df_test, df_train, y_test, y_train)
}
################################################################################
################################################################################


############## Build From Scratch Bayesian Model ###############################
if (!exists("fits_tranhypos")) {
############## Load Raw Data ###################################################
  df_path = 'keyGeneExpressionBySubtype.RData'
  raw_seq_path='/home/sam/scRNAseq/SCP509/expression/RGC_Atlas.csv'
  raw_clust_path="/home/sam/scRNAseq/SCP509/cluster/RGC_Atlas_coordinates.txt"
  
  if (file.exists(paste0(df_root,df_path))) {
    load(paste0(df_root,df_path))
  } else {
    # This code block reads in the full dataset from tran 2019
    expMatrix <- read.csv(raw_seq_path, head = TRUE) #, sep="\t",skip=1, col.names = c("ID", "UMAP1", "UMAP2", "Cluster"))
    rownames(expMatrix) <- expMatrix$GENE
    
    # Open and find the cluster specific sample IDs
    clusters <- read.csv(raw_clust_path, head = TRUE, sep="\t",skip=1, col.names = c("ID", "UMAP1", "UMAP2", "Cluster", "Experiment","SampleID", "BatchID","Background"))
    RGC_cluster_expression <- select(expMatrix, GENE)
    rownames(RGC_cluster_expression) <- RGC_cluster_expression$GENE
    
    # Manually enter in the Tran proposed classification scheme
    # Proposed marker sets from Tran (2019) Figure 1F
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
    
    # Add in the extra key genes I care about considering
    key_genes = c('Rbpms',# 'Opn1sw', 'Grm6',
                  'Zic1', 'Mafb','Etv1',
                  'Prkcq', 'Penk', 'Coch',
                  'Opn4', 'Meis2', 'Spp1',
                  'Neurod2', 'Kctd4',
                  'Lmo2',  'Cd24a')
    key_genes <- unique(c(unique_genes, key_genes))
  
    # Clean the cluster by ID pairs for merging
    clusters <- clusters %>%
      select(ID, Cluster) %>%
      mutate(ID = gsub("-", ".", ID))
      
    # Subset the data such that I have observations from each cell with genes as regressors and cluster assignment as a column
    df <- expMatrix %>%
      filter(GENE %in% key_genes) %>%
      select(-GENE) %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column('ID') %>%
      left_join(clusters, by='ID') %>%
      select(-ID) %>%
      mutate(Cluster = as.factor(Cluster))
    
    if (sum("4833423E24Rik" == key_genes) == 1) {
      key_genes = sub("4833423E24Rik" , "four833423E24Rik", key_genes)
      df <- df %>%
        mutate(four833423E24Rik = `4833423E24Rik`) %>%
        select(-`4833423E24Rik`)
    }
    
    # Save the list of frequency dataframes for easy loading elsewhere
    save(df, hypothesis, key_genes, file=paste0(df_root,df_path))
  }
############## Filter Raw Data ################# ###############################
# Condition data on RBPMS expression to build realistic models
  df <- df %>%
    filter(Rbpms != 0) 
  
  # Subset the data to only include observations from the 3 most prevelent categories
  obs.counts <- df %>%
    select(Cluster) %>%
    count(Cluster) %>%
    arrange(desc(n))
  
  topN <- obs.counts %>%
    top_n(10) %>%
    pull(Cluster)
  
  bottomN <- obs.counts %>%
    top_n(-10) %>%
    pull(Cluster)
  
  df_sub <- df %>%
    filter(Cluster %in% bottomN) #%>%
    # mutate(four833423E24Rik = `4833423E24Rik`) %>%
    # select(-`4833423E24Rik`)
  
  # Turn Cluster into a matrix
  y <- df_sub %>%
    select(Cluster) %>%
    mutate(ID = row_number(),
           value = 1) %>%
    spread(Cluster, value, fill = 0) %>%
    select(-ID) %>%
    as.matrix()
  
  
  # Clear memory prior to model construction
  rm(df, topN, df_path, df_root, raw_clust_path, raw_seq_path)
  gc()

############## Compute Bayesian Model ##########################################
# Model will be a multinomial logistic regression for the 3 categroical variables (Clusters)

#The multinomial family requires a matrix as input for the response, with one column per category. There are some examples of how to use the multinomal family 
#https://github.com/paul-buerkner/brms/issues/725
#https://discourse.mc-stan.org/t/example-with-family-multinomial/8707

# # In case the loop breaks early
# key_genes <- key_genes[!(key_genes %in% names(fits_1gene))]

# for (i in 1:length(bottomN)){
#   Cluster <- fits_tranhypos[[i]][[1]]
#   y <- fits_tranhypos[[i]][[5]]
#   X <- fits_tranhypos[[i]][[3]]
#   n_genes <- ncol(X) - (1+length(bottomN))
#   posteriors <- fit_sum$fixed$Estimate
# }
  
  # ### Loop through Tran hypothesis
  # n.cores <- parallel::detectCores()/5
  # #create the cluster
  # my.cluster <- parallel::makeCluster(
  #   n.cores, 
  #   type = "PSOCK"
  # )
  # #register it to be used by %dopar%
  # doParallel::registerDoParallel(cl = my.cluster)
  # 
  # fits_tranhypos <- foreach(i = 1:length(unique(df_sub$Cluster))) %dopar% {
  #   clusterHypothesis(unique(df_sub$Cluster)[i], df_sub, iter=100, chains=4)
  # }
  # parallel::stopCluster(cl = my.cluster)
  
  fits_tranhypos <- list()
  for (i in 1:length((unique(df_sub$Cluster)))) {
    fits_tranhypos[i] <- clusterHypothesis(unique(df_sub$Cluster)[i], df_sub, iter=100, chains=4)
  }
  
  save(fits_tranhypos, df_sub, hypothesis, file=save_path)
}

############## Compute Results from test of Models##############################
if (!exists("results_df")) {
  if (Workstation){
    ### Loop through Tran hypothesis
    n.cores <- parallel::detectCores()-5
    #create the cluster
    my.cluster <- parallel::makeCluster(
      n.cores, 
      type = "PSOCK"
    )
    #register it to be used by %dopar%
    doParallel::registerDoParallel(cl = my.cluster)
    IC_list <- foreach(i = 1:length(fits_tranhypos)) %dopar% {
      hyp_clust <- as.character(fits_tranhypos[[i]][[1]])
      WAIC <- waic(fits_tranhypos[[i]][[2]])
      LOOIC <- loo(fits_tranhypos[[i]][[2]])
      list(hyp_clust, WAIC, LOOIC)
    }
    parallel::stopCluster(cl = my.cluster)
    
    WAIC = numeric(length(IC_list))
    LOOIC = numeric(length(IC_list))
    hyp_clust = numeric(length(IC_list))
    for (i in 1:length(IC_list)) {
      hyp_clust[i] <- IC_list[[i]][[1]]
      WAIC[i] <- IC_list[[i]][[2]]$waic
      LOOIC[i] <- IC_list[[i]][[3]]$looic
    }
    IC_df <- tibble(Cluster = hyp_clust, WAIC, LOOIC)
  } else{ 
    IC_list <- list()
    for(i in 1:length(fits_tranhypos)) {
      hyp_clust <- as.character(fits_tranhypos[[i]][[1]])
      WAIC <- waic(fits_tranhypos[[i]][[2]])
      LOOIC <- loo(fits_tranhypos[[i]][[2]])
      IC_list[[i]] <- list(hyp_clust, WAIC, LOOIC)
    }
  }
  save(df_sub, hypothesis, IC_list, IC_df, file=save_path)
  
  for (i in 1:length(fits_tranhypos)) {
    fitted_i <- fitted(object = fits_tranhypos[[i]][[2]])
    fits_tranhypos[[i]][[7]] <- fitted_i
    fits_tranhypos[[i]][[3]]$size <- 1
    predict_i <- predict(object = fits_tranhypos[[i]][[2]],
                         newdata = fits_tranhypos[[i]][[3]])
    fits_tranhypos[[i]][[8]] <- predict_i
  }
  save(fits_tranhypos, df_sub, hypothesis, IC_list, IC_df, file=save_path)
  
  for (i in 1:length(fits_tranhypos)) {
    hyp_clust <- as.character(fits_tranhypos[[i]][[1]])
    hyp_y <- fits_tranhypos[[i]][[5]]
    est <- paste0("Estimate.",hyp_clust)
    test<-data.frame(predict_i) %>%
      select(contains(hyp_clust)) %>%
      mutate(y_hat = (!!sym(est) > 0.5)+0,
             y = hyp_y,
             Correct = (y_hat == y)+0)
    response <- data.frame(y = fits_tranhypos[[i]][[5]][,hyp_clust],
                           pred = test$y_hat) %>%
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
    fits_tranhypos[[i]][[9]] <- response
    fits_tranhypos[[i]][[10]] <- resp_sum
    print(hyp_clust)
    print(resp_sum)
  }
  save(fits_tranhypos, df_sub, hypothesis, IC_list, IC_df, file=save_path)
  
  
  results_df <- fits_tranhypos[[1]][[10]] %>%
    mutate(Cluster = fits_tranhypos[[1]][[1]] )
  for (i in 2:length(fits_tranhypos)) {
    temp <- fits_tranhypos[[i]][[10]] %>%
      mutate(Cluster = fits_tranhypos[[i]][[1]] )
    results_df <- rbind(results_df, temp)
  }
  results_df <- left_join(IC_df, results_df, by = 'Cluster')
  
  save(fits_tranhypos, df_sub, hypothesis, IC_list, IC_df, results_df, file=save_path)
}

############## Visualize Results ###############################################

if (Workstation) {
  hist_root = '/home/sam/scRNAseq/FreqTables/' # Workstation
} else {
  hist_root = '/home/sam/scRNAseq/FreqTables50/' # Laptop
}

hist_path = 'RGC_cluster_ExpressionFrequencies_allGenes.RData'
load(file=paste0(hist_root,hist_path))

Ns <- tibble(Subtype = RGC_cluster_hist[['N ID']],
             N = RGC_cluster_hist[['N']]) %>%
  mutate(p = N/RGC_cluster_hist[['N']][['Full']]) %>%
  filter(Subtype != 'Full')

Ns %>%
  ggplot(aes(x=reorder(Subtype, -p), y = p)) +
  geom_bar(stat='Identity') +
  ggtitle('Proportion of Subtypes') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90))


IC_df %>%
  gather("Metric", 'IC', -Cluster) %>%
  ggplot(aes(x=Cluster, y = IC, color=Metric)) +
  geom_point() +
  ggtitle("Information Criterion for Tran 2019 Proposals") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

clusts <- df_sub$Cluster
df_bubble <- df_sub %>%
  select(-Cluster) %>%
  scale(center=F) %>%
  as.data.frame() %>%
  mutate(Cluster = clusts) %>%
  group_by(Cluster) %>%
  summarise(across(everything(), mean)) %>%
  ungroup() %>%
  gather("Gene", "val", -Cluster)

clusts <- results_df$Cluster
df_bubble %>%
  filter(Cluster %in% clusts) %>%
  ggplot(aes(x=Gene, y=Cluster, size=val, color=val)) +
  geom_point() +
  scale_color_continuous(low = "yellow2", high = "red") +
  ggtitle("5 Rarest Classes Relative to All Proposed Genes")  +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "none") 

genes <- hypothesis %>%
  filter(Cluster %in% clusts) %>%
  select(Gene1, Gene2, Gene3) %>%
  gather("ID", "Gene") %>%
  drop_na() %>%
  pull(Gene) %>%
  unique()

df_bubble %>%
  filter(Gene %in% genes) %>%
  filter(Cluster %in% clusts) %>%
  ggplot(aes(x=Gene, y=Cluster, size=val, color=val)) +
  geom_point() +
  scale_color_continuous(low = "yellow2", high = "red") +
  ggtitle("5 Rarest Classes Relative to Specific Proposed Genes")  +  
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "none") 

results_df %>%
  ggplot(aes(x=Cluster, y = TPR)) +
  geom_col() +
  ggtitle("Ability of Bayesian Model To Correctly Classify Based on Tran Hypothesis") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "none") 


### Loop through 1 gene

# fits_1gene = list()

# n.cores <- parallel::detectCores()/5
# #create the cluster
# my.cluster <- parallel::makeCluster(
#   n.cores, 
#   type = "PSOCK"
# )
# #register it to be used by %dopar%
# doParallel::registerDoParallel(cl = my.cluster)
# 
# 
# clusterGenes <- function(gene, df_sub, dt, iter=5000, chains=4) {
# 
#   df_train <- df_sub[dt,]
#   df_train <- df_train[gene]
#   df_train$size = 1
#   df_train$Cluster <- with(df_train, cbind(y_train))
# 
#   brms::brm(Cluster | trials(size) ~ .,
#       data=df_train,
#       family = brms::multinomial(), #prior = prior,
#       chains = chains, iter = iter, cores = chains)
# }
# 
# 
# fits_1gene <- foreach(i = 1:length(key_genes)) %dopar% {
#   clusterGenes(key_genes[i], df_sub, dt, iter=50000, chains=4)
# }
# parallel::stopCluster(cl = my.cluster)
# save(fits_1gene, file='/home/sam/Classes/Stats/Bayes/fit1_b5clusts1gene_50000_multi.Rdata')

# load(save_path)

# WAIC_1gene = numeric(length(key_genes))
# LOOIC_1gene = numeric(length(key_genes))
# for (i in 1:length(key_genes)){
#   WAIC_1gene[i] <- waic(fits_1gene[[i]])$waic
#   LOOIC_1gene[i] <- loo(fits_1gene[[i]])$looic
# }
# IC_df_1gene <- tibble(gene = key_genes, WAIC, LOOIC)
# 
# IC_df_1gene %>%
#   gather("Metric", 'IC', -gene) %>%
#   ggplot(aes(x=gene, y = IC, color=Metric)) +
#   geom_point() +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90))






#prior <- prior()
# Full moel takes ~2 hours with 2000 samples (Chains are parallel)
# fit <- brm(Cluster | trials(size) ~ .,
#             data=df_train, 
#             family = multinomial(), #prior = prior,
#             chains = 4, iter = 100, cores = 18)
# 
# save(fit, file='/home/sam/Classes/Stats/Bayes/fit1_b5clusts64genes_10000_multi.Rdata')


fit_sum_df <- summary(fit)$fixed 

test <- fit_sum_df %>%
  mutate(believable = (`l-95% CI` >= 0) == (`u-95% CI` >= 0)) %>%
  filter(believable == T) %>%
  rownames_to_column("mu") %>%
  mutate(GENE = sub(".*_", "", mu))

genes <- test$GENE

best_waic = Inf
waics <- numeric(length(names(fits_1gene)))
i <- 1
for (gene in names(fits_1gene)) {
  waic_i <- waic(fits_1gene[[gene]])
  if (waic_i < best_waic) {
    best_waic <- waic_i
    goi <- gene
  }
  waics[i] <- waic_i
  i = i + 1
}

# Use WIC and/or LOOIC to assess model complexity and fit
# Use prediction error for both training and test set to decide on best models
# Set this up such that the regressor considered can easily be changed so we can easily cycle through all possible predictors and build up to the best set
# Use uniform (1/3, 1/3, 1/3) for prior
# Compare models with the starting hypothesis then add


