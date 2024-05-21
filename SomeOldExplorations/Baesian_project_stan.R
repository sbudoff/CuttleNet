library(tidyverse)
library(foreach)
library(rstan)
library(loo)

Workstation = T
B = 5000

if (Workstation) {
  df_root = '/home/sam/Classes/Stats/Bayes/assignments/' # Workstation
  # save_path = '/home/sam/Classes/Stats/Bayes/ALLclusts1hypo_50000samples_multiFit.Rdata' # Workstation
  save_path = '/home/sam/Classes/Stats/Bayes/fitSTAN_b5clusts1hypo_20000.Rdata' # Workstation
  save_path2 = '/home/sam/Classes/Stats/Bayes/fitSTAN_b5clustsNNhypo_20000.Rdata' # Workstation
  df_path = '/home/sam/Classes/Stats/Bayes/dfPLUShypothesis.Rdata'
  load(df_path)
} else {
  df_root = '/home/sam/Classes/bayesAssignments/' # Laptop
  save_path = '/home/sam/Classes/bayesAssignments/fit1_b5clusts1hypo_50000_multi.Rdata' # laptop
}

df_path = 'keyGeneExpressionBySubtype.RData'

load(paste0(df_root,df_path))


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
  top_n(-5) %>%
  pull(Cluster)

df_sub <- df %>%
  filter(Cluster %in% bottomN)

subtype <- as.numeric(as.factor(as.character(df_sub$Cluster)))
X <- as.matrix(df_sub[,!(names(df_sub) == "Cluster")])

# Organize data into list for the stan model
dat <- list()
dat[['k']] <- length(unique(subtype))
dat[['n_obs']] <- dim(X)[1]
dat[['n_genes']] <- dim(X)[2]
dat[['subtype']] <- subtype
dat[['X']] <- X


# Create model. 
mod = "
data {
  int<lower=2> k;  // number of sub types
  int<lower=0> n_obs;  // number of observations
  int<lower=1> n_genes;  // number of regressors
  int<lower=1> subtype[n_obs];     // response variable
  matrix[n_obs, n_genes] X;     // matrix of regressors
}
parameters {
  matrix[n_genes, k] beta;     // matrix of params
}
transformed parameters {
  matrix[n_obs, k] X_beta = X * beta;
}
model {
  to_vector(beta) ~ normal(0,100); // Skeptical prior
  
  for (n in 1:n_obs)
    subtype[n] ~ categorical_logit(X_beta[n]');
}
generated quantities {
  int <lower=1, upper=k> subtype_rep[n_obs];
  vector[n_obs] log_lik;
  for (i in 1:n_obs) {
    subtype_rep[i] = categorical_logit_rng(X_beta[i]');
    log_lik[i] = categorical_logit_lpmf(subtype[i] | X_beta[i]');
  }
}
"
#Compile the model
if (!file.exists("/home/sam/Classes/Stats/Bayes/project_stan_mod.rda")) {
  stan_mod_compiled = stan_model(model_code = mod)
  save(stan_mod_compiled, file = "/home/sam/Classes/Stats/Bayes/project_stan_mod.rda")
}
load("/home/sam/Classes/Stats/Bayes/project_stan_mod.rda")

# Prepare data for model fitting/testing
train.test.split <- function(df_sub, frac = 0.80){
  # #Shuffle the data
  dt <- sort(sample(nrow(df_sub), nrow(df_sub)*frac))
  
  subtype <- as.numeric(as.factor(as.character(df_sub$Cluster)))
  X <- as.matrix(df_sub[,!(names(df_sub) == "Cluster")])
  

  # Extract test set
  df_test <- data.frame(X[-dt,]) 
  y_test <- subtype[-dt]
  
  # Organize data into list for the stan model
  dat <- list()
  dat[['k']] <- length(unique(subtype))
  dat[['n_obs']] <- dim(X)[1]-length(y_test)
  dat[['n_genes']] <- dim(X)[2]
  dat[['subtype']] <- subtype[dt]
  dat[['X']] <- X[dt,]
  
  return(list(dat, df_test, y_test))
}
summaryCleaner <- function(df_genes, fit){
  gene_list <- colnames(select(df_genes, -Cluster)) 
  gene_inds <- which(gene_list %in% NN_mod_3)
  
  gene_list[gene_inds]
  
  sum_df <- data.frame(summary(fit)) %>%
    rownames_to_column('Parameter') %>%
    mutate(Gene = gene_list[as.numeric(stringr::str_extract(Parameter, "\\d+"))], 
           Mean = summary.mean, 
           SD = summary.sd,
           q02.5 = summary.2.5.,
           q97.5 = summary.97.5.,
           n_eff = summary.n_eff,
           R_hat = summary.Rhat) %>%
    select(Gene, Parameter, Mean, SD, q02.5, q97.5, n_eff, R_hat)
  
  param_df <- sum_df %>%
    filter(as.numeric(stringr::str_extract(Parameter, "\\d+")) %in% gene_inds,
           stringr::str_starts(Parameter, "beta"))
  
  rep_df <- sum_df %>%
    filter(stringr::str_starts(Parameter, "subtype_rep"))
  
  return(list(param_df, rep_df))
}

# draw samples from the model with all genes
tts_list <- train.test.split(df_sub)
dat <- tts_list[[1]]
df_test  <- tts_list[[2]]
y_test <- tts_list[[3]]
if (!file.exists("/home/sam/Classes/Stats/Bayes/project_full_stan.rda")) {
  full_fit = sampling(stan_mod_compiled, data = dat,
                      iter = B, chains = 4, cores = 4)
  full_fit_sum <- summary(full_fit)
  ll_full = extract_log_lik(full_fit, merge_chains = FALSE)
  save(full_fit_sum, ll_full, file = "/home/sam/Classes/Stats/Bayes/project_full_stan.rda")
}
load("/home/sam/Classes/Stats/Bayes/project_full_stan.rda")
################################################################################
tranGenes <- function(hypothesis, clust_ID) {
  # Subset hypothesis data frame
  hypothesis_i <- subset(hypothesis, Cluster == clust_ID)
  
  # Select Gene1, Gene2, and Gene3 columns and reshape to long format
  genes_hi <- reshape(hypothesis_i[, c("Gene1", "Gene2", "Gene3")], 
                      direction = "long", 
                      varying = list(c("Gene1", "Gene2", "Gene3")),
                      idvar = "id", 
                      times = c("Gene1", "Gene2", "Gene3"), 
                      v.names = "Gene")
  # pull genes as vector
  genes_hi <- na.omit(genes_hi)$Gene
  
  return(genes_hi)
}

clusterHypothesis <- function(genes_hi, stan_mod_compiled, df_sub, ID, iter=5000, chains=4) {
  # Remove rows with missing values and extract Gene column as a vector
  
  
  # Set all non-relevant genes to zro to drop them from the model
  df_hi <- df_sub
  null_genes <- names(df_sub)
  null_genes <- null_genes[!(null_genes %in% c("Cluster", genes_hi))]
  df_hi[, null_genes] = 0
  
  tts_list <- train.test.split(df_hi)
  dat <- tts_list[[1]]
  df_test  <- tts_list[[2]]
  y_test <- tts_list[[3]]
  
  fit <- rstan::sampling(stan_mod_compiled, data = dat,
                      iter = iter, chains = chains, cores = chains)
  fit_sum <- summary(fit)
  ll <- loo::extract_log_lik(fit, merge_chains = FALSE)
  
  out <- list()
  out[["ID"]] <- ID
  out[["fit"]] <- fit
  out[["summary"]] <- fit_sum
  out[["loglik"]] <- ll
  out[["X_test"]] <- df_test
  out[["y_test"]] <- y_test
  return(out)
}


if (!file.exists(save_path)) {
  fits_tranhypos <- list()
  for (i in 1:length((unique(df_sub$Cluster)))) {
    clust_ID = unique(df_sub$Cluster)[i]
    genes_hi_tran <- tranGenes(hypothesis, clust_ID)
    fits_tranhypos[[i]] <- clusterHypothesis(genes_hi_tran, stan_mod_compiled, df_sub, clust_ID, iter=B, chains=4)
  }
  
  save(fits_tranhypos, df_sub, hypothesis, file=save_path)
}


# Fit neural network proposed models
df_path2 = '/home/sam/scRNAseq/RNAscope/RNAscope/RGC_df354.Rdata'

load(df_path2)
df_350 <- df %>%
  filter(Rbpms != 0) 

# Subset the data to only include observations from the 3 most prevelent categories
obs.counts <- df_350 %>%
  select(Cluster) %>%
  count(Cluster) %>%
  arrange(desc(n))

bottomN <- obs.counts %>%
  top_n(-5) %>%
  pull(Cluster)

df_sub_350 <- df_350 %>%
  filter(Cluster %in% bottomN)

# Models to test discovered via NN
NN_mod_1 = c('Fes', 'Gnb3', 'Igfbp5', 'Isl2', 'Kcnip4', 'Kit', 'Meis2', 'Neurod2', 'Pcdh20',
             'Pcdh7', 'Penk', 'Prokr1', 'Prom1', 'Vit', 'Wscd1', 'Zic1', 'four833423E24Rik')

NN_mod_2 = c('Amigo2', 'Ano3', 'Calca', 'Col12a1', 'Col25a1', 'Fes', 'Gabrg3', 'Gad1', 'Gnb3', 
             'Hcn1', 'Igfbp5', 'Isl2', 'Kcnip4', 'Kcnk2', 'Kit', 'Maf', 'Meis2', 'Mmp17', 'Necab2',
             'Nefh', 'Neurod2', 'Pcdh11x', 'Pcdh20', 'Pcdh7', 'Penk', 'Plpp4', 'Pou3f1', 'Pou6f2', 
             'Prlr', 'Prokr1', 'Prom1', 'Pvalb', 'Qpct', 'Slc7a11', 'Sv2b', 'Syndig1l', 'Syt6', 
             'Tfap2d', 'Trpm1', 'Vit', 'Wscd1', 'Zic1', 'four833423E24Rik')

NN_mod_3 = c('Amigo2', 'Ano3', 'Calca', 'Col12a1', 'Col25a1', 'Fes', 'Gabrg3', 'Gad1', 'Gnb3', 
             'Hcn1', 'Igfbp5', 'Isl2', 'Kcnip4', 'Kcnk2', 'Kit', 'Maf', 'Meis2', 'Mmp17', 'Necab2',
             'Nefh', 'Neurod2', 'Pcdh11x', 'Pcdh20', 'Pcdh7', 'Penk', 'Plpp4', 'Pou3f1', 'Pou6f2', 
             'Prlr', 'Prokr1', 'Prom1', 'Pvalb', 'Qpct', 'Slc7a11', 'Sv2b', 'Syndig1l', 'Syt6', 
             'Tfap2d', 'Trpm1', 'Vit', 'Wscd1', 'Zic1', 'four833423E24Rik', 
             'Bhlhe22', 'Stxbp6', 'Cdhr1', 'Adcyap1')

NN_mod_4 = c('Ano3', 'Apela', 'Bhlhe22', 'Calca', 'Car8', 'Chrm2', 'Col12a1', 'Col25a1',
             'Ebf3', 'Fes', 'Gabrg3', 'Gngt1', 'Hap1', 'Igfbp5', 'Isl1', 'Isl2', 'Kcnip4',
             'Kcnk2', 'Kit', 'Lxn', 'Lypd1', 'Maf', 'Meis2', 'Mmp17', 'Necab2', 'Nefh',
             'Neurod2', 'Nr4a2', 'Pcdh11x', 'Pcdh20', 'Pcdh7', 'Penk', 'Pou3f1', 'Pou6f2',
             'Prlr', 'Prokr1', 'Ptgds', 'RP23-407N2.2', 'Slc17a7', 'Syndig1l', 'Syt6',
             'Tbx20', 'Tfap2d', 'Vit', 'Wscd1', 'Zic1', 'four833423E24Rik')

NN_mod <- c(NN_mod_1, NN_mod_2, NN_mod_3, NN_mod_4)

fits_NNhypos <- list()
for (i in 1:4) {
  fits_NNhypos[[i]] <- clusterHypothesis(NN_mod[i], stan_mod_compiled, df_sub_350, ID = paste0('NN mod ', i), iter=B, chains=4)
}

save(fits_NNhypos, df_sub, hypothesis, file=save_path2)




clean_sum <- summaryCleaner(df_350, fits_NNhypos[[1]][[2]])
param_df <- clean_sum[[1]]
rep_df <- clean_sum[[2]]

# # Set up the parallel backend
# library(doParallel)
# library(foreach)
# n.cores <- parallel::detectCores()/5
# #create the cluster
# my.cluster <- parallel::makeCluster(
#   n.cores,
#   type = "PSOCK"
# )
# #register it to be used by %dopar%
# doParallel::registerDoParallel(cl = my.cluster)
# 
# # Set up a scratch directory for your intermediate files
# intermediate_directory <- 'INTER_DIR'
# if (!dir.exists(intermediate_directory)) {
#   dir.create(intermediate_directory)
# }
# 
# # Run your parallel loop
# foreach(i = 1:length((unique(df_sub$Cluster))), .export = ls(environment())) %dopar% {
#   
#   # Create a unique filename for each interation of the parallel loop
#   each_filename <- paste0('RESULT_', as.character(i), '.rda') 
#   each_filepath <- file.path(intermediate_directory, each_filename)
#   
#   # If the file exists, skip to the next iteration
#   if (file.exists(each_filepath)) {
#     next
#   }
#   
#   # Otherwise, run your code
#   each_result <- clusterHypothesis(unique(df_sub$Cluster)[i], stan_mod_compiled, df_sub, iter=100, chains=4)
#   
#   # Save the result individually
#   save(each_result, file = each_filepath)
#   
#   # OR, save the contents of an environment:
#   save(list = ls(environment()), file = each_filepath)
#   
# }
# parallel::stopCluster(cl = my.cluster)