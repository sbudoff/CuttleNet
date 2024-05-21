library(tidyverse)

save_path <- '/home/sam/scRNAseq/Xenium/NeurIPS/ML_df.RData'

# Check if the file exists
if (file.exists(save_path)) {
  print("The file exists.")
} else {
  load('/home/sam/scRNAseq/FreqTables/ExpressingExpressionMats.RData')

  # Load RGCs
  expMatrix_RGC <- read.csv("/home/sam/scRNAseq/SCP509/expression/RGC_Atlas.csv", head = TRUE)
  rownames(expMatrix_RGC) <- expMatrix_RGC$GENE
  expMatrix_RGC <- expMatrix_RGC[useful_genes,]
  expMatrix_RGC <- na.omit(expMatrix_RGC)
  
  geneNames <- expMatrix_RGC$GENE  # Store gene names before removal
  expMatrix_RGC$GENE <- NULL  # Remove the gene column from dataframe
  
  expMatrix_RGC <- t(expMatrix_RGC)
  colnames(expMatrix_RGC) <- geneNames
  
  clusters <- read.csv("/home/sam/scRNAseq/SCP509/cluster/RGC_Atlas_coordinates.txt", head = TRUE, sep="\t",skip=1, col.names = c("ID", "UMAP1", "UMAP2", "Cluster", "Experiment","SampleID", "BatchID","Background"))
  clusters$ID <- gsub("-", ".", clusters$ID)
  
  RGC_df <- as.data.frame(expMatrix_RGC)
  rownames(RGC_df) <- rownames(expMatrix_RGC)
  RGC_df$Label <- clusters$Cluster[match(rownames(RGC_df), clusters$ID)]
  RGC_df$Label <- as.factor(RGC_df$Label)
  
  
  # Load BCs
  expMatrix_BC <- read.csv("/home/sam/scRNAseq/SCP3/expression/exp_matrix.txt", head = TRUE, sep="\t")
  rownames(expMatrix_BC) <- expMatrix_BC$GENE
  expMatrix_BC <- expMatrix_BC[useful_genes,]
  expMatrix_BC <- na.omit(expMatrix_BC)
  
  geneNames <- expMatrix_BC$GENE  # Store gene names before removal
  expMatrix_BC$GENE <- NULL  # Remove the gene column from dataframe
  
  expMatrix_BC <- t(expMatrix_BC)
  colnames(expMatrix_BC) <- geneNames
  
  clusters <- read.csv("/home/sam/scRNAseq/SCP3/metadata/clust_retinal_bipolar.txt", head = TRUE, sep="\t",skip=1, col.names = c("ID", "Cluster", "subcluster"))
  clusters$ID <- gsub("-", ".", clusters$ID)
  
  BC_df <- as.data.frame(expMatrix_BC)
  rownames(BC_df) <- rownames(expMatrix_BC)
  BC_df$Label <- clusters$Cluster[match(rownames(BC_df), clusters$ID)]
  BC_df$Label <- as.factor(BC_df$Label)
  
  # Load ACs
  expMatrix_AC <- read.csv("/home/sam/scRNAseq/MouseAC_gene_expression_matrix.csv", head = TRUE, sep="\t")
  rownames(expMatrix_AC) <- expMatrix_AC$GENE
  expMatrix_AC <- expMatrix_AC[useful_genes,]
  expMatrix_AC <- na.omit(expMatrix_AC)
  
  geneNames <- expMatrix_AC$GENE  # Store gene names before removal
  expMatrix_AC$GENE <- NULL  # Remove the gene column from dataframe
  
  expMatrix_AC <- t(expMatrix_AC)
  colnames(expMatrix_AC) <- geneNames
  
  clusters <- read.csv("/home/sam/scRNAseq/MouseAC_clusterfile.txt", head = TRUE, sep="\t",skip=1, col.names = c("ID", "UMAP1", "UMAP2", "Cluster"))
  clusters$ID <- gsub("-", ".", clusters$ID)
  
  AC_df <- as.data.frame(expMatrix_AC)
  rownames(AC_df) <- rownames(expMatrix_AC)
  AC_df$Label <- clusters$Cluster[match(rownames(AC_df), clusters$ID)]
  AC_df$Label <- as.factor(AC_df$Label)
  
  # Function to add missing columns from useful_genes to a dataframe and fill with zeros
  add_missing_columns <- function(df, columns) {
    missing_cols <- setdiff(columns, names(df))
    df[missing_cols] <- 0  # Assign zero to all missing columns
    df <- df[, columns]  # Ensure consistent column order
    return(df)
  }
  
  # Apply this function to each dataframe
  useful_genes <- c("Label", useful_genes)
  RGC_df <- add_missing_columns(RGC_df, useful_genes)
  AC_df <- add_missing_columns(AC_df, useful_genes)
  BC_df <- add_missing_columns(BC_df, useful_genes)
  
  # Combine the dataframes
  combined_df <- bind_rows(RGC_df, AC_df, BC_df)
  # Save the combined_df as an RData file
  save(combined_df, file = save_path)
  # Remove all objects from the global environment
  rm(list = ls())
  
}

seed <- 36
data_path <- '/home/sam/scRNAseq/Xenium/NeurIPS/ML_df.RData'
save_path <- paste0('/home/sam/scRNAseq/Xenium/NeurIPS/ML_results_',seed,'.RData')
load(data_path)
load(save_path)


library(caret)
library(rpart)
library(rpart.plot)
set.seed(seed)
shuffledDF <- combined_df[sample(nrow(combined_df)), ]
# shuffledDF <- shuffledDF[1:500,]
shuffledDF <- shuffledDF %>%
  mutate(Label = factor(Label))


split <- createDataPartition(shuffledDF$Label, p = 0.8, list = FALSE)
trainingSet <- shuffledDF[split, ]
testingSet <- shuffledDF[-split, ]
table(trainingSet$Label) 

# Function to modify gene names to ensure compatibility
normalize_gene_names <- function(names) {
  # Replace hyphens with underscores and prepend 'gene' if it starts with a numeric
  names <- gsub("-", "_", names)  # Replace hyphen with underscore
  names <- ifelse(grepl("^[0-9]", names), paste("gene", names, sep="_"), names)
  return(names)
}

# Apply this function to column names of your datasets
colnames(trainingSet) <- normalize_gene_names(colnames(trainingSet))
colnames(testingSet) <- normalize_gene_names(colnames(testingSet))


################################################################################
############### Decision Tree ##################################################
################################################################################

model_DT <- rpart(Label ~ ., data = trainingSet, method = "class")
predictions_DT <- predict(model_DT, testingSet, type = "class")
cMat_DT <- confusionMatrix(predictions_DT, testingSet$Label)

importance <- as.data.frame(model_DT$variable.importance)
colnames(importance) <- "importance"
importance$Gene <- rownames(importance)
rownames(importance) <- NULL

# Assuming the correct column is named correctly
if (nrow(importance) > 300) {
  top_genes_DT <- importance$Gene[1:300]

  trainingSet_reduced <- trainingSet[, c(top_genes_DT, "Label")]
  testingSet_reduced <- testingSet[, c(top_genes_DT, "Label")]

  model_DT <- rpart(Label ~ ., data = trainingSet_reduced, method = "class")
  predictions_DT <- predict(model_DT, testingSet_reduced, type = "class")
  cMat_DT <- confusionMatrix(predictions_reduced, testingSet_reduced$Label)
} else {
  print("Less than 300 genes needed by Decision Tree")
}
# Save the current results
save(model_DT, cMat_DT, file = save_path)

# ################################################################################
# ############### Random Forest ##################################################
# ################################################################################
if (!require("randomForest")) install.packages("randomForest", dependencies = TRUE)
library(randomForest)

# Train the initial Random Forest model
model_rf <- randomForest(Label ~ ., data = trainingSet, ntree = 500, importance = TRUE)

# Extract and sort the feature importance
importance_rf <- importance(model_rf)
feature_importance_rf <- data.frame(Feature = rownames(importance_rf), Importance = importance_rf[, "MeanDecreaseGini"])
sorted_features_rf <- feature_importance_rf[order(-feature_importance_rf$Importance), ]

predictions_rf <- predict(model_rf, testingSet)
cMat_rf <- confusionMatrix(predictions_rf, testingSet$Label)

# Check if there are enough features for subsetting
if (nrow(sorted_features_rf) >= 300) {
  top_300_features_rf <- head(sorted_features_rf, 300)
  top_features_names_rf <- top_300_features_rf$Feature

  # Subset training and testing sets with top 300 features
  trainingSet_rf_reduced <- trainingSet[, c(top_features_names_rf, "Label")]
  testingSet_rf_reduced <- testingSet[, c(top_features_names_rf, "Label")]

  # Retrain Random Forest on the reduced feature set
  model_rf <- randomForest(Label ~ ., data = trainingSet_rf_reduced, ntree = 300)
  predictions_rf <- predict(model_rf, testingSet_rf_reduced)
  cMat_rf <- confusionMatrix(predictions_rf, testingSet_rf_reduced$Label)
  } else {
  print("Less than 300 genes were identified as important by Random Forest.")
}
# Save the current results
save(model_DT, model_rf, cMat_rf, cMat_DT, file = save_path)

################################################################################
############### SVM ######### ##################################################
################################################################################
if (!require("caret")) install.packages("caret", dependencies = TRUE)
library(caret)
if (!require("e1071")) install.packages("e1071", dependencies = TRUE)
library(e1071)

svm_model <- svm(Label ~ ., data = trainingSet[, -ncol(trainingSet)], type = "C-classification", kernel = "linear")
summary(svm_model)

sv_indices <- svm_model$index
coefs <- abs(svm_model$coefs)  # Get absolute values of coefficients
weights <- svm_model$coefs  # Original weights

# Aggregate weights by feature across all one-vs-one classifiers
# If the structure of svm_model does not directly give you this, you might need
# to manually handle the extraction based on your specific SVM training
coefs_abs_sum <- apply(abs(weights), 2, sum)  # Sum across rows for each feature

# Assuming coefs_abs_sum is correctly calculated
if (length(coefs_abs_sum) > 0 && all(names(coefs_abs_sum) %in% colnames(trainingSet))) {
  # Get top feature names, ensuring 'Label' is retained
  top_features_svm <- colnames(trainingSet[, -ncol(trainingSet)])[order(-coefs_abs_sum)][1:length(coefs_abs_sum)]

  if (length(top_features_svm) > 300) {
    top_features_svm <- top_features_svm[1:300]
  }

  # Verify the 'Label' column is included in the subsetting process
  if (!"Label" %in% top_features_svm) {
    top_features_svm <- c(top_features_svm, "Label")
  }

  # Subset the training and testing datasets to include these top features along with 'Label'
  trainingSet_svm_reduced <- trainingSet[, top_features_svm, drop = FALSE]
  testingSet_svm_reduced <- testingSet[, top_features_svm, drop = FALSE]

  # Re-run the SVM model
  svm_model <- svm(Label ~ ., data = trainingSet_svm_reduced, type = "C-classification", kernel = "linear")
  predictions_svm <- predict(svm_model, testingSet_svm_reduced)
  cMat_svm <- confusionMatrix(predictions_svm, testingSet_svm_reduced$Label)

  print(cMat_svm)
} else {
  print("No significant features were identified by the SVM model or feature names mismatch.")
}
# Save the current results
save(model_DT, model_rf, svm_model, cMat_rf, cMat_DT, cMat_svm, file = save_path)

################################################################################
############### GBM ############################################################
################################################################################
if (!require("gbm")) install.packages("gbm", dependencies = TRUE)
if (!require("caret")) install.packages("caret", dependencies = TRUE)
library(gbm)
library(caret)

set.seed(seed)
# Train the initial GBM model
model_gbm <- gbm(Label ~ ., data = trainingSet, distribution = "multinomial", n.trees = 300, interaction.depth = 3, shrinkage = 0.01, n.minobsinnode = 10, verbose = TRUE)

# Extract and sort the feature importance
importance_gbm <- summary(model_gbm, n.trees = 300)  # Ensure n.trees is correctly referred
top_300_features_gbm <- head(importance_gbm, 300)

# Check if there are enough features for subsetting
top_features_names_gbm <- rownames(top_300_features_gbm)
top_features_names_gbm <- normalize_gene_names(top_features_names_gbm)  # Normalize feature names

# Subset training and testing sets with top 300 features
trainingSet_gbm_reduced <- trainingSet[, c(top_features_names_gbm, "Label")]
testingSet_gbm_reduced <- testingSet[, c(top_features_names_gbm, "Label")]

# Retrain GBM on the reduced feature set
model_gbm <- gbm(Label ~ ., data = trainingSet_gbm_reduced, distribution = "multinomial", n.trees = 300, interaction.depth = 3, shrinkage = 0.01, n.minobsinnode = 10, verbose = TRUE)
predictions_gbm_reduced <- predict(model_gbm, testingSet_gbm_reduced, n.trees = 500, type = "response")
predictions_matrix <- predictions_gbm_reduced[, , 1]

# Get the index of the maximum value in each row
max_indices <- max.col(predictions_matrix, ties.method = "first")
# Map these indices to the column names of the prediction matrix
predicted_class_labels <- colnames(predictions_gbm_reduced)[max_indices]

# Convert to factor with the same levels as the original Label column
predicted_class_labels <- factor(predicted_class_labels, levels = levels(testingSet_gbm_reduced$Label))

# Convert predictions to factors and calculate the confusion matrix
# Compute confusion matrix
cMat_gbm <- confusionMatrix(predicted_class_labels, testingSet_gbm_reduced$Label)

# Save the current results
save(model_DT, model_rf, svm_model, model_gbm, cMat_rf, cMat_DT, cMat_svm,cMat_gbm, file = save_path)

#################################################################################################

# Function to extract TPR for each label from a confusion matrix
extract_tpr <- function(cm) {
  tpr <- cm$byClass[, "Sensitivity"]
  return(tpr)
}

# Calculate TPR for each model
tpr_dt <- extract_tpr(cMat_DT)
tpr_rf <- extract_tpr(cMat_rf)
tpr_svm <- extract_tpr(cMat_svm)
tpr_gbm <- extract_tpr(cMat_gbm)

# Create a dataframe for each
df_dt <- data.frame(Label=names(tpr_dt), DecisionTree=tpr_dt)
df_rf <- data.frame(Label=names(tpr_rf), RandomForest=tpr_rf)
df_svm <- data.frame(Label=names(tpr_svm), SVM=tpr_svm)
df_gbm <- data.frame(Label=names(tpr_gbm), GBM=tpr_gbm)

# Merge all dataframes into one
df_tpr <- Reduce(function(x, y) merge(x, y, by="Label", all=TRUE), list(df_dt, df_rf, df_svm, df_gbm))

# Replace NA with 0 if any, because NA can appear if a label was not present in the test set of a particular model
df_tpr[is.na(df_tpr)] <- 0

# View the final dataframe
print(df_tpr)


# Function to extract Precision for each label from a confusion matrix
extract_precision <- function(cm) {
  precision <- cm$byClass[, "Pos Pred Value"]
  return(precision)
}

# Calculate Precision for each model
precision_dt <- extract_precision(cMat_DT)
precision_rf <- extract_precision(cMat_rf)
precision_svm <- extract_precision(cMat_svm)
precision_gbm <- extract_precision(cMat_gbm)

# Create a dataframe for each
df_dt <- data.frame(Label=names(precision_dt), DecisionTree=precision_dt)
df_rf <- data.frame(Label=names(precision_rf), RandomForest=precision_rf)
df_svm <- data.frame(Label=names(precision_svm), SVM=precision_svm)
df_gbm <- data.frame(Label=names(precision_gbm), GBM=precision_gbm)

# Merge all dataframes into one
df_precision <- Reduce(function(x, y) merge(x, y, by="Label", all=TRUE), list(df_dt, df_rf, df_svm, df_gbm))

# Replace NA with 0 if any, because NA can appear if a label was not present in the test set of a particular model
df_precision[is.na(df_precision)] <- 0

# View the final dataframe
print(df_precision)


# Define the threshold
threshold <- 0.9

# Assuming df_tpr and df_precision are already defined and populated as described in previous steps

# Match the labels and methods in both dataframes just to ensure they align correctly
df_combined <- merge(df_tpr, df_precision, by = "Label", suffixes = c("_tpr", "_prec"))

# Function to check if both values exceed the threshold
check_threshold <- function(tpr, prec) {
  return(tpr > threshold & prec > threshold)
}

# Apply the function to each pair of TPR and Precision columns for each method
df_combined$DecisionTree <- mapply(check_threshold, df_combined$DecisionTree_tpr, df_combined$DecisionTree_prec)
df_combined$RandomForest <- mapply(check_threshold, df_combined$RandomForest_tpr, df_combined$RandomForest_prec)
df_combined$SVM <- mapply(check_threshold, df_combined$SVM_tpr, df_combined$SVM_prec)
df_combined$GBM <- mapply(check_threshold, df_combined$GBM_tpr, df_combined$GBM_prec)

# Select only the columns with the final boolean results
df_final <- df_combined[, c("Label", "DecisionTree", "RandomForest", "SVM", "GBM")]

# View the final dataframe
print(df_final)

# Assuming df_final is already defined and populated as previously described

# Sum the boolean values to count labels above threshold for each ML model
summary_counts <- data.frame(
  DecisionTree = sum(df_final$DecisionTree),
  RandomForest = sum(df_final$RandomForest),
  SVM = sum(df_final$SVM),
  GBM = sum(df_final$GBM)
)

# Transpose the dataframe to match the requested format: models as rows, count as a column
summary_counts <- t(summary_counts)
colnames(summary_counts) <- c("Count")

# Convert to a dataframe if needed (since transposing results in a matrix)
summary_counts <- as.data.frame(summary_counts)

# View the summary
print(summary_counts)

# Save the current results
save(summary_counts, df_final, model_DT, model_rf, svm_model, model_gbm, cMat_rf, cMat_DT, cMat_svm,cMat_gbm, file = save_path)
