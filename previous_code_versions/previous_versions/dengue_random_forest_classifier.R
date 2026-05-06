###Dengue patient classifier using Random Forest with anti-overfitting measures
###If needed install packages before running (use commented code below):
###install.packages(c("tidyverse","caret","randomForest","pROC","ggplot2","gridExtra","Boruta","ranger"))

library(tidyverse)
library(caret)
library(randomForest)
library(ranger)
library(pROC)
library(ggplot2)
library(gridExtra)
library(Boruta)
library(treeshap)

###Random seed for reproducibility
set.seed(123)

###Data preparation (basically subsetting and modifying column names)

df <- read.csv("C:/dengue/LFQ_DENV_AP_CP_imputed_corrected.csv", header = TRUE)
df2 <- read.csv("C:/dengue/DENV_filtered_imputed_dataset.csv",
                header = TRUE
                )
#TODO: Classifiers for AP vs RP
#Transpose columns to have Group as columns,
#Genes as index and LFQ as values

df2 <- df2 %>% filter(Group %in% c("AP","RP"))

df_wide <- df2 %>%
  pivot_wider(names_from = Group, values_from = LFQ)

colnames(df) <- c("Genes","patient_id","timepoint","Run","class","LFQ")

###If required, subset data (perhaps not required in this case since data was already
###subsetted) (uncomment only if required) (also can subset for multiple timepoints)

df_sub <- df %>% filter(timepoint %in% c("AP"))

# Pivot to wide format to get the genes on the columns and the patients in the rows

df_wide <- df_sub %>%
  select(patient_id, timepoint, Genes, LFQ, class) %>%
  pivot_wider(
    id_cols = c(patient_id, timepoint, class),
    names_from = Genes,
    values_from = LFQ
  )

#Check dimensions if required
cat("Data dimensions:", nrow(df_wide), "observations x", ncol(df_wide), "columns\n\n")

### 2.Check for data leakage
###Why? Because data that "leaks" into the training set from outside during training can
###give us an inflated performance

# Check for patients with multiple timepoints
patient_counts <- df_wide %>%
  group_by(patient_id) %>%
  summarise(
    n_observations = n(),
    n_timepoints = n_distinct(timepoint),
    classes = paste(unique(class), collapse = ", ")
  ) %>%
  arrange(desc(n_observations))

cat("Total unique patients:", n_distinct(df_wide$patient_id), "\n")
cat("Total observations:", nrow(df_wide), "\n")
cat("Patients with multiple timepoints:", sum(patient_counts$n_timepoints > 1), "\n\n")

# Show patients with multiple observations
if(sum(patient_counts$n_timepoints > 1) > 0) {
  cat("The following patients have multiple timepoints:\n")
  print(patient_counts %>% filter(n_timepoints > 1))
  cat("Must split at PATIENT level to avoid leakage\n\n")
}

# Check class consistency per patient
class_check <- df_wide %>%
  group_by(patient_id) %>%
  summarise(n_classes = n_distinct(class)) %>%
  filter(n_classes > 1)

if(nrow(class_check) > 0) {
  cat("Some patients have inconsistent class labels!\n")
  print(class_check)
  cat("\n")
} else {
  cat("All patients have consistent class labels across timepoints\n\n")
}

##Creating a patient-level stratified train/test split

# Get unique patients with their class
patient_classes <- df_wide %>%
  group_by(patient_id) %>%
  summarise(class = first(class)) %>%
  ungroup()

cat("Class distribution across patients:\n")
print(table(patient_classes$class))
cat("\n")

# Stratified split at PATIENT level
train_patients_idx <- createDataPartition(
  patient_classes$class,
  p = 0.75,  # 75-25 split
  list = FALSE
)[,1]

train_patients <- patient_classes$patient_id[train_patients_idx]
test_patients <- patient_classes$patient_id[-train_patients_idx]

cat("Training set:", length(train_patients), "patients\n")
cat("Test set:", length(test_patients), "patients\n\n")

# Split data based on patient membership
train_data <- df_wide %>% filter(patient_id %in% train_patients)
test_data <- df_wide %>% filter(patient_id %in% test_patients)

cat("Training observations:", nrow(train_data), "\n")
cat("Test observations:", nrow(test_data), "\n")
cat("Training class distribution:\n")
print(table(train_data$class))
cat("Test class distribution:\n")
print(table(test_data$class))
cat("\n")

###Checking data leakage again: make sure there are no overlapping patients
###If this fails, go back to line #53-57

overlap <- intersect(train_patients, test_patients)
if(length(overlap) > 0) {
  stop("Patients found in both train and test sets!")
} else {
  cat("No patient overlap between train and test sets\n\n")
}

###Creating the train and test matrices
###In case there are any errors in the dimensions, then we need to check
###data structures that were passed

y_train <- factor(train_data$class)
y_test <- factor(test_data$class)

X_train <- train_data %>% select(-patient_id, -timepoint, -class) %>% as.matrix()
X_test <- test_data %>% select(-patient_id, -timepoint, -class) %>% as.matrix()

# Check for missing values
na_count_train <- sum(is.na(X_train))
na_count_test <- sum(is.na(X_test))

cat("Missing values in training set:", na_count_train,
    "(", round(100*na_count_train/length(X_train), 2), "%)\n")
cat("Missing values in test set:", na_count_test,
    "(", round(100*na_count_test/length(X_test), 2), "%)\n")

# Impute missing values with 0
X_train[is.na(X_train)] <- 0
X_test[is.na(X_test)] <- 0

cat("Missing values imputed to 0\n\n")

###FEATURE SELECTION TO REDUCE OVERFITTING
###Using multiple methods to select the most important features

cat("=== FEATURE SELECTION ===\n")
cat("Original number of features:", ncol(X_train), "\n")

# Method 1: Variance-based filtering (remove low-variance features)
feature_var <- apply(X_train, 2, var)
variance_threshold <- quantile(feature_var, 0.1)  # Keep top 90% by variance
high_var_features <- feature_var > variance_threshold
X_train_var <- X_train[, high_var_features]
X_test_var <- X_test[, high_var_features]

cat("After variance filtering:", ncol(X_train_var), "features\n")

# Method 2: Correlation-based filtering (remove highly correlated features)
cor_matrix <- cor(X_train_var)
high_cor_pairs <- which(abs(cor_matrix) > 0.95 & cor_matrix != 1, arr.ind = TRUE)
features_to_remove <- unique(high_cor_pairs[, 2])
if(length(features_to_remove) > 0) {
  X_train_cor <- X_train_var[, -features_to_remove]
  X_test_cor <- X_test_var[, -features_to_remove]
} else {
  X_train_cor <- X_train_var
  X_test_cor <- X_test_var
}

cat("After correlation filtering:", ncol(X_train_cor), "features\n")

# Method 3: Random Forest importance (select top features)
rf_data <- data.frame(y = y_train, X_train_cor)
rf_model <- randomForest(y ~ ., data = rf_data, ntree = 100, importance = TRUE)
rf_importance <- randomForest::importance(rf_model)[, "MeanDecreaseGini"]
top_features_rf <- names(sort(rf_importance, decreasing = TRUE)[1:min(100, length(rf_importance))])

X_train_rf <- X_train_cor[, top_features_rf]
X_test_rf <- X_test_cor[, top_features_rf]

cat("After Random Forest selection:", ncol(X_train_rf), "features\n")

# Method 4: Boruta feature selection (if computationally feasible)
if(ncol(X_train_rf) <= 200) {  # Only run Boruta if not too many features
  cat("Running Boruta feature selection...\n")
  boruta_data <- data.frame(y = y_train, X_train_rf)
  boruta_result <- Boruta(y ~ ., data = boruta_data, doTrace = 0, maxRuns = 50)
  confirmed_features <- names(boruta_result$finalDecision[boruta_result$finalDecision == "Confirmed"])
  
  if(length(confirmed_features) > 0) {
    X_train_final <- X_train_rf[, confirmed_features]
    X_test_final <- X_test_rf[, confirmed_features]
    cat("After Boruta selection:", ncol(X_train_final), "features\n")
  } else {
    X_train_final <- X_train_rf
    X_test_final <- X_test_rf
    cat("Boruta found no confirmed features, using RF selection\n")
  }
} else {
  X_train_final <- X_train_rf
  X_test_final <- X_test_rf
  cat("Skipping Boruta due to high dimensionality\n")
}

cat("Final feature count:", ncol(X_train_final), "\n\n")

###Now, we need to handle class imbalance because the number of
###data-points in serious and subclinical are highly unequal
###We think using SMOTE-like upsampling could help with this

#First check the training counts before balancing to get a better idea of
#the scale of the imbalance

train_tbl <- table(y_train)
cat("Class counts before balancing:\n")
print(train_tbl)
cat("\n")

if(diff(range(train_tbl)) > 0) {
  maj_class <- names(train_tbl)[which.max(train_tbl)]
  min_class <- names(train_tbl)[which.min(train_tbl)]
  n_maj <- max(train_tbl)
  
  idx_maj <- which(y_train == maj_class)
  idx_min <- which(y_train == min_class)
  
  #The sampling happens below
  
  # Upsample minority class
  idx_min_up <- sample(idx_min, size = n_maj, replace = TRUE)
  idx_bal <- c(idx_maj, idx_min_up)
  
  X_train_bal <- X_train_final[idx_bal, ]
  y_train_bal <- y_train[idx_bal]
  
  # Shuffle to mix classes
  perm <- sample(nrow(X_train_bal))
  X_train_bal <- X_train_bal[perm, ]
  y_train_bal <- y_train_bal[perm]
  
  cat("Class counts after upsampling:\n")
  print(table(y_train_bal))
  cat("\n")
} else {
  cat("Classes already balanced\n\n")
  X_train_bal <- X_train_final
  y_train_bal <- y_train
}

###Feature scaling (optional for Random Forest, but included for consistency)
###Random Forest is generally robust to feature scaling, but we include it for consistency

# Calculate scaling parameters ONLY on training data
scaler_center <- apply(X_train_bal, 2, mean)
scaler_scale <- apply(X_train_bal, 2, sd)
scaler_scale[scaler_scale == 0] <- 1  # Avoid division by zero

# Apply scaling
X_train_scaled <- scale(X_train_bal, center = scaler_center, scale = scaler_scale)
X_test_scaled <- scale(X_test_final, center = scaler_center, scale = scaler_scale)

cat("Features scaled using training set parameters\n")
cat("Training set dimensions:", nrow(X_train_scaled), "x", ncol(X_train_scaled), "\n")
cat("Test set dimensions:", nrow(X_test_scaled), "x", ncol(X_test_scaled), "\n\n")

###RANDOM FOREST HYPERPARAMETER TUNING
###Using comprehensive hyperparameter tuning for Random Forest

cat("=== RANDOM FOREST HYPERPARAMETER TUNING ===\n")
cat("Finding optimal Random Forest parameters...\n\n")

# Prepare data for Random Forest
rf_data_scaled <- data.frame(y = y_train_bal, X_train_scaled, check.names = FALSE)

# Define parameter grid for tuning
param_grid <- expand.grid(
  ntree = c(100, 200, 300, 500),
  mtry = c(sqrt(ncol(X_train_scaled)), log2(ncol(X_train_scaled)), ncol(X_train_scaled)/3),
  nodesize = c(1, 3, 5, 10),
  maxnodes = c(10, 20, 50, NULL)
)

# Remove invalid maxnodes values
param_grid$maxnodes[is.na(param_grid$maxnodes)] <- NULL

# Sample a subset of parameters for computational efficiency
set.seed(123)
param_subset <- param_grid[sample(nrow(param_grid), min(30, nrow(param_grid))), ]

cat("Testing", nrow(param_subset), "parameter combinations...\n")

# Cross-validation for parameter tuning
cv_results <- data.frame()

for(i in 1:nrow(param_subset)) {
  params <- list(
    ntree = param_subset$ntree[i],
    mtry = param_subset$mtry[i],
    nodesize = param_subset$nodesize[i],
    maxnodes = param_subset$maxnodes[i]
  )
  
  # Cross-validation using caret
  cv_model <- train(
    y ~ .,
    data = rf_data_scaled,
    method = "rf",
    trControl = trainControl(
      method = "cv",
      number = 5,
      classProbs = TRUE,
      summaryFunction = twoClassSummary
    ),
    tuneGrid = data.frame(mtry = params$mtry),
    ntree = params$ntree,
    nodesize = params$nodesize,
    maxnodes = params$maxnodes,
    importance = TRUE
  )
  
  best_auc <- max(cv_model$results$ROC)
  
  cv_results <- rbind(cv_results, data.frame(
    ntree = params$ntree,
    mtry = params$mtry,
    nodesize = params$nodesize,
    maxnodes = params$maxnodes,
    best_auc = best_auc
  ))
  
  if(i %% 5 == 0) {
    cat("  Processed", i, "/", nrow(param_subset), "combinations\n")
  }
}

# Select best parameters
best_idx <- which.max(cv_results$best_auc)
best_params <- cv_results[best_idx, ]

cat("\nOPTIMAL PARAMETERS:\n")
print(best_params)
cat("\n")

# Train final model with best parameters
rf_model_final <- randomForest(
  y ~ .,
  data = rf_data_scaled,
  ntree = best_params$ntree,
  mtry = best_params$mtry,
  nodesize = best_params$nodesize,
  maxnodes = best_params$maxnodes,
  importance = TRUE,
  proximity = TRUE,
  keep.inbag = TRUE
)

###Now, we do the predictions on the test set which we split earlier

# Predict probabilities
pred_probs <- predict(rf_model_final, X_test_scaled, type = "prob")[, 2]

cat("Predicted probability statistics:\n")
cat("  Min:", round(min(pred_probs), 4), "\n")
cat("  Q1:", round(quantile(pred_probs, 0.25), 4), "\n")
cat("  Median:", round(median(pred_probs), 4), "\n")
cat("  Mean:", round(mean(pred_probs), 4), "\n")
cat("  Q3:", round(quantile(pred_probs, 0.75), 4), "\n")
cat("  Max:", round(max(pred_probs), 4), "\n\n")

###ROC analysis
###Here tp = true positive, fp = false positive,
###tn = true negative and fn = false negative


# ROC curve
roc_obj <- roc(response = y_test, predictor = pred_probs, quiet = TRUE)
auc_val <- auc(roc_obj)
cat("Test Set AUC:", round(auc_val, 4), "\n\n")

# Calculate 95% CI for AUC
ci_obj <- ci.auc(roc_obj)
ci_lower <- ci_obj[1]
ci_upper <- ci_obj[3]

# Manual calculation of Youden's index
# Get all thresholds and their corresponding sens/spec
all_coords <- coords(roc_obj, x = "all", ret = c("threshold", "sensitivity", "specificity"))

# Calculate Youden's J for each threshold
#This value is used to find the optimal threshold for the model, which is the
#point that is furthest from the diagonal line (random guessing)
youden_j <- all_coords$sensitivity + all_coords$specificity - 1

# Find the threshold that maximizes Youden's J
best_idx <- which.max(youden_j)
best_thresh <- all_coords$threshold[best_idx]
best_sens <- all_coords$sensitivity[best_idx]
best_spec <- all_coords$specificity[best_idx]

# Calculate PPV and NPV manually
y_test_binary <- as.numeric(y_test == levels(y_test)[2])
pred_binary <- as.numeric(pred_probs > best_thresh)

tp <- sum(pred_binary == 1 & y_test_binary == 1)
fp <- sum(pred_binary == 1 & y_test_binary == 0)
tn <- sum(pred_binary == 0 & y_test_binary == 0)
fn <- sum(pred_binary == 0 & y_test_binary == 1)

best_ppv <- if((tp + fp) > 0) tp / (tp + fp) else NA
best_npv <- if((tn + fn) > 0) tn / (tn + fn) else NA

cat("Optimal threshold (Youden's J):", round(best_thresh, 4), "\n")
cat("  Youden's J:", round(youden_j[best_idx], 3), "\n")
cat("  Sensitivity:", round(best_sens, 3), "\n")
cat("  Specificity:", round(best_spec, 3), "\n")
if(!is.na(best_ppv)) {
  cat("  PPV:", round(best_ppv, 3), "\n")
}
if(!is.na(best_npv)) {
  cat("  NPV:", round(best_npv, 3), "\n")
}
cat("\nConfusion matrix at optimal threshold:\n")
cat("  TP:", tp, " FP:", fp, "\n")
cat("  FN:", fn, " TN:", tn, "\n\n")

# Predict classes using optimal threshold
pred_class <- ifelse(pred_probs > best_thresh,
                     levels(y_test)[2],
                     levels(y_test)[1])
pred_class <- factor(pred_class, levels = levels(y_test))

###Cross validation using leave-one-patient-out
###Exclude one patient each time from the test set to check the impact on
###model performance

unique_patients <- unique(df_wide$patient_id)
n_patients <- length(unique_patients)

cat("=== LEAVE-ONE-PATIENT-OUT CROSS-VALIDATION ===\n")
cat("Performing LOPO-CV on", n_patients, "patients...\n")

lopo_predictions <- list()
lopo_true_labels <- list()
lopo_patient_ids <- list()

for(i in seq_along(unique_patients)) {
  
  test_patient <- unique_patients[i]
  
  if(i %% 5 == 0) {
    cat("  Processed", i, "/", n_patients, "patients\n")
  }
  
  # Debug information for problematic folds
  if(i <= 3) {
    cat("    Debug - Patient", test_patient, ": Original features =", ncol(X_train), "\n")
  }
  
  
  # Split data
  train_fold <- df_wide %>% filter(patient_id != test_patient)
  test_fold <- df_wide %>% filter(patient_id == test_patient)
  
  # Prepare matrices
  y_train_fold <- factor(train_fold$class)
  y_test_fold <- factor(test_fold$class)
  
  X_train_fold <- train_fold %>%
    select(-patient_id, -timepoint, -class) %>%
    as.matrix()
  X_test_fold <- test_fold %>%
    select(-patient_id, -timepoint, -class) %>%
    as.matrix()
  
  X_train_fold[is.na(X_train_fold)] <- 0
  X_test_fold[is.na(X_test_fold)] <- 0
  
  # Apply same feature selection pipeline to each fold
  # Variance filtering
  feature_var_fold <- apply(X_train_fold, 2, var)
  variance_threshold_fold <- quantile(feature_var_fold, 0.1)
  high_var_features_fold <- feature_var_fold > variance_threshold_fold
  X_train_fold_var <- X_train_fold[, high_var_features_fold, drop = FALSE]
  X_test_fold_var <- X_test_fold[, high_var_features_fold, drop = FALSE]
  
  # Correlation filtering
  if(ncol(X_train_fold_var) > 1) {
    cor_matrix_fold <- cor(X_train_fold_var)
    high_cor_pairs_fold <- which(abs(cor_matrix_fold) > 0.95 & cor_matrix_fold != 1, arr.ind = TRUE)
    features_to_remove_fold <- unique(high_cor_pairs_fold[, 2])
    if(length(features_to_remove_fold) > 0) {
      X_train_fold_cor <- X_train_fold_var[, -features_to_remove_fold, drop = FALSE]
      X_test_fold_cor <- X_test_fold_var[, -features_to_remove_fold, drop = FALSE]
    } else {
      X_train_fold_cor <- X_train_fold_var
      X_test_fold_cor <- X_test_fold_var
    }
  } else {
    X_train_fold_cor <- X_train_fold_var
    X_test_fold_cor <- X_test_fold_var
  }
  
  # Random Forest selection (top 50 features for computational efficiency)
  if(ncol(X_train_fold_cor) > 1) {
    rf_data_fold <- data.frame(y = y_train_fold, X_train_fold_cor)
    rf_model_fold <- randomForest(y ~ ., data = rf_data_fold, ntree = 50, importance = TRUE)
    rf_importance_fold <- randomForest::importance(rf_model_fold)[, "MeanDecreaseGini"]
    top_features_rf_fold <- names(sort(rf_importance_fold, decreasing = TRUE)[1:min(50, length(rf_importance_fold))])
    
    # Ensure features exist in both training and test sets
    common_features <- intersect(top_features_rf_fold, colnames(X_test_fold_cor))
    
    if(length(common_features) > 0) {
      X_train_fold_final <- X_train_fold_cor[, common_features, drop = FALSE]
      X_test_fold_final <- X_test_fold_cor[, common_features, drop = FALSE]
      
      # Debug information for first few folds
      if(i <= 3) {
        cat("    Debug - After RF selection: Features =", length(common_features), "\n")
      }
    } else {
      # Fallback: use all available features if no common features found
      X_train_fold_final <- X_train_fold_cor
      X_test_fold_final <- X_test_fold_cor
      
      if(i <= 3) {
        cat("    Debug - Using all features after RF (no common features found)\n")
      }
    }
  } else {
    X_train_fold_final <- X_train_fold_cor
    X_test_fold_final <- X_test_fold_cor
  }
  
  # Scale using training fold parameters
  scaler_c <- apply(X_train_fold_final, 2, mean)
  scaler_s <- apply(X_train_fold_final, 2, sd)
  scaler_s[scaler_s == 0] <- 1
  
  X_train_fold_sc <- scale(X_train_fold_final, center = scaler_c, scale = scaler_s)
  X_test_fold_sc <- scale(X_test_fold_final, center = scaler_c, scale = scaler_s)
  
  # Handle class imbalance within each LOPO fold
  train_tbl_fold <- table(y_train_fold)
  
  if(diff(range(train_tbl_fold)) > 0) {
    maj_class_fold <- names(train_tbl_fold)[which.max(train_tbl_fold)]
    min_class_fold <- names(train_tbl_fold)[which.min(train_tbl_fold)]
    n_maj_fold <- max(train_tbl_fold)
    
    idx_maj_fold <- which(y_train_fold == maj_class_fold)
    idx_min_fold <- which(y_train_fold == min_class_fold)
    
    # Upsample minority class
    idx_min_up_fold <- sample(idx_min_fold, size = n_maj_fold, replace = TRUE)
    idx_bal_fold <- c(idx_maj_fold, idx_min_up_fold)
    
    X_train_fold_bal <- X_train_fold_sc[idx_bal_fold, ]
    y_train_fold_bal <- y_train_fold[idx_bal_fold]
    
    # Shuffle to mix classes
    perm_fold <- sample(nrow(X_train_fold_bal))
    X_train_fold_bal <- X_train_fold_bal[perm_fold, ]
    y_train_fold_bal <- y_train_fold_bal[perm_fold]
  } else {
    X_train_fold_bal <- X_train_fold_sc
    y_train_fold_bal <- y_train_fold
  }
  
  # Train Random Forest model using best parameters from full training set
  rf_data_fold <- data.frame(y = y_train_fold_bal, X_train_fold_bal)
  
  rf_fold <- randomForest(
    y ~ .,
    data = rf_data_fold,
    ntree = best_params$ntree,
    mtry = best_params$mtry,
    nodesize = best_params$nodesize,
    maxnodes = best_params$maxnodes,
    importance = TRUE
  )
  
  # Predict
  pred_prob_fold <- predict(rf_fold, X_test_fold_sc, type = "prob")[, 2]
  
  # Store results
  lopo_predictions[[i]] <- pred_prob_fold
  lopo_true_labels[[i]] <- as.character(y_test_fold)
  lopo_patient_ids[[i]] <- rep(test_patient, length(y_test_fold))
}

cat("LOPO-CV complete\n\n")

# Aggregate results
all_preds <- unlist(lopo_predictions)
all_labels <- factor(unlist(lopo_true_labels))
all_patient_ids <- unlist(lopo_patient_ids)

# Compute LOPO-CV AUC
lopo_roc <- roc(all_labels, all_preds, quiet = TRUE)
lopo_auc <- auc(lopo_roc)


cat("LOPO-CV Results:\n")

cat("LOPO-CV AUC:", round(lopo_auc, 4), "\n")
cat("Number of predictions:", length(all_preds), "\n")
cat("Class distribution:\n")
print(table(all_labels))
cat("\n")

###Feature importance analysis
###Here we check which genes are more important for the final prediction

###Feature importance

importance_matrix <- randomForest::importance(rf_model_final)
importance_df <- data.frame(
  Feature = rownames(importance_matrix),
  MeanDecreaseGini = importance_matrix[, "MeanDecreaseGini"],
  MeanDecreaseAccuracy = importance_matrix[, "MeanDecreaseAccuracy"]
) %>%
  arrange(desc(MeanDecreaseGini))

cat("Number of selected features:", nrow(importance_df), "\n\n")

if(nrow(importance_df) > 0) {
  cat("Top 20 features by importance:\n")
  print(head(importance_df, 20))
  cat("\n")
  
  # Count positive and negative coefficients (for Random Forest, we look at feature importance)
  cat("Feature importance analysis completed\n\n")
}

###SHAP value analysis using treeshap
cat("=== SHAP VALUE ANALYSIS (severe vs subclinical) ===\n")

shap_training_df <- as.data.frame(X_train_scaled, check.names = FALSE)
colnames(shap_training_df) <- colnames(X_train_scaled)

rf_model_unified <- treeshap::randomForest.unify(rf_model_final, data = shap_training_df)
rf_model_unified <- treeshap::set_reference_dataset(rf_model_unified, shap_training_df)
shap_result <- treeshap::treeshap(rf_model_unified, shap_training_df)

shap_importance <- colMeans(abs(shap_result$shaps))
shap_importance_df <- data.frame(
  Feature = names(shap_importance),
  MeanAbsSHAP = shap_importance
) %>%
  arrange(desc(MeanAbsSHAP))

cat("Top genes ranked by mean |SHAP| values:\n")
print(head(shap_importance_df, 20))
cat("\n")

if(nrow(shap_importance_df) > 0) {
  top_shap_features <- head(shap_importance_df, 20)
  
  shap_plot <- ggplot(top_shap_features, aes(x = reorder(Feature, MeanAbsSHAP), y = MeanAbsSHAP)) +
    geom_col(fill = "#4DBBD5") +
    coord_flip() +
    labs(
      title = "Top 20 SHAP Feature Importances (Severe vs Subclinical)",
      x = "Gene",
      y = "Mean |SHAP|"
    ) +
    theme_minimal(base_size = 12)
  
  print(shap_plot)
}

###Checking for overfitting

# Predict on training set
pred_probs_train <- predict(rf_model_final, X_train_scaled, type = "prob")[, 2]

# Get a subset that matches original training patients (before upsampling)
# Use only original training indices
original_train_idx <- 1:nrow(X_train_final)
pred_probs_train_orig <- predict(rf_model_final, scale(X_train_final, center = scaler_center, scale = scaler_scale), type = "prob")[, 2]

train_roc <- roc(y_train, pred_probs_train_orig, quiet = TRUE)
train_auc <- auc(train_roc)

cat("=== MODEL PERFORMANCE SUMMARY ===\n")
cat("Training Set AUC:", round(train_auc, 4), "\n")
cat("Test Set AUC:", round(auc_val, 4), "\n")
cat("LOPO-CV AUC:", round(lopo_auc, 4), "\n")
cat("Difference (Train - Test):", round(train_auc - auc_val, 4), "\n")
cat("Difference (Train - LOPO-CV):", round(train_auc - lopo_auc, 4), "\n\n")

if(abs(train_auc - auc_val) > 0.15) {
  cat("Large gap between train and test AUC suggests overfitting\n\n")
} else if(abs(train_auc - lopo_auc) > 0.15) {
  cat("Large gap between train and LOPO-CV AUC suggests overfitting\n\n")
} else {
  cat("Train, test, and CV performance are similar - overfitting controlled!\n\n")
}

###Visualizing the various metrics and results

# Set up plotting layout
par(mfrow = c(2, 2))

# Plot 1: ROC Curve
plot(roc_obj, col = "blue", lwd = 2,
     main = paste0("ROC Curve (AUC = ", round(auc_val, 3), ")"),
     cex.main = 1.2)
abline(a = 0, b = 1, lty = 2, col = "gray")
legend("bottomright",
       legend = c(paste0("Test AUC = ", round(auc_val, 3)),
                  paste0("95% CI: [", round(ci_lower, 3), ", ", round(ci_upper, 3), "]")),
       bty = "n")

# Plot 2: LOPO-CV ROC
plot(lopo_roc, col = "darkgreen", lwd = 2,
     main = paste0("LOPO-CV ROC (AUC = ", round(lopo_auc, 3), ")"),
     cex.main = 1.2)
abline(a = 0, b = 1, lty = 2, col = "gray")

# Plot 3: Hyperparameter tuning results
plot(cv_results$best_auc, 
     type = "b", pch = 19, col = "purple", lwd = 2,
     xlab = "Parameter Combination", ylab = "CV AUC",
     main = "Random Forest Hyperparameter Tuning",
     cex.main = 1.2)
abline(h = best_params$best_auc, lty = 2, col = "red")
grid()

# Plot 4: Train vs Test comparison
plot(c(1, 2, 3), c(train_auc, auc_val, lopo_auc),
     type = "b", pch = 19, cex = 2, col = c("red", "blue", "darkgreen"),
     ylim = c(0.5, 1.0), xlim = c(0.5, 3.5),
     xaxt = "n", xlab = "", ylab = "AUC",
     main = "Model Performance Comparison",
     cex.main = 1.2)
axis(1, at = c(1, 2, 3), labels = c("Train", "Test", "LOPO-CV"))
abline(h = 0.5, lty = 2, col = "gray")
grid()

par(mfrow = c(1, 1))

# Plot 5: Feature importance (top 20)
if(nrow(importance_df) > 0) {
  top_features <- head(importance_df, 20)
  
  p1 <- ggplot(top_features, aes(x = reorder(Feature, MeanDecreaseGini), y = MeanDecreaseGini)) +
    geom_col(fill = "#4DBBD5") +
    coord_flip() +
    labs(title = "Top 20 Feature Importance (Random Forest)",
         subtitle = paste0("Selected from ", nrow(importance_df), " features"),
         x = "Gene", y = "Mean Decrease Gini") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
  
  print(p1)
}

# Plot 6: Probability distribution
prob_df <- data.frame(
  Probability = pred_probs,
  True_Class = as.character(y_test),
  Patient_ID = test_data$patient_id
)

p2 <- ggplot(prob_df, aes(x = Probability, fill = True_Class)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = best_thresh, linetype = "dashed",
             color = "red", linewidth = 1) +
  labs(title = "Predicted Probability Distribution",
       subtitle = paste0("Test set (", length(test_patients), " patients, ",
                         nrow(test_data), " observations)"),
       x = "Predicted Probability",
       fill = "True Class") +
  theme_minimal(base_size = 12) +
  scale_fill_manual(values = c("severe" = "#E64B35", "subclinical" = "#4DBBD5")) +
  annotate("text", x = best_thresh, y = Inf,
           label = paste0("Threshold = ", round(best_thresh, 3)),
           hjust = -0.1, vjust = 2, color = "red", size = 3.5)

print(p2)

# Plot 7: Prediction confidence by patient
prob_summary <- prob_df %>%
  group_by(Patient_ID, True_Class) %>%
  summarise(
    Mean_Prob = mean(Probability),
    SD_Prob = sd(Probability),
    N_Obs = n(),
    .groups = "drop"
  ) %>%
  arrange(Mean_Prob)

p3 <- ggplot(prob_summary, aes(x = reorder(Patient_ID, Mean_Prob),
                               y = Mean_Prob, color = True_Class)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = pmax(0, Mean_Prob - SD_Prob),
                    ymax = pmin(1, Mean_Prob + SD_Prob)),
                width = 0.2) +
  geom_hline(yintercept = best_thresh, linetype = "dashed", color = "red") +
  labs(title = "Per-Patient Prediction Confidence",
       subtitle = "Mean probability ± SD across timepoints",
       x = "Patient ID", y = "Predicted Probability",
       color = "True Class") +
  theme_minimal(base_size = 10) +
  scale_color_manual(values = c("severe" = "#E64B35", "subclinical" = "#4DBBD5")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p3)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Optimal hyperparameters:\n")
print(best_params)
cat("Final LOPO-CV AUC:", round(lopo_auc, 4), "\n")
cat("Features reduced from", ncol(X_train), "to", ncol(X_train_final), "\n")
