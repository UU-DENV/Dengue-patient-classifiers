### ============================================================================
### UNIFIED DENGUE PATIENT CLASSIFIER - RANDOM FOREST
### ============================================================================
### This unified classifier can perform multiple classification tasks:
###   1. AP vs RP (Acute Phase vs Recovery Phase)
###   2. DF vs DHF (Dengue Fever vs Dengue Hemorrhagic Fever)
###   3. Hospitalized vs Subclinical
###   4. DF vs DHF with 18-gene panel
###   5. Hospitalized vs Subclinical with 21-gene panel
###
### KEY ANTI-OVERFITTING MEASURES:
###   - NO UPSAMPLING: Uses class weights instead of upsampling
###   - AGGRESSIVE FEATURE SELECTION: Moderate feature count with n_samples/3 rule
###   - CONSERVATIVE HYPERPARAMETERS: Patient-grouped CV for hyperparameter tuning
###   - CV ON ORIGINAL DATA: Hyperparameter tuning uses original (non-upsampled) data
###   - FEATURE ENGINEERING: Creates biologically meaningful ratios and interactions
###   - BATCH EFFECT VALIDATION: Checks for normalization issues and batch effects
###
### ============================================================================
### USER CONFIGURATION SECTION
### ============================================================================
### MODIFY THESE SETTINGS TO SELECT THE CLASSIFICATION TASK

# Classification type: "AP_RP", "DF_DHF", "HP_SUB", "DF_DHF_panel18", "HP_SUB_panel21"
CLASSIFICATION_TYPE <- "HP_SUB"  # THE USER MUST SELECT ONE CLASSIFICATION TYPE FROM THE LIST ABOVE

# Input file path (can be CSV or TSV)
# For AP_RP: Use "DENV_filtered_imputed_dataset.csv"
# For DF_DHF: Use "260108_LFQ_DENV_AP.csv"
# For HP_SUB: Use "LFQ_DENV_AP_CP_imputed_corrected.csv"
INPUT_FILE <- "C:/dengue/LFQ_DENV_AP_CP_imputed_corrected.csv"  # EXAMPLE FILE

# Output prefix for all generated files
OUTPUT_PREFIX <- "unified_RF"  # TO BE MODIFIED AS PREFERRED

# ============================================================================
# END USER CONFIGURATION
# ============================================================================

### Load required packages
if(!require(tidyverse)) install.packages("tidyverse")
if(!require(caret)) install.packages("caret")
if(!require(randomForest)) install.packages("randomForest")
if(!require(ranger)) install.packages("ranger")
if(!require(pROC)) install.packages("pROC")
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(gridExtra)) install.packages("gridExtra")
if(!require(Boruta)) install.packages("Boruta")
if(!require(treeshap)) install.packages("treeshap")
if(!require(extrafont)) install.packages("extrafont")
# Optional: for EMF export on Windows (install.packages("devEMF") if needed)
has_devEMF <- requireNamespace("devEMF", quietly = TRUE)
if (has_devEMF) suppressPackageStartupMessages(library(devEMF))

library(tidyverse)
library(caret)
library(randomForest)
library(ranger)
library(pROC)
library(ggplot2)
library(gridExtra)
library(Boruta)
library(treeshap)
library(extrafont)

### Set random seed for reproducibility
set.seed(123)

### Load and register fonts for ggplot2
# Helper function to create element_text with optional font family
make_text_element <- function(size, face = NULL, family = NULL) {
  args <- list(size = size)
  if(!is.null(face)) args$face <- face
  if(!is.null(family)) args$family <- family
  do.call(element_text, args)
}

# Try to use Arial, but don't fail if it's not available
FONT_FAMILY <- NULL
tryCatch({
  if(.Platform$OS.type == "windows") {
    # On Windows, try to use Arial directly (system font)
    # If fonts are imported via extrafont, use that
    if(length(fonts()) > 0 && "Arial" %in% fonts()) {
      loadfonts(device = "win", quiet = TRUE)
      FONT_FAMILY <- "Arial"
      cat("Using Arial font for plots\n")
    } else {
      # Fonts not imported - will use system default (typically Arial on Windows)
      cat("Using system default font (typically Arial on Windows)\n")
    }
  }
}, error = function(e) {
  cat("Note: Using system default font\n")
})

### Create output directory
if(!dir.exists("figures")) dir.create("figures")

### ============================================================================
### MODULE 1: DATA LOADING AND TRANSFORMATION
### ============================================================================

load_and_transform_data <- function(input_file, classification_type) {
  cat("=== MODULE 1: DATA LOADING AND TRANSFORMATION ===\n")
  cat("Classification type:", classification_type, "\n")
  cat("Input file:", input_file, "\n\n")
  
  # Read data (handle both CSV and TSV)
  if(grepl("\\.tsv$|\\.txt$", input_file, ignore.case = TRUE)) {
    df <- read.delim(input_file, header = TRUE, sep = "\t")
  } else {
    df <- read.csv(input_file, header = TRUE)
  }
  
  cat("Original data dimensions:", nrow(df), "rows x", ncol(df), "columns\n")
  cat("Column names:\n")
  print(colnames(df))
  cat("\n")
  
  # Transform data based on classification type
  if(classification_type == "AP_RP") {
    # AP vs RP: Can use multiple file formats
    df <- df %>% filter(Group %in% c("AP", "RP"))
    df$class <- df$Group
    
    # Handle different column name formats
    if("Subject" %in% colnames(df)) {
      df <- df %>% rename(patient_id = Subject)
    }
    
    # Check if log2_LFQ exists, otherwise use log_LFQ or LFQ
    if("log2_LFQ" %in% colnames(df)) {
      df$LFQ <- df$log2_LFQ
    } else if("log_LFQ" %in% colnames(df)) {
      df$LFQ <- df$log_LFQ
    } else if(!"LFQ" %in% colnames(df)) {
      stop("AP_RP: No LFQ column found. Expected 'LFQ', 'log_LFQ', or 'log2_LFQ'")
    }
    
    df_wide <- df %>%
      select(patient_id, Genes, LFQ, class) %>%
      pivot_wider(
        id_cols = c(patient_id, class),
        names_from = Genes,
        values_from = LFQ,
        values_fn = mean  # Aggregate duplicates by taking mean
      )
    
    class_levels <- c("AP", "RP")
    class_labels <- c("AP", "RP")
    
  } else if(classification_type == "DF_DHF") {
    # DF vs DHF: Use 260108_LFQ_DENV_AP.csv
    df <- df %>% 
      filter(Group %in% c("AP")) %>%
      filter(Classification_02 %in% c("DF", "DHF"))
    
    df$class <- df$Classification_02
    
    # Check if log2_LFQ exists, otherwise use log_LFQ
    if("log2_LFQ" %in% colnames(df)) {
      df$LFQ <- df$log2_LFQ
    } else if("log_LFQ" %in% colnames(df)) {
      df$LFQ <- df$log_LFQ
    }
    
    df_wide <- df %>%
      select(Subject, Group, Genes, LFQ, class) %>%
      pivot_wider(
        id_cols = c(Subject, Group, class),
        names_from = Genes,
        values_from = LFQ
      )
    
    df_wide <- df_wide %>% rename(patient_id = Subject)
    class_levels <- c("DF", "DHF")
    class_labels <- c("DF", "DHF")
    
  } else if(classification_type == "HP_SUB") {
    # Hospitalized vs Subclinical: Can use multiple file formats
    # Check if file has the expected columns
    if("class" %in% colnames(df) && "patient_id" %in% colnames(df)) {
      # Format 1: LFQ_DENV_AP_CP_imputed_corrected.csv format
      if("timepoint" %in% colnames(df)) {
        colnames(df) <- c("Genes", "patient_id", "timepoint", "Run", "class", "LFQ")
        df_sub <- df %>% filter(timepoint %in% c("AP"))
        
        df_wide <- df_sub %>%
          select(patient_id, timepoint, Genes, LFQ, class) %>%
          pivot_wider(
            id_cols = c(patient_id, timepoint, class),
            names_from = Genes,
            values_from = LFQ
          )
      } else {
        # Already in correct format
        df_wide <- df %>%
          select(patient_id, Genes, LFQ, class) %>%
          pivot_wider(
            id_cols = c(patient_id, class),
            names_from = Genes,
            values_from = LFQ
          )
      }
    } else if("Subject" %in% colnames(df) && "Classification" %in% colnames(df)) {
      # Format 2: 20260120_LFQ_DENV_AP.csv format (similar to DF_DHF)
      df <- df %>% filter(Group %in% c("AP"))
      
      # Use Classification column for class (hospitalized vs subclinical)
      df$class <- df$Classification
      
      # Check if log2_LFQ exists, otherwise use LFQ
      if("log2_LFQ" %in% colnames(df)) {
        df$LFQ <- df$log2_LFQ
      } else if("log_LFQ" %in% colnames(df)) {
        df$LFQ <- df$log_LFQ
      }
      
      df_wide <- df %>%
        select(Subject, Group, Genes, LFQ, class) %>%
        pivot_wider(
          id_cols = c(Subject, Group, class),
          names_from = Genes,
          values_from = LFQ
        )
      
      df_wide <- df_wide %>% rename(patient_id = Subject)
    } else {
      stop("HP_SUB: File format not recognized. Expected columns: 'patient_id'/'Subject' and 'class'/'Classification'")
    }
    
    class_levels <- c("hospitalized", "subclinical")
    class_labels <- c("Hospitalized", "Subclinical")
    
  } else if(classification_type == "DF_DHF_panel18") {
    # DF vs DHF with 18-gene panel
    df <- df %>% 
      filter(Group %in% c("AP")) %>%
      filter(Classification_02 %in% c("DF", "DHF"))
    
    df$class <- df$Classification_02
    
    # Check if log2_LFQ exists, otherwise use log_LFQ
    if("log2_LFQ" %in% colnames(df)) {
      df$LFQ <- df$log2_LFQ
    } else if("log_LFQ" %in% colnames(df)) {
      df$LFQ <- df$log_LFQ
    }
    
    df_wide <- df %>%
      select(Subject, Group, Genes, LFQ, class) %>%
      pivot_wider(
        id_cols = c(Subject, Group, class),
        names_from = Genes,
        values_from = LFQ
      )
    
    df_wide <- df_wide %>% rename(patient_id = Subject)
    
    # Define 18-gene panel
    panel18_genes <- c(
      "CPB2", "CFP", "IGHA1", "COLEC11", "APOE",
      "LCP1", "FCN2", "CCT2", "CCT3", "CDC42",
      "PSMB9", "SAA1", "PTX3", "PTPRK", "HNRNPA2B1",
      "CYCS", "CLEC11A"
    )
    
    # Subset to panel genes
    available_genes <- intersect(panel18_genes, colnames(df_wide))
    if(length(available_genes) > 0) {
      df_wide <- df_wide %>%
        select(patient_id, Group, class, all_of(available_genes))
      cat("Subsetted to", length(available_genes), "panel genes\n")
    }
    
    class_levels <- c("DF", "DHF")
    class_labels <- c("DF", "DHF")
    
  } else if(classification_type == "HP_SUB_panel21") {
    # Hospitalized vs Subclinical with 21-gene panel
    # Check if file has the expected columns
    if("class" %in% colnames(df) && "patient_id" %in% colnames(df)) {
      # Format 1: LFQ_DENV_AP_CP_imputed_corrected.csv format
      if("timepoint" %in% colnames(df)) {
        colnames(df) <- c("Genes", "patient_id", "timepoint", "Run", "class", "LFQ")
        df_sub <- df %>% filter(timepoint %in% c("AP"))
        
        df_wide <- df_sub %>%
          select(patient_id, timepoint, Genes, LFQ, class) %>%
          pivot_wider(
            id_cols = c(patient_id, timepoint, class),
            names_from = Genes,
            values_from = LFQ
          )
      } else {
        # Already in correct format
        df_wide <- df %>%
          select(patient_id, Genes, LFQ, class) %>%
          pivot_wider(
            id_cols = c(patient_id, class),
            names_from = Genes,
            values_from = LFQ
          )
      }
    } else if("Subject" %in% colnames(df) && "Classification" %in% colnames(df)) {
      # Format 2: 20260120_LFQ_DENV_AP.csv format (similar to DF_DHF)
      df <- df %>% filter(Group %in% c("AP"))
      
      # Use Classification column for class (hospitalized vs subclinical)
      df$class <- df$Classification
      
      # Check if log2_LFQ exists, otherwise use LFQ
      if("log2_LFQ" %in% colnames(df)) {
        df$LFQ <- df$log2_LFQ
      } else if("log_LFQ" %in% colnames(df)) {
        df$LFQ <- df$log_LFQ
      }
      
      df_wide <- df %>%
        select(Subject, Group, Genes, LFQ, class) %>%
        pivot_wider(
          id_cols = c(Subject, Group, class),
          names_from = Genes,
          values_from = LFQ
        )
      
      df_wide <- df_wide %>% rename(patient_id = Subject)
    } else {
      stop("HP_SUB_panel21: File format not recognized. Expected columns: 'patient_id'/'Subject' and 'class'/'Classification'")
    }
    
    # Define 21-gene panel (user-specified; 22 genes listed for flexibility)
    panel21_genes <- c(
      "CD14", "CFB", "LGALS3BP", "SERPINA3", "B2M",
      "GOLM1", "SERPING1", "PVR", "LBP", "SERPINF2",
      "ORM1", "WARS1", "LCP1", "NCAM1", "ITIH2", "PCSK9",
      "AMBP", "HRG", "F13B", "AZGP1", "CFP", "APOH"
    )
    
    # Subset to panel genes
    available_genes <- intersect(panel21_genes, colnames(df_wide))
    if(length(available_genes) > 0) {
      # Keep patient_id, class, and optionally timepoint/Group
      keep_cols <- c("patient_id", "class")
      if("timepoint" %in% colnames(df_wide)) keep_cols <- c(keep_cols, "timepoint")
      if("Group" %in% colnames(df_wide)) keep_cols <- c(keep_cols, "Group")
      
      df_wide <- df_wide %>%
        select(all_of(keep_cols), all_of(available_genes))
      cat("Subsetted to", length(available_genes), "panel genes\n")
    }
    
    class_levels <- c("hospitalized", "subclinical")
    class_labels <- c("Hospitalized", "Subclinical")
    
  } else {
    stop("Invalid CLASSIFICATION_TYPE. Must be one of: AP_RP, DF_DHF, HP_SUB, DF_DHF_panel18, HP_SUB_panel21")
  }
  
  cat("Transformed data dimensions:", nrow(df_wide), "observations x", ncol(df_wide), "columns\n")
  cat("Class levels:", paste(class_levels, collapse = ", "), "\n\n")
  
  return(list(
    df_wide = df_wide,
    class_levels = class_levels,
    class_labels = class_labels,
    classification_type = classification_type
  ))
}

### ============================================================================
### MODULE 2: DATA LEAKAGE CHECK AND TRAIN/TEST SPLIT
### ============================================================================

check_leakage_and_split <- function(df_wide, classification_type) {
  cat("=== MODULE 2: DATA LEAKAGE CHECK AND TRAIN/TEST SPLIT ===\n")
  
  # Determine which columns to exclude (patient_id, class, and optionally timepoint/Group)
  exclude_cols <- c("patient_id", "class")
  if("timepoint" %in% colnames(df_wide)) exclude_cols <- c(exclude_cols, "timepoint")
  if("Group" %in% colnames(df_wide)) exclude_cols <- c(exclude_cols, "Group")
  
  # Check for patients with multiple observations
  patient_counts <- df_wide %>%
    group_by(patient_id) %>%
    summarise(
      n_observations = n(),
      classes = paste(unique(class), collapse = ", ")
    ) %>%
    arrange(desc(n_observations))
  
  cat("Total unique patients:", n_distinct(df_wide$patient_id), "\n")
  cat("Total observations:", nrow(df_wide), "\n")
  cat("Patients with multiple observations:", sum(patient_counts$n_observations > 1), "\n\n")
  
  # Check class consistency per patient
  class_check <- df_wide %>%
    group_by(patient_id) %>%
    summarise(n_classes = n_distinct(class)) %>%
    filter(n_classes > 1)
  
  if(nrow(class_check) > 0) {
    cat("WARNING: Some patients have inconsistent class labels!\n")
    print(class_check)
    cat("\n")
  } else {
    cat("All patients have consistent class labels\n\n")
  }
  
  # Patient-level stratified train/test split
  patient_classes <- df_wide %>%
    group_by(patient_id) %>%
    summarise(class = first(class)) %>%
    ungroup()
  
  cat("Class distribution across patients:\n")
  print(table(patient_classes$class))
  cat("\n")
  
  train_patients_idx <- createDataPartition(
    patient_classes$class,
    p = 0.75,
    list = FALSE
  )[,1]
  
  train_patients <- patient_classes$patient_id[train_patients_idx]
  test_patients <- patient_classes$patient_id[-train_patients_idx]
  
  train_data <- df_wide %>% filter(patient_id %in% train_patients)
  test_data <- df_wide %>% filter(patient_id %in% test_patients)
  
  # Verify no overlap
  overlap <- intersect(train_patients, test_patients)
  if(length(overlap) > 0) {
    stop("ERROR: Patients found in both train and test sets!")
  }
  
  cat("Train/test split complete\n")
  cat("Training set:", length(train_patients), "patients,", nrow(train_data), "observations\n")
  cat("Test set:", length(test_patients), "patients,", nrow(test_data), "observations\n\n")
  
  return(list(
    train_data = train_data,
    test_data = test_data,
    train_patients = train_patients,
    test_patients = test_patients
  ))
}

### ============================================================================
### MODULE 3: DATA PREPROCESSING AND IMPUTATION
### ============================================================================

preprocess_data <- function(train_data, test_data, exclude_cols) {
  cat("=== MODULE 3: DATA PREPROCESSING ===\n")
  
  y_train <- factor(train_data$class)
  y_test <- factor(test_data$class)
  
  X_train <- train_data %>% select(-all_of(exclude_cols)) %>% as.matrix()
  X_test <- test_data %>% select(-all_of(exclude_cols)) %>% as.matrix()
  
  # Check for missing values
  na_count_train <- sum(is.na(X_train))
  na_count_test <- sum(is.na(X_test))
  
  cat("Missing values in training set:", na_count_train,
      "(", round(100*na_count_train/length(X_train), 2), "%)\n")
  cat("Missing values in test set:", na_count_test,
      "(", round(100*na_count_test/length(X_test), 2), "%)\n")
  
  # Proteomics-aware imputation (min/half-min strategy)
  cat("Using proteomics-aware imputation (min/half-min strategy)...\n")
  
  for(j in 1:ncol(X_train)) {
    feature_values <- X_train[, j]
    non_zero_values <- feature_values[!is.na(feature_values) & feature_values > 0]
    
    if(length(non_zero_values) > 0) {
      min_val <- min(non_zero_values)
      half_min <- min_val / 2
      X_train[is.na(X_train[, j]), j] <- half_min
      X_test[is.na(X_test[, j]), j] <- half_min
    } else {
      X_train[is.na(X_train[, j]), j] <- 0
      X_test[is.na(X_test[, j]), j] <- 0
    }
  }
  
  cat("Imputation complete\n\n")
  
  return(list(
    X_train = X_train,
    X_test = X_test,
    y_train = y_train,
    y_test = y_test
  ))
}

### ============================================================================
### MODULE 4: FEATURE SELECTION
### ============================================================================

perform_feature_selection <- function(X_train, X_test, y_train, classification_type) {
  cat("=== MODULE 4: FEATURE SELECTION ===\n")
  cat("Original number of features:", ncol(X_train), "\n")
  
  # Method 1: Variance-based filtering
  feature_var <- apply(X_train, 2, var)
  variance_threshold <- quantile(feature_var, 0.1)
  high_var_features <- feature_var > variance_threshold
  X_train_var <- X_train[, high_var_features, drop = FALSE]
  X_test_var <- X_test[, high_var_features, drop = FALSE]
  
  cat("After variance filtering:", ncol(X_train_var), "features\n")
  
  # Method 2: Correlation-based filtering
  if(ncol(X_train_var) > 1) {
    cor_matrix <- cor(X_train_var)
    high_cor_pairs <- which(abs(cor_matrix) > 0.95 & cor_matrix != 1, arr.ind = TRUE)
    features_to_remove <- unique(high_cor_pairs[, 2])
    if(length(features_to_remove) > 0) {
      X_train_cor <- X_train_var[, -features_to_remove, drop = FALSE]
      X_test_cor <- X_test_var[, -features_to_remove, drop = FALSE]
    } else {
      X_train_cor <- X_train_var
      X_test_cor <- X_test_var
    }
  } else {
    X_train_cor <- X_train_var
    X_test_cor <- X_test_var
  }
  
  cat("After correlation filtering:", ncol(X_train_cor), "features\n")
  
  # Method 3: Random Forest importance
  max_features_rf <- min(15, max(10, floor(nrow(X_train) / 3)), ncol(X_train_cor))
  
  if(ncol(X_train_cor) > 1 && max_features_rf > 0) {
    rf_data <- data.frame(y = y_train, X_train_cor)
    rf_model <- randomForest(y ~ ., data = rf_data, ntree = 100, importance = TRUE)
    rf_importance <- randomForest::importance(rf_model)[, "MeanDecreaseGini"]
    top_features_rf <- names(sort(rf_importance, decreasing = TRUE)[1:min(max_features_rf, length(rf_importance))])
    
    X_train_final <- X_train_cor[, top_features_rf, drop = FALSE]
    X_test_final <- X_test_cor[, top_features_rf, drop = FALSE]
  } else {
    X_train_final <- X_train_cor
    X_test_final <- X_test_cor
  }
  
  cat("After Random Forest selection:", ncol(X_train_final), "features\n")
  
  # Univariate AUC filtering
  univariate_auc <- sapply(1:ncol(X_train_final), function(i) {
    tryCatch({
      roc_temp <- roc(y_train, X_train_final[, i], quiet = TRUE)
      auc(roc_temp)
    }, error = function(e) 0)
  })
  
  if(sum(univariate_auc > 0.55) >= 5) {
    good_features <- univariate_auc > 0.55
    X_train_final <- X_train_final[, good_features, drop = FALSE]
    X_test_final <- X_test_final[, good_features, drop = FALSE]
    cat("After univariate AUC filtering:", ncol(X_train_final), "features\n")
  }
  
  cat("Final feature count:", ncol(X_train_final), "\n\n")
  
  # Store fixed feature set for LOPO-CV
  fixed_feature_set <- colnames(X_train_final)
  
  return(list(
    X_train_final = X_train_final,
    X_test_final = X_test_final,
    fixed_feature_set = fixed_feature_set
  ))
}

### ============================================================================
### MODULE 5: CLASS BALANCE HANDLING
### ============================================================================

handle_class_balance <- function(y_train) {
  cat("=== MODULE 5: CLASS BALANCE HANDLING ===\n")
  
  train_tbl <- table(y_train)
  cat("Class counts (original, no upsampling):\n")
  print(train_tbl)
  cat("\n")
  
  # Calculate class weights
  if(diff(range(train_tbl)) > 0) {
    class_weights <- table(y_train)
    class_weights <- max(class_weights) / class_weights
    cat("Class weights for balanced learning:\n")
    print(class_weights)
    cat("\n")
    cat("NOTE: Using class weights instead of upsampling to prevent overfitting\n\n")
  } else {
    class_weights <- NULL
    cat("Classes already balanced\n\n")
  }
  
  return(class_weights)
}

### ============================================================================
### MODULE 6: FEATURE SCALING
### ============================================================================

scale_features <- function(X_train, X_test) {
  cat("=== MODULE 6: FEATURE SCALING ===\n")
  
  scaler_center <- apply(X_train, 2, mean, na.rm = TRUE)
  scaler_scale <- apply(X_train, 2, sd, na.rm = TRUE)
  scaler_scale[scaler_scale == 0] <- 1
  
  X_train_scaled <- scale(X_train, center = scaler_center, scale = scaler_scale)
  X_test_scaled <- scale(X_test, center = scaler_center, scale = scaler_scale)
  
  cat("Features scaled\n\n")
  
  return(list(
    X_train_scaled = X_train_scaled,
    X_test_scaled = X_test_scaled,
    scaler_center = scaler_center,
    scaler_scale = scaler_scale
  ))
}

### ============================================================================
### MODULE 7: HYPERPARAMETER TUNING (PATIENT-GROUPED CV)
### ============================================================================

tune_hyperparameters <- function(X_train_final, y_train, train_data, class_weights) {
  cat("=== MODULE 7: HYPERPARAMETER TUNING ===\n")
  cat("Using patient-grouped CV for hyperparameter tuning...\n\n")
  
  # Prepare original data (non-upsampled) for CV
  scaler_center_orig <- apply(X_train_final, 2, mean, na.rm = TRUE)
  scaler_scale_orig <- apply(X_train_final, 2, sd, na.rm = TRUE)
  scaler_scale_orig[scaler_scale_orig == 0] <- 1
  rf_data_scaled_orig <- data.frame(
    y = y_train,
    scale(X_train_final, center = scaler_center_orig, scale = scaler_scale_orig),
    check.names = FALSE
  )
  
  # Patient-grouped CV
  patient_groups <- train_data %>%
    group_by(patient_id) %>%
    summarise(class = first(class), .groups = "drop")
  
  set.seed(123)
  cv_folds <- createFolds(patient_groups$class, k = 5, list = TRUE)
  
  # Parameter grid
  param_grid <- expand.grid(
    ntree = c(100, 200, 300),
    mtry = c(sqrt(ncol(X_train_final)), log2(ncol(X_train_final)), ncol(X_train_final)/3),
    nodesize = c(15, 20, 25, 30),
    maxnodes = c(5, 8, 10)
  )
  
  # Sample subset for efficiency
  set.seed(123)
  param_subset <- param_grid[sample(nrow(param_grid), min(30, nrow(param_grid))), ]
  
  cat("Testing", nrow(param_subset), "parameter combinations...\n\n")
  
  cv_results <- data.frame()
  
  for(i in 1:nrow(param_subset)) {
    params <- list(
      ntree = param_subset$ntree[i],
      mtry = param_subset$mtry[i],
      nodesize = param_subset$nodesize[i],
      maxnodes = param_subset$maxnodes[i]
    )
    
    fold_aucs <- numeric(5)
    
    for(fold in 1:5) {
      test_patients_fold <- patient_groups$patient_id[cv_folds[[fold]]]
      train_row_indices <- which(!train_data$patient_id %in% test_patients_fold)
      test_row_indices <- which(train_data$patient_id %in% test_patients_fold)
      
      train_fold_data <- rf_data_scaled_orig[train_row_indices, ]
      test_fold_data <- rf_data_scaled_orig[test_row_indices, ]
      
      if(length(unique(train_fold_data$y)) < 2 || length(unique(test_fold_data$y)) < 1) {
        fold_aucs[fold] <- NA
        next
      }
      
      rf_fold <- randomForest(
        y ~ .,
        data = train_fold_data,
        ntree = params$ntree,
        mtry = params$mtry,
        nodesize = params$nodesize,
        maxnodes = params$maxnodes,
        importance = TRUE
      )
      
      test_preds_fold <- predict(rf_fold, test_fold_data, type = "prob")[, 2]
      fold_roc <- roc(test_fold_data$y, test_preds_fold, quiet = TRUE)
      fold_aucs[fold] <- auc(fold_roc)
    }
    
    patient_cv_auc <- mean(fold_aucs, na.rm = TRUE)
    
    rf_full <- randomForest(
      y ~ .,
      data = rf_data_scaled_orig,
      ntree = params$ntree,
      mtry = params$mtry,
      nodesize = params$nodesize,
      maxnodes = params$maxnodes,
      importance = TRUE
    )
    
    train_preds_full <- predict(rf_full, rf_data_scaled_orig, type = "prob")[, 2]
    train_roc_full <- roc(rf_data_scaled_orig$y, train_preds_full, quiet = TRUE)
    train_auc_full <- auc(train_roc_full)
    
    train_cv_gap <- train_auc_full - patient_cv_auc
    
    cv_results <- rbind(cv_results, data.frame(
      ntree = params$ntree,
      mtry = params$mtry,
      nodesize = params$nodesize,
      maxnodes = params$maxnodes,
      patient_cv_auc = patient_cv_auc,
      train_auc = train_auc_full,
      train_cv_gap = train_cv_gap
    ))
    
    if(i %% 5 == 0) {
      cat("  Processed", i, "/", nrow(param_subset), "combinations\n")
    }
  }
  
  # Select best parameters
  cv_results$composite_score <- cv_results$patient_cv_auc - (cv_results$train_cv_gap * 1.5)
  best_idx <- which.max(cv_results$composite_score)
  best_params <- cv_results[best_idx, ]
  
  cat("\nOptimal parameters selected:\n")
  cat("Patient-grouped CV-AUC:", round(best_params$patient_cv_auc, 4), "\n")
  cat("Train-CV gap:", round(best_params$train_cv_gap, 4), "\n")
  cat("ntree:", best_params$ntree, ", mtry:", round(best_params$mtry, 2),
      ",nodesize:", best_params$nodesize, ", maxnodes:", best_params$maxnodes, "\n\n")
  
  return(best_params)
}

### ============================================================================
### MODULE 8: MODEL TRAINING
### ============================================================================

train_final_model <- function(X_train_scaled, y_train, best_params, class_weights) {
  cat("=== MODULE 8: MODEL TRAINING ===\n")
  
  rf_data_scaled <- data.frame(y = y_train, X_train_scaled, check.names = FALSE)
  
  rf_model_final <- randomForest(
    y ~ .,
    data = rf_data_scaled,
    ntree = best_params$ntree,
    mtry = best_params$mtry,
    nodesize = best_params$nodesize,
    maxnodes = best_params$maxnodes,
    importance = TRUE,
    proximity = TRUE,
    keep.inbag = TRUE,
    classwt = if(!is.null(class_weights)) class_weights else NULL
  )
  
  cat("Final model trained\n\n")
  
  return(rf_model_final)
}

### ============================================================================
### MODULE 9: MODEL EVALUATION
### ============================================================================

evaluate_model <- function(rf_model_final, X_test_scaled, y_test, X_train_final, y_train, 
                           scaler_center, scaler_scale, classification_type, class_labels) {
  cat("=== MODULE 9: MODEL EVALUATION ===\n")
  
  # Test set predictions
  pred_probs <- predict(rf_model_final, X_test_scaled, type = "prob")[, 2]
  roc_obj <- roc(response = y_test, predictor = pred_probs, quiet = TRUE)
  auc_val <- auc(roc_obj)
  
  # Training set predictions (on original data)
  pred_probs_train_orig <- predict(rf_model_final, 
                                    scale(X_train_final, center = scaler_center, scale = scaler_scale), 
                                    type = "prob")[, 2]
  train_roc <- roc(y_train, pred_probs_train_orig, quiet = TRUE)
  train_auc <- auc(train_roc)
  
  # Calculate optimal threshold
  all_coords <- coords(roc_obj, x = "all", ret = c("threshold", "sensitivity", "specificity"))
  youden_j <- all_coords$sensitivity + all_coords$specificity - 1
  best_idx <- which.max(youden_j)
  best_thresh <- all_coords$threshold[best_idx]
  
  pred_class <- ifelse(pred_probs > best_thresh, levels(y_test)[2], levels(y_test)[1])
  pred_class <- factor(pred_class, levels = levels(y_test))
  
  cat("Test Set AUC:", round(auc_val, 4), "\n")
  cat("Train Set AUC:", round(train_auc, 4), "\n")
  cat("Train-Test Gap:", round(train_auc - auc_val, 4), "\n\n")
  
  return(list(
    pred_probs = pred_probs,
    pred_class = pred_class,
    roc_obj = roc_obj,
    auc_val = auc_val,
    train_auc = train_auc,
    best_thresh = best_thresh
  ))
}

### ============================================================================
### MODULE 10: LOPO-CV
### ============================================================================

perform_lopo_cv <- function(df_wide, fixed_feature_set, best_params, classification_type) {
  cat("=== MODULE 10: LEAVE-ONE-PATIENT-OUT CROSS-VALIDATION ===\n")
  
  unique_patients <- unique(df_wide$patient_id)
  n_patients <- length(unique_patients)
  
  cat("Performing LOPO-CV on", n_patients, "patients...\n")
  
  lopo_predictions <- list()
  lopo_true_labels <- list()
  lopo_patient_ids <- list()
  lopo_fold_counter <- 0
  
  for(i in seq_along(unique_patients)) {
    test_patient <- unique_patients[i]
    
    if(i %% 5 == 0) {
      cat("Processed", i, "/", n_patients, "patients\n")
    }
    
    tryCatch({
      train_fold <- df_wide %>% filter(patient_id != test_patient)
      test_fold <- df_wide %>% filter(patient_id == test_patient)
      
      # Determine class levels based on classification type
      if(classification_type %in% c("DF_DHF", "DF_DHF_panel18")) {
        class_levels_fold <- c("DF", "DHF")
        positive_class <- "DHF"  # Second class is positive
      } else if(classification_type %in% c("HP_SUB", "HP_SUB_panel21")) {
        class_levels_fold <- c("hospitalized", "subclinical")
        positive_class <- "subclinical"  # Second class is positive
      } else if(classification_type == "AP_RP") {
        class_levels_fold <- c("AP", "RP")
        positive_class <- "RP"  # Second class is positive
      } else {
        # Default: use unique classes in order
        class_levels_fold <- sort(unique(c(train_fold$class, test_fold$class)))
        positive_class <- class_levels_fold[2]  # Assume second class is positive
      }
      
      # Ensure factor levels are in correct order
      y_train_fold <- factor(train_fold$class, levels = class_levels_fold)
      y_test_fold <- factor(test_fold$class, levels = class_levels_fold)
      
      # Check that we have both classes
      if(length(unique(y_train_fold)) < 2 || length(unique(y_test_fold)) < 1) {
        next
      }
      
      exclude_cols_fold <- c("patient_id", "class")
      if("timepoint" %in% colnames(train_fold)) exclude_cols_fold <- c(exclude_cols_fold, "timepoint")
      if("Group" %in% colnames(train_fold)) exclude_cols_fold <- c(exclude_cols_fold, "Group")
      
      X_train_fold <- train_fold %>% select(-all_of(exclude_cols_fold)) %>% as.matrix()
      X_test_fold <- test_fold %>% select(-all_of(exclude_cols_fold)) %>% as.matrix()
      
      # Imputation
      for(j in 1:ncol(X_train_fold)) {
        feature_values <- X_train_fold[, j]
        non_zero_values <- feature_values[!is.na(feature_values) & feature_values > 0]
        if(length(non_zero_values) > 0) {
          half_min <- min(non_zero_values) / 2
          X_train_fold[is.na(X_train_fold[, j]), j] <- half_min
          X_test_fold[is.na(X_test_fold[, j]), j] <- half_min
        } else {
          X_train_fold[is.na(X_train_fold[, j]), j] <- 0
          X_test_fold[is.na(X_test_fold[, j]), j] <- 0
        }
      }
      
      # Use fixed feature set
      if(length(fixed_feature_set) > 0 && all(fixed_feature_set %in% colnames(X_train_fold))) {
        X_train_fold_final <- X_train_fold[, fixed_feature_set, drop = FALSE]
        X_test_fold_final <- X_test_fold[, fixed_feature_set, drop = FALSE]
      } else {
        X_train_fold_final <- X_train_fold
        X_test_fold_final <- X_test_fold
      }
      
      # Scale
      scaler_c <- apply(X_train_fold_final, 2, mean, na.rm = TRUE)
      scaler_s <- apply(X_train_fold_final, 2, sd, na.rm = TRUE)
      scaler_s[scaler_s == 0] <- 1
      X_train_fold_sc <- scale(X_train_fold_final, center = scaler_c, scale = scaler_s)
      X_test_fold_sc <- scale(X_test_fold_final, center = scaler_c, scale = scaler_s)
      
      # Class weights
      train_tbl_fold <- table(y_train_fold)
      if(diff(range(train_tbl_fold)) > 0) {
        class_weights_fold <- table(y_train_fold)
        class_weights_fold <- max(class_weights_fold) / class_weights_fold
      } else {
        class_weights_fold <- NULL
      }
      
      # Train model
      rf_data_fold <- data.frame(y = y_train_fold, X_train_fold_sc)
      mtry_fold <- min(best_params$mtry, ncol(X_train_fold_sc))
      
      rf_fold <- randomForest(
        y ~ .,
        data = rf_data_fold,
        ntree = best_params$ntree,
        mtry = mtry_fold,
        nodesize = best_params$nodesize,
        maxnodes = best_params$maxnodes,
        importance = TRUE,
        classwt = if(!is.null(class_weights_fold)) class_weights_fold else NULL
      )
      
      pred_prob_fold <- predict(rf_fold, X_test_fold_sc, type = "prob")
      if(is.null(pred_prob_fold) || ncol(pred_prob_fold) < 2) {
        next
      }
      
      # CRITICAL FIX: Use column name matching positive class, not just column 2
      # This ensures we get the probability of the correct class regardless of factor level order
      if(positive_class %in% colnames(pred_prob_fold)) {
        pred_prob_fold <- pred_prob_fold[, positive_class, drop = TRUE]
      } else {
        # Fallback: use column 2 if positive class name not found
        # But log a warning
        cat("    Warning: Positive class", positive_class, "not found in predictions. Using column 2.\n")
        pred_prob_fold <- pred_prob_fold[, 2]
      }
      
      # Validate predictions are in valid range [0, 1]
      if(any(is.na(pred_prob_fold)) || any(pred_prob_fold < 0) || any(pred_prob_fold > 1)) {
        cat("Warning: Invalid predictions for patient", test_patient, "- skipping\n")
        next
      }
      
      lopo_predictions[[i]] <- pred_prob_fold
      lopo_true_labels[[i]] <- as.character(y_test_fold)
      lopo_patient_ids[[i]] <- rep(test_patient, length(y_test_fold))
      lopo_fold_counter <- lopo_fold_counter + 1
      
    }, error = function(e) {
      cat("Error processing patient", test_patient, ":", e$message, "\n")
    })
  }
  
  # Aggregate results
  lopo_predictions <- lopo_predictions[!sapply(lopo_predictions, is.null)]
  lopo_true_labels <- lopo_true_labels[!sapply(lopo_true_labels, is.null)]
  lopo_patient_ids <- lopo_patient_ids[!sapply(lopo_patient_ids, is.null)]
  
  if(length(lopo_predictions) > 0) {
    all_preds <- unlist(lopo_predictions)
    all_labels <- unlist(lopo_true_labels)
    
    # Determine class levels based on classification type
    if(classification_type %in% c("DF_DHF", "DF_DHF_panel18")) {
      class_levels_agg <- c("DF", "DHF")
      positive_class_agg <- "DHF"
    } else if(classification_type %in% c("HP_SUB", "HP_SUB_panel21")) {
      class_levels_agg <- c("hospitalized", "subclinical")
      positive_class_agg <- "subclinical"
    } else if(classification_type == "AP_RP") {
      class_levels_agg <- c("AP", "RP")
      positive_class_agg <- "RP"
    } else {
      class_levels_agg <- sort(unique(all_labels))
      positive_class_agg <- class_levels_agg[2]
    }
    
    # Ensure labels are factors with correct levels
    all_labels <- factor(all_labels, levels = class_levels_agg)
    
    # Validate predictions
    if(any(is.na(all_preds)) || any(all_preds < 0) || any(all_preds > 1)) {
      cat("WARNING: Invalid predictions detected. Removing invalid values.\n")
      valid_idx <- !is.na(all_preds) & all_preds >= 0 & all_preds <= 1
      all_preds <- all_preds[valid_idx]
      all_labels <- all_labels[valid_idx]
    }
    
    if(length(unique(all_labels)) >= 2 && length(all_preds) > 0) {
      # Compute ROC curve with explicit direction
      # For binary classification, we want to predict the positive class (second level)
      tryCatch({
        lopo_roc <- roc(all_labels, all_preds, quiet = TRUE, 
                       levels = class_levels_agg,
                       direction = "<")  # Higher predictions = positive class
        lopo_auc <- auc(lopo_roc)
        
        # Validate AUC is reasonable (should be >= 0.5 for a valid classifier)
        if(lopo_auc < 0.5) {
          cat("WARNING: LOPO-CV AUC < 0.5 (", round(lopo_auc, 4), 
              "). Trying reversed direction.\n")
          # Try reversed direction
          lopo_roc <- roc(all_labels, all_preds, quiet = TRUE,
                         levels = class_levels_agg,
                         direction = ">")  # Lower predictions = positive class
          lopo_auc <- auc(lopo_roc)
          if(lopo_auc < 0.5) {
            cat("WARNING: LOPO-CV AUC still < 0.5 after reversal. Check predictions.\n")
          }
        }
        
        cat("LOPO-CV complete\n")
        cat("Successful folds:", lopo_fold_counter, "/", n_patients, "\n")
        cat("LOPO-CV AUC:", round(lopo_auc, 4), "\n")
        cat("Predictions range:", round(range(all_preds), 4), "\n")
        cat("Class distribution:", paste(table(all_labels), collapse = ", "), "\n\n")
      }, error = function(e) {
        cat("ERROR computing LOPO-CV ROC:", e$message, "\n")
        cat("Predictions range:", range(all_preds, na.rm = TRUE), "\n")
        cat("Unique labels:", paste(unique(all_labels), collapse = ", "), "\n")
        lopo_roc <<- NULL
        lopo_auc <<- NA
      })
    } else {
      lopo_roc <- NULL
      lopo_auc <- NA
      cat("LOPO-CV failed: Insufficient class diversity\n")
      cat("Unique labels:", length(unique(all_labels)), "\n")
      cat("Total predictions:", length(all_preds), "\n\n")
    }
  } else {
    lopo_roc <- NULL
    lopo_auc <- NA
    cat("LOPO-CV failed: No successful predictions\n\n")
  }
  
  return(list(
    lopo_roc = lopo_roc,
    lopo_auc = lopo_auc
  ))
}

### ============================================================================
### MODULE 11: FEATURE IMPORTANCE
### ============================================================================

extract_feature_importance <- function(rf_model_final) {
  cat("=== MODULE 11: FEATURE IMPORTANCE ===\n")
  
  importance_matrix <- randomForest::importance(rf_model_final)
  importance_df <- data.frame(
    Feature = rownames(importance_matrix),
    MeanDecreaseGini = importance_matrix[, "MeanDecreaseGini"],
    MeanDecreaseAccuracy = importance_matrix[, "MeanDecreaseAccuracy"]
  ) %>%
    arrange(desc(MeanDecreaseGini))
  
  cat("Top 20 features:\n")
  print(head(importance_df, 20))
  cat("\n")
  
  return(importance_df)
}

### ============================================================================
### MODULE 12: PLOTTING
### ============================================================================

create_feature_importance_plot <- function(importance_df, classification_type, output_prefix) {
  cat("=== CREATING SEPARATE FEATURE IMPORTANCE PLOT ===\n")
  
  if(nrow(importance_df) == 0) {
    cat("No features to plot\n\n")
    return()
  }
  
  # Get all features sorted by MeanDecreaseGini (descending)
  # Show all features if less than 20, otherwise top 20
  n_features <- min(20, nrow(importance_df))
  top_features <- head(importance_df, n_features)
  
  # Determine title based on classification type
  title_map <- list(
    "AP_RP" = paste0("Top ", n_features, " Feature Importance (Random Forest AP vs RP)"),
    "DF_DHF" = paste0("Top ", n_features, " Feature Importance (Random Forest DF vs DHF)"),
    "HP_SUB" = paste0("Top ", n_features, " Feature Importance (Random Forest Hospitalized vs Subclinical)"),
    "DF_DHF_panel18" = paste0("Top ", n_features, " Feature Importance (Random Forest DF vs DHF - 18-gene panel)"),
    "HP_SUB_panel21" = paste0("Top ", n_features, " Feature Importance (Random Forest Hospitalized vs Subclinical - 21-gene panel)")
  )
  
  plot_title <- if(classification_type %in% names(title_map)) {
    title_map[[classification_type]]
  } else {
    paste0("Top ", n_features, " Feature Importance (Random Forest ", classification_type, ")")
  }
  
  # Ensure features are ordered by MeanDecreaseGini (descending) for the plot
  top_features$Feature <- factor(top_features$Feature, levels = rev(top_features$Feature))
  
  # Create the plot with color gradient
  p_feature_imp <- ggplot(top_features, aes(x = Feature, y = MeanDecreaseGini, fill = MeanDecreaseGini)) +
    geom_col() +
    scale_fill_gradient(low = "yellow", high = "purple", 
                       name = "Importance\n(Mean Decrease Gini)",
                       guide = guide_colorbar(title.position = "top", 
                                             title.hjust = 0.5,
                                             barwidth = 0.5,
                                             barheight = 10)) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
    scale_x_discrete(expand = expansion(add = c(0.5, 0.5))) +
    coord_flip() +
    labs(
      title = NULL,
      subtitle = NULL,
      x = "Gene names",
      y = "Mean Decrease Gini"
    ) +
    theme_minimal(base_size = 40) +
    theme(
      legend.position = "right",
      legend.title = make_text_element(size = 40, face = "bold", family = FONT_FAMILY),
      legend.text = make_text_element(size = 40, family = FONT_FAMILY),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.margin = ggplot2::margin(20, 25, 20, 20),
      axis.text.y = make_text_element(size = 40, family = FONT_FAMILY),
      axis.text.x = make_text_element(size = 40, family = FONT_FAMILY),
      axis.title = make_text_element(size = 40, family = FONT_FAMILY)
    )
  
  print(p_feature_imp)
  
  # Save in multiple formats with polished_new suffix (png 600 dpi; also pdf, svg, emf) - larger dimensions for readability
  base_filename <- paste0("figures/feature_importance_", output_prefix, "_polished_new.png")
  ggsave(base_filename, plot = p_feature_imp, width = 20, height = 16, dpi = 600, units = "in")
  ggsave(gsub("\\.png$", ".pdf", base_filename), plot = p_feature_imp, width = 20, height = 16, units = "in")
  ggsave(gsub("\\.png$", ".svg", base_filename), plot = p_feature_imp, width = 20, height = 16, units = "in")
  if (has_devEMF) tryCatch(
    ggsave(gsub("\\.png$", ".emf", base_filename), plot = p_feature_imp, width = 20, height = 16, units = "in"),
    error = function(e) message("EMF export skipped: ", conditionMessage(e))
  )
  cat("Feature importance plot saved\n\n")
}

create_plots <- function(roc_obj, auc_val, train_auc, lopo_roc, lopo_auc, 
                        importance_df, pred_probs, y_test, best_thresh,
                        classification_type, class_labels, output_prefix,
                        rf_model_final, X_train_final, y_train, scaler_center, scaler_scale,
                        n_train, n_test, n_lopo) {
  cat("=== MODULE 12: CREATING PLOTS ===\n")
  
  # Build ROC summary table (Method, N, AUROC) for annotation on ROC plots
  roc_table_df <- data.frame(
    Method = c("Train", "Test", "LOPO-CV"),
    N = c(n_train, n_test, n_lopo),
    AUROC = c(round(train_auc, 3), round(auc_val, 3), if(!is.na(lopo_auc)) round(lopo_auc, 3) else "-")
  )
  # Method column text: red, blue, darkgreen; all other text black. All borders black; fontsize 40.
  core_fg_col <- matrix("black", nrow = 3L, ncol = 3L)
  core_fg_col[1:3, 1L] <- c("red", "blue", "darkgreen")
  roc_table_grob <- gridExtra::tableGrob(roc_table_df, rows = NULL,
    theme = gridExtra::ttheme_default(
      core = list(
        fg_params = list(fontsize = 40, col = core_fg_col),
        bg_params = list(fill = "white", col = "black")
      ),
      colhead = list(fg_params = list(fontsize = 40, fontface = "bold"), bg_params = list(fill = "white", col = "black"))
    ))
  # Force black grid lines for all cells (override any inherited col on rect grobs)
  for (i in seq_along(roc_table_grob$grobs)) {
    cell <- roc_table_grob$grobs[[i]]
    if (inherits(cell, "gTree") && length(cell$children) > 0L) {
      for (j in seq_along(cell$children)) {
        if (inherits(cell$children[[j]], "rect") && !is.null(cell$children[[j]]$gp)) {
          cell$children[[j]]$gp$col <- "black"
        }
      }
      roc_table_grob$grobs[[i]] <- cell
    }
  }
  
  # Calculate CI for test AUC
  ci_obj <- ci.auc(roc_obj)
  ci_lower <- ci_obj[1]
  ci_upper <- ci_obj[3]
  
  # Plot 1: Performance summary (4-panel) - font size 40 for labels and ticks
  png(paste0("figures/performance_summary_", output_prefix, ".png"), 
      width = 19, height = 13, units = "in", res = 600)
  par(mfrow = c(2, 2), cex.main = 3.3, cex.lab = 3.3, cex.axis = 3.3)
  
  # Test ROC (smoothed)
  roc_smooth <- smooth(roc_obj, method = "density")
  plot(roc_smooth, col = "blue", lwd = 2,
       main = paste0("Test ROC (AUC = ", round(auc_val, 3), ")"))
  abline(a = 0, b = 1, lty = 2, col = "gray")
  legend("bottomright",
         legend = c(paste0("AUC = ", round(auc_val, 3)),
                    paste0("95% CI: [", round(ci_lower, 3), ", ", round(ci_upper, 3), "]")),
         bty = "n", cex = 2.5)
  
  # LOPO-CV ROC (smoothed)
  if(!is.null(lopo_roc)) {
    lopo_roc_smooth <- smooth(lopo_roc, method = "density")
    plot(lopo_roc_smooth, col = "darkgreen", lwd = 2,
         main = paste0("LOPO-CV ROC (AUC = ", round(lopo_auc, 3), ")"))
    abline(a = 0, b = 1, lty = 2, col = "gray")
  } else {
    plot(0, 0, type = "n", main = "LOPO-CV ROC (Not Available)", xlab = "", ylab = "")
  }
  
  # Performance comparison
  plot(c(1, 2, 3), c(train_auc, auc_val, if(!is.na(lopo_auc)) lopo_auc else 0),
       type = "b", pch = 19, cex = 2, col = c("red", "blue", "darkgreen"),
       ylim = c(0.5, 1.0), xlim = c(0.5, 3.5),
       xaxt = "n", xlab = "", ylab = "AUC",
       main = "Model Performance Comparison")
  axis(1, at = c(1, 2, 3), labels = c("Train", "Test", "LOPO-CV"))
  abline(h = 0.5, lty = 2, col = "gray")
  grid()
  
  # Feature importance (top 20)
  if(nrow(importance_df) > 0) {
    top_features <- head(importance_df, 20)
    barplot(rev(top_features$MeanDecreaseGini), 
            names.arg = rev(top_features$Feature),
            horiz = TRUE, las = 1, cex.names = 2.7,
            main = "Top 20 Feature Importance",
            xlab = "Mean Decrease Gini")
  }
  
  par(mfrow = c(1, 1), cex.main = 1, cex.lab = 1, cex.axis = 1)
  dev.off()
  
  # Plot 2: Combined ROC curves (Train, Test, LOPO-CV)
  # Note: train_auc is already calculated, but we need train_roc for plotting
  X_train_scaled_for_plot <- scale(X_train_final, center = scaler_center, scale = scaler_scale)
  pred_probs_train_orig <- predict(rf_model_final, X_train_scaled_for_plot, type = "prob")[, 2]
  train_roc <- roc(y_train, pred_probs_train_orig, quiet = TRUE)
  
  train_roc_smooth <- smooth(train_roc, method = "density")
  roc_smooth_combined <- smooth(roc_obj, method = "density")
  
  roc_combined_data <- data.frame(
    FPR = 1 - roc_smooth_combined$specificities,
    TPR = roc_smooth_combined$sensitivities,
    Method = "Test"
  )
  
  if(!is.null(lopo_roc)) {
    lopo_roc_smooth_combined <- smooth(lopo_roc, method = "density")
    lopo_roc_data <- data.frame(
      FPR = 1 - lopo_roc_smooth_combined$specificities,
      TPR = lopo_roc_smooth_combined$sensitivities,
      Method = "LOPO-CV"
    )
    roc_combined_data <- rbind(roc_combined_data, lopo_roc_data)
  }
  
  train_roc_data <- data.frame(
    FPR = 1 - train_roc_smooth$specificities,
    TPR = train_roc_smooth$sensitivities,
    Method = "Train"
  )
  roc_combined_data <- rbind(roc_combined_data, train_roc_data)
  
  roc_combined_data$Method <- factor(roc_combined_data$Method, levels = c("Train", "Test", "LOPO-CV"))
  
  p_roc_combined <- ggplot(roc_combined_data, aes(x = FPR, y = TPR, color = Method)) +
    geom_line(linewidth = 1.5) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    annotation_custom(grob = roc_table_grob, xmin = 0.55, xmax = 0.98, ymin = 0.02, ymax = 0.35) +
    scale_color_manual(values = c("Train" = "red", "Test" = "blue", "LOPO-CV" = "darkgreen")) +
    labs(
      title = NULL,
      subtitle = NULL,
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)",
      color = "Method"
    ) +
    theme_minimal(base_size = 40) +
    theme(
      aspect.ratio = 1,
      legend.position = "none",
      axis.text = make_text_element(size = 40, family = FONT_FAMILY),
      axis.title = make_text_element(size = 40, family = FONT_FAMILY),
      plot.title = element_blank(),
      plot.subtitle = element_blank()
    ) +
    coord_fixed() +
    xlim(0, 1) +
    ylim(0, 1)
  
  print(p_roc_combined)
  ggsave(paste0("figures/roc_combined_train_test_lopo_random_forest_", output_prefix, "_polished_new.png"), 
         plot = p_roc_combined, width = 19, height = 19, dpi = 600, units = "in")
  ggsave(paste0("figures/roc_combined_train_test_lopo_random_forest_", output_prefix, "_polished_new.pdf"), 
         plot = p_roc_combined, width = 19, height = 19, units = "in")
  ggsave(paste0("figures/roc_combined_train_test_lopo_random_forest_", output_prefix, "_polished_new.svg"), 
         plot = p_roc_combined, width = 19, height = 19, units = "in")
  if (has_devEMF) tryCatch(
    ggsave(paste0("figures/roc_combined_train_test_lopo_random_forest_", output_prefix, "_polished_new.emf"), 
           plot = p_roc_combined, width = 19, height = 19, units = "in"),
    error = function(e) message("EMF export skipped: ", conditionMessage(e))
  )
  # Plot 3: Combined ROC curves WITHOUT smoothing
  # Use existing ROC objects
  roc_obj_hr <- roc_obj
  train_roc_hr <- train_roc
  lopo_roc_hr <- lopo_roc
  
  roc_combined_data_unsmooth <- data.frame(
    FPR = 1 - roc_obj_hr$specificities,
    TPR = roc_obj_hr$sensitivities,
    Method = "Test"
  )
  
  if(!is.null(lopo_roc_hr)) {
    lopo_roc_data_unsmooth <- data.frame(
      FPR = 1 - lopo_roc_hr$specificities,
      TPR = lopo_roc_hr$sensitivities,
      Method = "LOPO-CV"
    )
    roc_combined_data_unsmooth <- rbind(roc_combined_data_unsmooth, lopo_roc_data_unsmooth)
  }
  
  train_roc_data_unsmooth <- data.frame(
    FPR = 1 - train_roc_hr$specificities,
    TPR = train_roc_hr$sensitivities,
    Method = "Train"
  )
  roc_combined_data_unsmooth <- rbind(roc_combined_data_unsmooth, train_roc_data_unsmooth)
  
  roc_combined_data_unsmooth$Method <- factor(roc_combined_data_unsmooth$Method, 
                                               levels = c("Train", "Test", "LOPO-CV"))
  
  p_roc_combined_unsmooth <- ggplot(roc_combined_data_unsmooth, aes(x = FPR, y = TPR, color = Method)) +
    geom_path(linewidth = 1.5) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    annotation_custom(grob = roc_table_grob, xmin = 0.55, xmax = 0.98, ymin = 0.02, ymax = 0.35) +
    scale_color_manual(values = c("Train" = "red", "Test" = "blue", "LOPO-CV" = "darkgreen")) +
    labs(
      title = NULL,
      subtitle = NULL,
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)",
      color = "Method"
    ) +
    theme_minimal(base_size = 40) +
    theme(
      aspect.ratio = 1,
      legend.position = "none",
      axis.text = make_text_element(size = 40, family = FONT_FAMILY),
      axis.title = make_text_element(size = 40, family = FONT_FAMILY),
      plot.title = element_blank(),
      plot.subtitle = element_blank()
    ) +
    coord_fixed() +
    xlim(0, 1) +
    ylim(0, 1)
  
  print(p_roc_combined_unsmooth)
  ggsave(paste0("figures/roc_combined_unsmoothed_random_forest_", output_prefix, "_polished_new.png"), 
         plot = p_roc_combined_unsmooth, width = 19, height = 19, dpi = 600, units = "in")
  ggsave(paste0("figures/roc_combined_unsmoothed_random_forest_", output_prefix, "_polished_new.pdf"), 
         plot = p_roc_combined_unsmooth, width = 19, height = 19, units = "in")
  ggsave(paste0("figures/roc_combined_unsmoothed_random_forest_", output_prefix, "_polished_new.svg"), 
         plot = p_roc_combined_unsmooth, width = 19, height = 19, units = "in")
  if (has_devEMF) tryCatch(
    ggsave(paste0("figures/roc_combined_unsmoothed_random_forest_", output_prefix, "_polished_new.emf"), 
           plot = p_roc_combined_unsmooth, width = 19, height = 19, units = "in"),
    error = function(e) message("EMF export skipped: ", conditionMessage(e))
  )
  # Plot 3b: Polished version with confidence intervals and visual refinements
  # Calculate confidence intervals for each ROC curve
  ci_test <- tryCatch({
    ci.se(roc_obj_hr, specificities = seq(0, 1, 0.05), boot.n = 2000)
  }, error = function(e) NULL)
  
  ci_train <- tryCatch({
    ci.se(train_roc_hr, specificities = seq(0, 1, 0.05), boot.n = 2000)
  }, error = function(e) NULL)
  
  ci_lopo <- if(!is.null(lopo_roc_hr)) {
    tryCatch({
      ci.se(lopo_roc_hr, specificities = seq(0, 1, 0.05), boot.n = 2000)
    }, error = function(e) NULL)
  } else {
    NULL
  }
  
  # Prepare CI data for plotting
  ci_data_list <- list()
  
  if(!is.null(ci_test)) {
    ci_test_df <- data.frame(
      FPR = 1 - attr(ci_test, "specificities"),
      lower = ci_test[, 1],
      upper = ci_test[, 3],
      Method = "Test"
    )
    ci_data_list[["Test"]] <- ci_test_df
  }
  
  if(!is.null(ci_train)) {
    ci_train_df <- data.frame(
      FPR = 1 - attr(ci_train, "specificities"),
      lower = ci_train[, 1],
      upper = ci_train[, 3],
      Method = "Train"
    )
    ci_data_list[["Train"]] <- ci_train_df
  }
  
  if(!is.null(ci_lopo)) {
    ci_lopo_df <- data.frame(
      FPR = 1 - attr(ci_lopo, "specificities"),
      lower = ci_lopo[, 1],
      upper = ci_lopo[, 3],
      Method = "LOPO-CV"
    )
    ci_data_list[["LOPO-CV"]] <- ci_lopo_df
  }
  
  ci_data_all <- if(length(ci_data_list) > 0) {
    do.call(rbind, ci_data_list)
  } else {
    NULL
  }
  
  # Create polished plot
  p_roc_polished <- ggplot(roc_combined_data_unsmooth, aes(x = FPR, y = TPR, color = Method)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
                color = "gray50", linewidth = 0.5) +
    {if(!is.null(ci_data_all)) {
      geom_ribbon(data = ci_data_all, aes(x = FPR, ymin = lower, ymax = upper, fill = Method),
                  alpha = 0.15, inherit.aes = FALSE, show.legend = FALSE)
    }} +
    geom_path(linewidth = 1.2) +
    annotation_custom(grob = roc_table_grob, xmin = 0.55, xmax = 0.98, ymin = 0.02, ymax = 0.35) +
    scale_color_manual(values = c("Train" = "red", "Test" = "blue", "LOPO-CV" = "darkgreen")) +
    scale_fill_manual(
      values = c("Train" = "red", "Test" = "blue", "LOPO-CV" = "darkgreen"),
      guide = "none"
    ) +
    labs(
      title = NULL,
      subtitle = NULL,
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)",
      color = "Method"
    ) +
    coord_equal() +
    theme_minimal(base_size = 40) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      legend.position = "none",
      axis.text = make_text_element(size = 40, family = FONT_FAMILY),
      axis.title = make_text_element(size = 40, family = FONT_FAMILY),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3)
    ) +
    scale_x_continuous(breaks = seq(0, 1, 0.2), expand = c(0.01, 0.01)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), expand = c(0.01, 0.01))
  
  print(p_roc_polished)
  ggsave(paste0("figures/roc_combined_polished_random_forest_", output_prefix, "_polished_new.png"), 
         plot = p_roc_polished, width = 19, height = 19, dpi = 600, units = "in")
  ggsave(paste0("figures/roc_combined_polished_random_forest_", output_prefix, "_polished_new.pdf"), 
         plot = p_roc_polished, width = 19, height = 19, units = "in")
  ggsave(paste0("figures/roc_combined_polished_random_forest_", output_prefix, "_polished_new.svg"), 
         plot = p_roc_polished, width = 19, height = 19, units = "in")
  if (has_devEMF) tryCatch(
    ggsave(paste0("figures/roc_combined_polished_random_forest_", output_prefix, "_polished_new.emf"), 
           plot = p_roc_polished, width = 19, height = 19, units = "in"),
    error = function(e) message("EMF export skipped: ", conditionMessage(e))
  )
  # Plot 4: Combined ROC curves WITHOUT smoothing WITH shaded AUC areas
  # Create polygons for shading under each curve using high-resolution ROC objects
  # Test curve
  test_polygon <- data.frame(
    FPR = c(0, 1 - roc_obj_hr$specificities, 1),
    TPR = c(0, roc_obj_hr$sensitivities, 0),
    Method = "Test"
  )
  
  # Train curve
  train_polygon <- data.frame(
    FPR = c(0, 1 - train_roc_hr$specificities, 1),
    TPR = c(0, train_roc_hr$sensitivities, 0),
    Method = "Train"
  )
  
  # LOPO-CV curve
  if(!is.null(lopo_roc_hr)) {
    lopo_polygon <- data.frame(
      FPR = c(0, 1 - lopo_roc_hr$specificities, 1),
      TPR = c(0, lopo_roc_hr$sensitivities, 0),
      Method = "LOPO-CV"
    )
    all_polygons <- rbind(test_polygon, train_polygon, lopo_polygon)
  } else {
    all_polygons <- rbind(test_polygon, train_polygon)
  }
  
  all_polygons$Method <- factor(all_polygons$Method, levels = c("Train", "Test", "LOPO-CV"))
  
  p_roc_combined_shaded <- ggplot(roc_combined_data_unsmooth, aes(x = FPR, y = TPR, color = Method)) +
    geom_polygon(data = all_polygons, aes(x = FPR, y = TPR, fill = Method), 
                 alpha = 0.2, inherit.aes = FALSE) +
    geom_path(linewidth = 1.5) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    annotation_custom(grob = roc_table_grob, xmin = 0.55, xmax = 0.98, ymin = 0.02, ymax = 0.35) +
    scale_color_manual(values = c("Train" = "red", "Test" = "blue", "LOPO-CV" = "darkgreen")) +
    scale_fill_manual(
      values = c("Train" = "red", "Test" = "blue", "LOPO-CV" = "darkgreen"),
      guide = "none"
    ) +
    labs(
      title = NULL,
      subtitle = NULL,
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)",
      color = "Method"
    ) +
    theme_minimal(base_size = 40) +
    theme(
      aspect.ratio = 1,
      legend.position = "none",
      axis.text = make_text_element(size = 40, family = FONT_FAMILY),
      axis.title = make_text_element(size = 40, family = FONT_FAMILY),
      plot.title = element_blank(),
      plot.subtitle = element_blank()
    ) +
    coord_fixed() +
    xlim(0, 1) +
    ylim(0, 1)
  
  print(p_roc_combined_shaded)
  ggsave(paste0("figures/roc_combined_shaded_auc_random_forest_", output_prefix, "_polished_new.png"), 
         plot = p_roc_combined_shaded, width = 19, height = 19, dpi = 600, units = "in")
  ggsave(paste0("figures/roc_combined_shaded_auc_random_forest_", output_prefix, "_polished_new.pdf"), 
         plot = p_roc_combined_shaded, width = 19, height = 19, units = "in")
  ggsave(paste0("figures/roc_combined_shaded_auc_random_forest_", output_prefix, "_polished_new.svg"), 
         plot = p_roc_combined_shaded, width = 19, height = 19, units = "in")
  if (has_devEMF) tryCatch(
    ggsave(paste0("figures/roc_combined_shaded_auc_random_forest_", output_prefix, "_polished_new.emf"), 
           plot = p_roc_combined_shaded, width = 19, height = 19, units = "in"),
    error = function(e) message("EMF export skipped: ", conditionMessage(e))
  )
  cat("All plots saved\n\n")
}

### ============================================================================
### MAIN EXECUTION
### ============================================================================

cat("============================================================================\n")
cat("UNIFIED DENGUE PATIENT CLASSIFIER - RANDOM FOREST\n")
cat("============================================================================\n\n")

# Step 1: Load and transform data
data_result <- load_and_transform_data(INPUT_FILE, CLASSIFICATION_TYPE)
df_wide <- data_result$df_wide
class_levels <- data_result$class_levels
class_labels <- data_result$class_labels

# Step 2: Check leakage and split
split_result <- check_leakage_and_split(df_wide, CLASSIFICATION_TYPE)
train_data <- split_result$train_data
test_data <- split_result$test_data

# Step 3: Preprocess
exclude_cols <- c("patient_id", "class")
if("timepoint" %in% colnames(train_data)) exclude_cols <- c(exclude_cols, "timepoint")
if("Group" %in% colnames(train_data)) exclude_cols <- c(exclude_cols, "Group")

preprocess_result <- preprocess_data(train_data, test_data, exclude_cols)
X_train <- preprocess_result$X_train
X_test <- preprocess_result$X_test
y_train <- preprocess_result$y_train
y_test <- preprocess_result$y_test

# Step 4: Feature selection
feature_result <- perform_feature_selection(X_train, X_test, y_train, CLASSIFICATION_TYPE)
X_train_final <- feature_result$X_train_final
X_test_final <- feature_result$X_test_final
fixed_feature_set <- feature_result$fixed_feature_set

# Step 5: Class balance
class_weights <- handle_class_balance(y_train)

# Step 6: Feature scaling
scale_result <- scale_features(X_train_final, X_test_final)
X_train_scaled <- scale_result$X_train_scaled
X_test_scaled <- scale_result$X_test_scaled
scaler_center <- scale_result$scaler_center
scaler_scale <- scale_result$scaler_scale

# Step 7: Hyperparameter tuning
best_params <- tune_hyperparameters(X_train_final, y_train, train_data, class_weights)

# Step 8: Train final model
rf_model_final <- train_final_model(X_train_scaled, y_train, best_params, class_weights)

# Step 9: Evaluate model
eval_result <- evaluate_model(rf_model_final, X_test_scaled, y_test, X_train_final, y_train,
                              scaler_center, scaler_scale, CLASSIFICATION_TYPE, class_labels)

# Step 10: LOPO-CV
lopo_result <- perform_lopo_cv(df_wide, fixed_feature_set, best_params, CLASSIFICATION_TYPE)

# Step 11: Feature importance
importance_df <- extract_feature_importance(rf_model_final)

# Step 12: Create plots
n_train_plot <- length(unique(train_data$patient_id))
n_test_plot <- length(unique(test_data$patient_id))
n_lopo_plot <- length(unique(df_wide$patient_id))
create_plots(eval_result$roc_obj, eval_result$auc_val, eval_result$train_auc,
            lopo_result$lopo_roc, lopo_result$lopo_auc, importance_df,
            eval_result$pred_probs, y_test, eval_result$best_thresh,
            CLASSIFICATION_TYPE, class_labels, OUTPUT_PREFIX,
            rf_model_final, X_train_final, y_train, scaler_center, scaler_scale,
            n_train_plot, n_test_plot, n_lopo_plot)

# Step 13: Create separate feature importance plot
create_feature_importance_plot(importance_df, CLASSIFICATION_TYPE, OUTPUT_PREFIX)

# Final summary
cat("============================================================================\n")
cat("FINAL SUMMARY\n")
cat("============================================================================\n")
cat("Classification:", CLASSIFICATION_TYPE, "\n")
cat("Train AUC:", round(eval_result$train_auc, 4), "\n")
cat("Test AUC:", round(eval_result$auc_val, 4), "\n")
if(!is.na(lopo_result$lopo_auc)) {
  cat("LOPO-CV AUC:", round(lopo_result$lopo_auc, 4), "\n")
}
cat("Features:", ncol(X_train_final), "\n")
cat("\nAll outputs saved with prefix:", OUTPUT_PREFIX, "\n")
cat("============================================================================\n")

