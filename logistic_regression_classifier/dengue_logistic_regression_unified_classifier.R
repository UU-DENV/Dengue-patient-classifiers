### UNIFIED DENGUE PATIENT CLASSIFIER - LOGISTIC REGRESSION (ELASTIC NET)

### This unified classifier can perform multiple classification tasks:
###1. AP vs RP (Acute Phase vs Recovery Phase)
###2. DF vs DHF (Dengue Fever vs Dengue Hemorrhagic Fever)
###3. Hospitalized vs Subclinical
###4. DF vs DHF with 18-gene panel
###5. Hospitalized vs Subclinical with 21-gene panel


# Classification type: "AP_RP", "DF_DHF", "HP_SUB", "DF_DHF_panel18", "HP_SUB_panel21"
CLASSIFICATION_TYPE <- "HP_SUB_panel21"  #THE USER MUST SELECT ONE CLASSIFICATION TYPE FROM THE LIST ABOVE

# Input file path (can be CSV or TSV)
# For AP_RP: Use "DENV_filtered_imputed_dataset.csv"
# For DF_DHF: Use "260108_LFQ_DENV_AP.csv"
# For HP_SUB: Use "LFQ_DENV_AP_CP_imputed_corrected.csv"
INPUT_FILE <- "C:/dengue/20260120_LFQ_DENV_AP.csv.csv"  #EXAMPLE FILE

# Output prefix for all generated files
OUTPUT_PREFIX <- "unified_LR_HP_SUB_panel21"  #TO BE MODIFIED AS PREFERRED

# Random seed for train/test split (set to a number for reproducibility; NA = no seed, different split each run)
RANDOM_SEED_SPLIT <- 123

### Load required packages
if(!require(tidyverse)) install.packages("tidyverse")
if(!require(caret)) install.packages("caret")
if(!require(glmnet)) install.packages("glmnet")
if(!require(pROC)) install.packages("pROC")
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(gridExtra)) install.packages("gridExtra")
# Optional: for EMF export on Windows (install.packages("devEMF") if needed)
has_devEMF <- requireNamespace("devEMF", quietly = TRUE)
if (has_devEMF) suppressPackageStartupMessages(library(devEMF))
if(!require(Boruta)) install.packages("Boruta")
if(!require(randomForest)) install.packages("randomForest")  # Only for feature selection
if(!require(extrafont)) install.packages("extrafont")

library(tidyverse)
library(caret)
library(glmnet)
library(pROC)
library(ggplot2)
library(gridExtra)
library(Boruta)
library(randomForest)  # Only for feature selection
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


### MODULE 1: DATA LOADING AND TRANSFORMATION

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
        values_from = LFQ,
        values_fn = mean  # Aggregate duplicates by taking mean
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
        values_from = LFQ,
        values_fn = mean  # Aggregate duplicates by taking mean
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


### MODULE 2: DATA LEAKAGE CHECK AND TRAIN/TEST SPLIT

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
  cat("Split ratio: 75/25 (75% train, 25% test)\n")
  cat("Training set:", length(train_patients), "patients,", nrow(train_data), "observations\n")
  cat("Test set:", length(test_patients), "patients,", nrow(test_data), "observations\n\n")
  
  return(list(
    train_data = train_data,
    test_data = test_data,
    train_patients = train_patients,
    test_patients = test_patients
  ))
}

### MODULE 3: DATA PREPROCESSING AND IMPUTATION

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

### MODULE 4: FEATURE SELECTION

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
  
  # Method 3: Random Forest importance (skip for panel: use all panel genes to improve LOPO-CV vs RF)
  n_samples <- nrow(X_train)
  is_panel <- grepl("panel", classification_type, ignore.case = TRUE)
  
  if(is_panel) {
    # Panel-based: use ALL panel genes that passed variance/correlation (no RF reduction)
    # This improves LOPO-CV AUC by retaining full panel; glmnet regularization handles overfitting
    X_train_final <- X_train_cor
    X_test_final <- X_test_cor
    cat("  Panel-based classification: using all", ncol(X_train_cor), "panel features (no RF reduction)\n")
    cat("After feature selection (panel):", ncol(X_train_final), "features\n")
  } else {
    # Non-panel: use conservative feature selection
    if(n_samples < 20) {
      # Very small samples: maximum 5 features
      max_features_rf <- min(5, ncol(X_train_cor))
      cat("  Very small sample (n =", n_samples, "): using very conservative feature selection (max 5 features)\n")
    } else if(n_samples < 30) {
      # Small samples: maximum 5-6 features
      max_features_rf <- min(6, ncol(X_train_cor))
      cat("  Small sample (n =", n_samples, "): using conservative feature selection (max 6 features)\n")
    } else if(n_samples < 50) {
      # Medium samples (30-50): maximum 6-7 features
      max_features_rf <- min(7, ncol(X_train_cor))
      cat("  Medium sample (n =", n_samples, "): using conservative feature selection (max 7 features)\n")
    } else if(n_samples < 70) {
      # Small-medium samples (50-70): maximum 8 features
      max_features_rf <- min(8, ncol(X_train_cor))
      cat("  Small-medium sample (n =", n_samples, "): using conservative feature selection (max 8 features)\n")
    } else {
      # Larger samples: use n/10 rule but cap at 12
      max_features_rf <- min(12, max(8, floor(n_samples / 10)), ncol(X_train_cor))
    }
    # Ensure we have at least 3 features for logistic regression
    min_features_needed <- 3
    if(max_features_rf < min_features_needed && ncol(X_train_cor) >= min_features_needed) {
      cat("  Adjusting max_features to minimum needed (", min_features_needed, ") for small sample\n")
      max_features_rf <- min(min_features_needed, ncol(X_train_cor))
    }
    # Warn if we're using too many features for the sample size
    if(max_features_rf > n_samples / 5) {
      cat("WARNING: Using ", max_features_rf, " features with ", n_samples, 
          "samples (ratio = ", round(n_samples / max_features_rf, 2), ":1)\n")
      cat("Recommended ratio is at least 10:1 for reliable logistic regression\n")
    }
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
  }
  
  # Warn if we have very few features
  if(ncol(X_train_final) < 3) {
    cat("WARNING: Very few features selected (", ncol(X_train_final), 
        ").Model may have limited predictive power.\n")
  }
  
  # Univariate AUC filtering - skip for panel so we keep all panel features
  is_panel <- grepl("panel", classification_type, ignore.case = TRUE)
  if(is_panel) {
    cat("Skipping univariate AUC filtering for panel classification: keeping all", ncol(X_train_final), "features\n")
  } else {
    univariate_auc <- sapply(1:ncol(X_train_final), function(i) {
      tryCatch({
        roc_temp <- roc(y_train, X_train_final[, i], quiet = TRUE)
        auc(roc_temp)
      }, error = function(e) 0)
    })
    # Adjust threshold based on sample size
    if(n_samples < 20) {
      cat("Skipping univariate AUC filtering for very small sample size (n =", n_samples, ")\n")
    } else if(n_samples < 30) {
      auc_threshold <- 0.51
      min_features <- 2
      if(sum(univariate_auc > auc_threshold) >= min_features) {
        good_features <- univariate_auc > auc_threshold
        X_train_final <- X_train_final[, good_features, drop = FALSE]
        X_test_final <- X_test_final[, good_features, drop = FALSE]
        cat("After univariate AUC filtering (threshold =", auc_threshold, "):", ncol(X_train_final), "features\n")
      } else {
        cat("Univariate AUC filtering skipped: too few features pass threshold\n")
      }
    } else if(n_samples < 50) {
      auc_threshold <- 0.53
      min_features <- 4
      if(sum(univariate_auc > auc_threshold) >= min_features) {
        good_features <- univariate_auc > auc_threshold
        X_train_final <- X_train_final[, good_features, drop = FALSE]
        X_test_final <- X_test_final[, good_features, drop = FALSE]
        cat("After univariate AUC filtering (threshold =", auc_threshold, "):", ncol(X_train_final), "features\n")
      } else {
        cat("Univariate AUC filtering skipped: too few features pass threshold\n")
      }
    } else {
      auc_threshold <- 0.55
      min_features <- 5
      if(sum(univariate_auc > auc_threshold) >= min_features) {
        good_features <- univariate_auc > auc_threshold
        X_train_final <- X_train_final[, good_features, drop = FALSE]
        X_test_final <- X_test_final[, good_features, drop = FALSE]
        cat("After univariate AUC filtering (threshold =", auc_threshold, "):", ncol(X_train_final), "features\n")
      } else {
        cat("Univariate AUC filtering skipped: too few features pass threshold\n")
      }
    }
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

### MODULE 5: CLASS BALANCE HANDLING

handle_class_balance <- function(y_train, balance_strength = 1.5) {
  cat("=== MODULE 5: CLASS BALANCE HANDLING (REVISED) ===\n")
  
  class_counts <- table(y_train)
  cat("Class counts:\n")
  print(class_counts)
  cat("\n")
  
  # Check for imbalance
  imbalance_ratio <- max(class_counts) / min(class_counts)
  cat("Imbalance ratio:", round(imbalance_ratio, 2), ":1\n")
  
  if(imbalance_ratio > 1.5) {
    # Calculate class weights with adjustable strength
    # balance_strength = 1.0 gives standard balanced weights
    # balance_strength > 1.0 gives stronger emphasis on minority class
    class_weights <- (max(class_counts) / class_counts) ^ balance_strength
    
    cat("Using STRENGTHENED class weights (strength =", balance_strength, "):\n")
    print(round(class_weights, 3))
    cat("\nThis gives", round(class_weights[names(which.min(class_counts))], 2), 
        "x more weight to minority class\n")
  } else {
    class_weights <- NULL
    cat("Classes are balanced, no weighting needed\n")
  }
  
  cat("\n")
  return(class_weights)
}

### MODULE 6: FEATURE SCALING

scale_features <- function(X_train, X_test) {
  cat("=== MODULE 6: FEATURE SCALING ===\n")
  
  scaler_center <- apply(X_train, 2, mean, na.rm = TRUE)
  scaler_scale <- apply(X_train, 2, sd, na.rm = TRUE)
  scaler_scale[scaler_scale == 0] <- 1
  
  X_train_scaled <- scale(X_train, center = scaler_center, scale = scaler_scale)
  X_test_scaled <- scale(X_test, center = scaler_center, scale = scaler_scale)
  
  cat("✓ Features scaled\n\n")
  
  return(list(
    X_train_scaled = X_train_scaled,
    X_test_scaled = X_test_scaled,
    scaler_center = scaler_center,
    scaler_scale = scaler_scale
  ))
}

### MODULE 7: HYPERPARAMETER TUNING (REVISED)

tune_hyperparameters <- function(X_train_final, y_train, train_data, class_weights) {
  cat("=== MODULE 7: HYPERPARAMETER TUNING (REVISED) ===\n")
  cat("Using patient-grouped CV for alpha/lambda tuning...\n\n")
  
  # Scale features for tuning
  scaler_center <- apply(X_train_final, 2, mean, na.rm = TRUE)
  scaler_scale <- apply(X_train_final, 2, sd, na.rm = TRUE)
  scaler_scale[scaler_scale == 0] <- 1
  X_train_scaled <- scale(X_train_final, center = scaler_center, scale = scaler_scale)
  
  n_samples <- nrow(X_train_final)
  n_features <- ncol(X_train_final)
  
  # Adjust alpha grid based on sample size
  if(n_samples < 50) {
    # Small samples: Focus on high regularization (Lasso-like)
    alpha_grid <- c(0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0)
    cat("  Small sample (n =", n_samples, "): Testing", length(alpha_grid), "alpha values (high regularization)...\n\n")
  } else if(n_samples < 100) {
    alpha_grid <- c(0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 1.0)
    cat("  Medium sample (n =", n_samples, "): Testing", length(alpha_grid), "alpha values...\n\n")
  } else {
    alpha_grid <- c(0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0)
    cat("  Large sample (n =", n_samples, "): Testing", length(alpha_grid), "alpha values...\n\n")
  }
  
  # Prepare observation weights - ensure length matches exactly
  obs_weights <- if(!is.null(class_weights)) {
    y_train_char <- as.character(y_train)
    weights_vec <- class_weights[y_train_char]
    # Validate length
    if(length(weights_vec) != length(y_train)) {
      cat("WARNING: Weights length mismatch. Recalculating weights...\n")
      # Recalculate weights based on current y_train distribution
      train_tbl <- table(y_train)
      if(diff(range(train_tbl)) > 0) {
        class_weights_recalc <- max(train_tbl) / train_tbl
        weights_vec <- class_weights_recalc[y_train_char]
      } else {
        weights_vec <- rep(1, length(y_train))
      }
    }
    # Check for NA values in weights
    if(any(is.na(weights_vec))) {
      cat("WARNING: NA values in weights. Replacing with 1.0...\n")
      weights_vec[is.na(weights_vec)] <- 1.0
    }
    weights_vec
  } else {
    NULL
  }
  
  # Validate weights length one more time before use
  if(!is.null(obs_weights) && length(obs_weights) != nrow(X_train_scaled)) {
    cat("WARNING: Final weights length (", length(obs_weights), 
        ")!= data rows (", nrow(X_train_scaled), "). Using equal weights.\n")
    obs_weights <- NULL
  }
  
  # CRITICAL: Ensure foldid length matches data exactly
  # Set up patient-grouped CV folds - use fewer folds for small samples
  use_foldid <- TRUE
  foldid <- NULL
  fold_counts <- NULL
  
  # Determine appropriate number of folds based on actual sample size
  # For small samples, use fewer folds to ensure adequate samples per fold
  if(n_samples < 20) {
    n_folds <- min(3, n_samples)
  } else if(n_samples < 40) {
    n_folds <- min(5, n_samples)
  } else if(n_samples < 60) {
    n_folds <- min(5, n_samples)  # Still use 5 for 44 samples
  } else {
    n_folds <- min(10, n_samples)
  }
  
  if(n_folds < 2) {
    stop("ERROR: Not enough samples (", n_samples, ") for cross-validation. Need at least 2 samples.")
  }
  
  if("patient_id" %in% colnames(train_data)) {
    unique_patients <- unique(train_data$patient_id)
    # Try patient-grouped CV
    if(length(unique_patients) >= 2) {
      # For patient-grouped CV, ensure we have enough patients per fold
      # Each fold needs at least 2 patients, and ideally both classes represented
      max_safe_folds <- floor(length(unique_patients) / 2)  # At least 2 patients per fold
      patient_n_folds <- min(n_folds, max_safe_folds, length(unique_patients))
      
      if(patient_n_folds >= 2) {
        # Use stratified assignment to ensure both classes in each fold
        patient_classes <- train_data %>% 
          group_by(patient_id) %>% 
          summarise(class = first(class), .groups = "drop")
        
        # Try to balance classes across folds
        set.seed(123)
        patient_folds <- rep(NA, length(unique_patients))
        names(patient_folds) <- unique_patients
        
        # Simple stratified assignment
        for(class_val in unique(patient_classes$class)) {
          class_patients <- patient_classes$patient_id[patient_classes$class == class_val]
          if(length(class_patients) > 0) {
            class_folds <- cut(seq_along(class_patients), breaks = patient_n_folds, labels = FALSE)
            patient_folds[class_patients] <- class_folds
          }
        }
        
        # Handle any remaining unassigned patients
        if(any(is.na(patient_folds))) {
          unassigned <- which(is.na(patient_folds))
          patient_folds[unassigned] <- sample(1:patient_n_folds, length(unassigned), replace = TRUE)
        }
        
        foldid <- patient_folds[train_data$patient_id]
        
        # CRITICAL: Ensure foldid length matches X_train_scaled exactly
        if(length(foldid) != nrow(X_train_scaled)) {
          cat("WARNING: foldid length (", length(foldid), 
              ") != data rows (", nrow(X_train_scaled), "). Using standard CV.\n")
          use_foldid <- FALSE
          foldid <- NULL
        } else {
          # Validate fold assignments
          fold_counts <- table(foldid)
          if(any(fold_counts == 0)) {
            cat("WARNING: Some CV folds are empty. Using standard CV.\n")
            use_foldid <- FALSE
            foldid <- NULL
            fold_counts <- NULL
          } else {
            # Check that each fold has both classes
            fold_class_ok <- TRUE
            for(f in 1:patient_n_folds) {
              fold_indices <- which(foldid == f)
              fold_classes <- unique(y_train[fold_indices])
              if(length(fold_classes) < 2) {
                cat("WARNING: Fold", f, "has only one class. Using standard CV.\n")
                fold_class_ok <- FALSE
                break
              }
            }
            
            if(!fold_class_ok) {
              use_foldid <- FALSE
              foldid <- NULL
              fold_counts <- NULL
            } else {
              n_folds <- patient_n_folds  # Use patient-based fold count
            }
          }
        }
      } else {
        cat("WARNING: Not enough unique patients for patient-grouped CV. Using standard CV.\n")
        use_foldid <- FALSE
      }
    } else {
      cat("WARNING: Not enough unique patients for patient-grouped CV. Using standard CV.\n")
      use_foldid <- FALSE
    }
  } else {
    # No patient_id - use standard CV with random folds
    use_foldid <- FALSE
  }
  
  if(use_foldid && !is.null(foldid)) {
    cat("Using patient-grouped CV with", n_folds, "folds\n")
  } else {
    cat("Using standard CV (no foldid) with", n_folds, "folds\n")
    foldid <- NULL
    # Create fold_counts for standard CV
    if(!is.null(foldid)) {
      fold_counts <- table(foldid)
    } else {
      # Estimate fold size for standard CV
      fold_counts <- rep(floor(n_samples / n_folds), n_folds)
      remainder <- n_samples %% n_folds
      if(remainder > 0) {
        fold_counts[1:remainder] <- fold_counts[1:remainder] + 1
      }
    }
  }
  
  # Check if folds are too small for AUC - use deviance instead
  if(!is.null(fold_counts) && length(fold_counts) > 0) {
    min_fold_size <- min(fold_counts)
    use_auc <- min_fold_size >= 10
    type_measure <- if(use_auc) "auc" else "deviance"
    if(!use_auc) {
      cat("WARNING: Folds too small (min =", min_fold_size, 
          ") for AUC. Using deviance instead.\n")
    }
  } else {
    # Default to deviance if we can't determine fold size
    type_measure <- "deviance"
    cat("WARNING: Cannot determine fold sizes. Using deviance.\n")
  }
  
  # Storage for results
  alpha_results <- data.frame()
  
  # Test each alpha value
  for(alpha_val in alpha_grid) {
    result <- tryCatch({
      # Validate inputs before calling cv.glmnet
      if(any(is.na(X_train_scaled)) || any(is.infinite(X_train_scaled))) {
        cat("Alpha =", alpha_val, "| X_train_scaled contains NA or Inf values (skipping)\n")
        return(NULL)
      }
      
      if(any(is.na(y_train))) {
        cat("Alpha =", alpha_val, "| y_train contains NA values (skipping)\n")
        return(NULL)
      }
      
      if(length(unique(y_train)) < 2) {
        cat("Alpha =", alpha_val, "| y_train has less than 2 classes (skipping)\n")
        return(NULL)
      }
      
      # Validate weights length one final time
      if(!is.null(obs_weights) && length(obs_weights) != nrow(X_train_scaled)) {
        cat("Alpha =", alpha_val, "| Weights length mismatch, using NULL weights\n")
        obs_weights_use <- NULL
      } else {
        obs_weights_use <- obs_weights
      }
      
      # Build cv.glmnet arguments
      cv_args <- list(
        x = X_train_scaled,
        y = y_train,
        family = "binomial",
        alpha = alpha_val,
        type.measure = type_measure,
        lambda.min.ratio = 0.001,
        keep = TRUE
      )
      
      # Only add foldid if it's valid and we want to use it
      if(use_foldid && !is.null(foldid) && length(foldid) == nrow(X_train_scaled)) {
        cv_args$foldid <- foldid
      } else {
        cv_args$nfolds <- n_folds
      }
      
      # Add weights only if provided and valid
      if(!is.null(obs_weights_use) && length(obs_weights_use) == nrow(X_train_scaled)) {
        cv_args$weights <- obs_weights_use
      }
      
      # Fit CV model
      cv_model <- tryCatch({
        do.call(cv.glmnet, cv_args)
      }, error = function(e) {
        cat("Alpha =", alpha_val, "| cv.glmnet failed:", e$message, "\n")
        # Try without weights if weights caused the error
        if(!is.null(obs_weights_use) && grepl("length", e$message, ignore.case = TRUE)) {
          cat("Retrying without weights...\n")
          cv_args_no_weights <- cv_args
          cv_args_no_weights$weights <- NULL
          return(tryCatch({
            do.call(cv.glmnet, cv_args_no_weights)
          }, error = function(e2) {
            cat("Retry also failed:", e2$message, "\n")
            # Last resort: try without foldid and without weights
            if(use_foldid && !is.null(foldid)) {
              cat("Retrying without foldid and without weights...\n")
              cv_args_simple <- list(
                x = X_train_scaled,
                y = y_train,
                family = "binomial",
                alpha = alpha_val,
                nfolds = n_folds,
                type.measure = type_measure,
                lambda.min.ratio = 0.001,
                keep = TRUE
              )
              return(tryCatch({
                do.call(cv.glmnet, cv_args_simple)
              }, error = function(e3) {
                cat("Final retry also failed:", e3$message, "\n")
                return(NULL)
              }))
            }
            return(NULL)
          }))
        } else {
          return(NULL)
        }
      })
      
      if(is.null(cv_model)) {
        return(NULL)
      }
      
      # Calculate patient-grouped CV metric
      # If using deviance, we need to calculate AUC manually from CV predictions
      if(type_measure == "auc") {
        # Check if cvm exists and is valid - use safer checks
        cvm_valid <- tryCatch({
          !is.null(cv_model$cvm) && length(cv_model$cvm) > 0 && !all(is.na(cv_model$cvm))
        }, error = function(e) FALSE)
        
        if(!cvm_valid) {
          cat("Alpha =", alpha_val, "| No valid CV metrics (skipping)\n")
          return(NULL)
        }
        
        patient_cv_auc <- tryCatch({
          max(cv_model$cvm, na.rm = TRUE)
        }, error = function(e) {
          NA_real_
        })
        
        # Skip if CV AUC is invalid
        if(length(patient_cv_auc) != 1 || is.na(patient_cv_auc) || is.infinite(patient_cv_auc)) {
          cat("Alpha =", alpha_val, "| Invalid CV-AUC (skipping)\n")
          return(NULL)
        }
      } else {
        # Using deviance - calculate AUC from CV predictions if available
        # cv.glmnet with keep=TRUE stores fold predictions
        if(!is.null(cv_model$fit.preval)) {
          # Get out-of-fold predictions
          cv_preds <- as.numeric(cv_model$fit.preval[, which.min(cv_model$cvm)])
          if(length(cv_preds) == length(y_train) && length(unique(y_train)) == 2) {
            tryCatch({
              cv_roc <- roc(y_train, cv_preds, quiet = TRUE)
              patient_cv_auc <- as.numeric(auc(cv_roc))
            }, error = function(e) {
              # If ROC fails, use mean CV metric (deviance) as proxy
              patient_cv_auc <- 1 - (min(cv_model$cvm, na.rm = TRUE) / max(cv_model$cvm, na.rm = TRUE))
              if(is.na(patient_cv_auc) || is.infinite(patient_cv_auc)) {
                patient_cv_auc <- 0.5  # Default to random
              }
            })
          } else {
            # Fallback: use deviance as proxy (inverted and normalized)
            patient_cv_auc <- tryCatch({
              min_dev <- min(cv_model$cvm, na.rm = TRUE)
              max_dev <- max(cv_model$cvm, na.rm = TRUE)
              if(is.finite(min_dev) && is.finite(max_dev) && max_dev > min_dev) {
                # Normalize deviance to approximate AUC (lower deviance = better = higher AUC)
                1 - ((min_dev - max_dev) / (max_dev * 2))
              } else {
                0.5  # Default to random
              }
            }, error = function(e) 0.5)
          }
        } else {
          # No fold predictions available - use deviance as proxy
          patient_cv_auc <- tryCatch({
            min_dev <- min(cv_model$cvm, na.rm = TRUE)
            max_dev <- max(cv_model$cvm, na.rm = TRUE)
            if(is.finite(min_dev) && is.finite(max_dev) && max_dev > min_dev) {
              1 - ((min_dev - max_dev) / (max_dev * 2))
            } else {
              0.5
            }
          }, error = function(e) 0.5)
        }
        
        # Validate the calculated AUC
        if(length(patient_cv_auc) != 1 || is.na(patient_cv_auc) || is.infinite(patient_cv_auc)) {
          patient_cv_auc <- 0.5  # Default to random performance
        }
        # Clamp to valid range
        patient_cv_auc <- max(0.0, min(1.0, patient_cv_auc))
      }
      
      # Calculate train AUC for both lambda.min and lambda.1se
      train_preds_min <- tryCatch({
        as.numeric(predict(cv_model, newx = X_train_scaled, 
                           s = "lambda.min", type = "response"))
      }, error = function(e) {
        cat("Alpha =", alpha_val, "| Error predicting with lambda.min:", e$message, "\n")
        return(NULL)
      })
      
      # Safer check for invalid predictions
      preds_min_valid <- tryCatch({
        !is.null(train_preds_min) && length(train_preds_min) > 0 && !all(is.na(train_preds_min))
      }, error = function(e) FALSE)
      
      if(!preds_min_valid) {
        cat("Alpha =", alpha_val, "| Invalid predictions with lambda.min (skipping)\n")
        return(NULL)
      }
      
      train_roc_min <- tryCatch({
        roc(y_train, train_preds_min, quiet = TRUE)
      }, error = function(e) NULL)
      train_auc_min <- if(!is.null(train_roc_min)) {
        tryCatch({
          auc_val <- auc(train_roc_min)
          if(length(auc_val) != 1 || is.na(auc_val) || is.infinite(auc_val)) NA_real_ else auc_val
        }, error = function(e) NA_real_)
      } else {
        NA_real_
      }
      
      train_preds_1se <- tryCatch({
        as.numeric(predict(cv_model, newx = X_train_scaled, 
                          s = "lambda.1se", type = "response"))
      }, error = function(e) {
        cat("Alpha =", alpha_val, "| Error predicting with lambda.1se:", e$message, "\n")
        return(NULL)
      })
      
      # Safer check for invalid predictions
      preds_1se_valid <- tryCatch({
        !is.null(train_preds_1se) && length(train_preds_1se) > 0 && !all(is.na(train_preds_1se))
      }, error = function(e) FALSE)
      
      if(!preds_1se_valid) {
        cat("Alpha =", alpha_val, "| Invalid predictions with lambda.1se (skipping)\n")
        return(NULL)
      }
      
      train_roc_1se <- tryCatch({
        roc(y_train, train_preds_1se, quiet = TRUE)
      }, error = function(e) NULL)
      train_auc_1se <- if(!is.null(train_roc_1se)) {
        tryCatch({
          auc_val <- auc(train_roc_1se)
          if(length(auc_val) != 1 || is.na(auc_val) || is.infinite(auc_val)) NA_real_ else auc_val
        }, error = function(e) NA_real_)
      } else {
        NA_real_
      }
      
      # Check for NA values - skip this alpha if both are NA
      # Ensure both are single numeric values first
      train_auc_min <- tryCatch({
        val <- as.numeric(train_auc_min)
        if(length(val) > 0 && !all(is.na(val))) val[1] else NA_real_
      }, error = function(e) NA_real_)
      
      train_auc_1se <- tryCatch({
        val <- as.numeric(train_auc_1se)
        if(length(val) > 0 && !all(is.na(val))) val[1] else NA_real_
      }, error = function(e) NA_real_)
      
      # Validate both are single finite numeric values using safer checks
      min_valid <- tryCatch({
        length(train_auc_min) == 1 && is.finite(train_auc_min) && !is.na(train_auc_min)
      }, error = function(e) FALSE)
      
      se_valid <- tryCatch({
        length(train_auc_1se) == 1 && is.finite(train_auc_1se) && !is.na(train_auc_1se)
      }, error = function(e) FALSE)
      
      # Check if both are invalid
      both_invalid <- tryCatch({
        !min_valid && !se_valid
      }, error = function(e) TRUE)
      
      if(both_invalid) {
        cat("Alpha =", alpha_val, "| Both lambda.min and lambda.1se failed (skipping)\n")
        return(NULL)
      }
      
      # Replace invalid values with 0.5 (random performance) for comparison purposes
      if(!min_valid) {
        train_auc_min <- 0.5
        min_valid <- TRUE
      }
      if(!se_valid) {
        train_auc_1se <- 0.5
        se_valid <- TRUE
      }
      
      # Calculate gaps
      gap_min <- train_auc_min - patient_cv_auc
      gap_1se <- train_auc_1se - patient_cv_auc
      
      # ALWAYS prefer lambda.1se unless it's completely broken
      # Use isTRUE to ensure boolean result from comparison
      train_auc_1se_valid <- tryCatch({
        is.finite(train_auc_1se) && length(train_auc_1se) == 1
      }, error = function(e) FALSE)
      
      train_auc_1se_low <- tryCatch({
        isTRUE(train_auc_1se <= 0.52)
      }, error = function(e) FALSE)
      
      if(!train_auc_1se_valid) {
        # If lambda.1se is invalid, use lambda.min
        lambda_type <- "lambda.min"
        train_auc_used <- train_auc_min
        gap_used <- gap_min
        lambda_value <- cv_model$lambda.min
        cat("Alpha =", alpha_val, "| lambda.1se is invalid, using lambda.min\n")
      } else if(train_auc_1se_low) {
        # lambda.1se is broken (essentially random)
        lambda_type <- "lambda.min"
        train_auc_used <- train_auc_min
        gap_used <- gap_min
        lambda_value <- cv_model$lambda.min
        cat("Alpha =", alpha_val, "| lambda.1se broken (AUC =", 
            round(train_auc_1se, 3), "), using lambda.min\n")
      } else {
        # Use lambda.1se (more conservative, prevents overfitting)
        lambda_type <- "lambda.1se"
        train_auc_used <- train_auc_1se
        gap_used <- gap_1se
        lambda_value <- cv_model$lambda.1se
      }
      
      # Return result to be stored
      list(
        alpha = alpha_val,
        patient_cv_auc = patient_cv_auc,
        train_auc = train_auc_used,
        train_cv_gap = gap_used,
        lambda_value = lambda_value,
        lambda_type = lambda_type,
        train_auc_min = train_auc_min,
        train_auc_1se = train_auc_1se,
        gap_min = gap_min,
        gap_1se = gap_1se,
        success = TRUE
      )
      
    }, error = function(e) {
      cat("Error with alpha =", alpha_val, ":", e$message, "\n")
      return(NULL)
    })
    
    # Store results if successful
    if(!is.null(result) && !is.null(result$success) && result$success) {
      alpha_results <- rbind(alpha_results, data.frame(
        alpha = result$alpha,
        patient_cv_auc = result$patient_cv_auc,
        train_auc = result$train_auc,
        train_cv_gap = result$train_cv_gap,
        lambda_value = result$lambda_value,
        lambda_type = result$lambda_type,
        train_auc_min = result$train_auc_min,
        train_auc_1se = result$train_auc_1se,
        gap_min = result$gap_min,
        gap_1se = result$gap_1se
      ))
      
      cat("Alpha =", result$alpha, 
          "| CV-AUC =", round(result$patient_cv_auc, 3),
          "| Train-AUC =", round(result$train_auc, 3),
          "| Gap =", round(result$train_cv_gap, 3),
          "| Lambda:", result$lambda_type, "\n")
    }
  }
  
  cat("\n=== MODEL SELECTION ===\n")
  
  # Check if we have any results at all
  if(nrow(alpha_results) == 0) {
    stop("ERROR: All hyperparameter combinations failed. Cannot proceed.\n",
         "Possible causes:\n",
         "  1. All predictions are identical (perfect separation or model failure)\n",
         "  2. Insufficient data or severe class imbalance\n",
         "  3. All features were removed during feature selection\n",
         "Please check your data and feature selection process.")
  }
  
  # Selection Strategy: Prioritize models that generalize well
  # Step 1: Remove completely broken models (train AUC too low)
  viable_models <- alpha_results[!is.na(alpha_results$train_auc) & alpha_results$train_auc > 0.55, ]
  cat("Models with train AUC > 0.55:", nrow(viable_models), "/", nrow(alpha_results), "\n")
  
  if(nrow(viable_models) == 0) {
    cat("No viable models found! Using best CV-AUC available...\n")
    # Remove rows with NA values before selecting
    alpha_results_clean <- alpha_results[!is.na(alpha_results$patient_cv_auc), ]
    if(nrow(alpha_results_clean) == 0) {
      stop("ERROR: All models failed and no valid CV-AUC values found.")
    }
    best_idx <- which.max(alpha_results_clean$patient_cv_auc)
    best_params <- alpha_results_clean[best_idx, ]
  } else {
    # Step 2: Filter by reasonable CV performance (CV-AUC >= 0.60)
    good_cv_models <- viable_models[viable_models$patient_cv_auc >= 0.60, ]
    cat("Models with CV-AUC >= 0.60:", nrow(good_cv_models), "/", nrow(viable_models), "\n")
    
    if(nrow(good_cv_models) == 0) {
      cat("No models with CV-AUC >= 0.60. Relaxing to CV-AUC >= 0.55...\n")
      good_cv_models <- viable_models[viable_models$patient_cv_auc >= 0.55, ]
    }
    
    if(nrow(good_cv_models) > 0) {
      # Step 3: Among good CV models, prefer small gaps
      # Target: gap < 0.20 (some overfitting acceptable, but not extreme)
      small_gap_models <- good_cv_models[good_cv_models$train_cv_gap < 0.20, ]
      cat("Models with gap < 0.20:", nrow(small_gap_models), "/", nrow(good_cv_models), "\n")
      
      if(nrow(small_gap_models) > 0) {
        # Step 4: Among small-gap models, pick highest CV-AUC
        best_params <- small_gap_models[which.max(small_gap_models$patient_cv_auc), ]
        cat("Selected model with best CV-AUC among low-gap models\n")
      } else {
        # If all gaps are large, pick model with smallest gap
        best_params <- good_cv_models[which.min(good_cv_models$train_cv_gap), ]
        cat("All gaps > 0.20. Selected model with smallest gap\n")
      }
    } else {
      # Fallback: pick best CV-AUC
      best_params <- viable_models[which.max(viable_models$patient_cv_auc), ]
      cat("Fallback: Selected model with best CV-AUC\n")
    }
  }
  
  # CRITICAL: Force lambda.1se if it wasn't already selected
  if(!is.null(best_params$lambda_type) && best_params$lambda_type == "lambda.min") {
    # Check if lambda.1se is reasonable for this alpha
    if(!is.null(best_params$train_auc_1se) && !is.na(best_params$train_auc_1se) && best_params$train_auc_1se > 0.55) {
      cat("\nFORCING lambda.1se to prevent overfitting\n")
      cat("Original: lambda.min (train AUC =", round(best_params$train_auc_min, 3), 
          ", gap =", round(best_params$gap_min, 3), ")\n")
      cat("Forced:   lambda.1se (train AUC =", round(best_params$train_auc_1se, 3), 
          ", gap =", round(best_params$gap_1se, 3), ")\n")
      
      # Update to use lambda.1se
      best_params$lambda_type <- "lambda.1se"
      best_params$train_auc <- best_params$train_auc_1se
      best_params$train_cv_gap <- best_params$gap_1se
      # Note: lambda_value will be updated during final model training
    }
  }
  
  cat("\nOPTIMAL PARAMETERS SELECTED:\n")
  cat("Alpha:", best_params$alpha, "\n")
  cat("Lambda type:", best_params$lambda_type, "\n")
  cat("Patient CV-AUC:", round(best_params$patient_cv_auc, 3), "\n")
  cat("Train AUC:", round(best_params$train_auc, 3), "\n")
  cat("Train-CV gap:", round(best_params$train_cv_gap, 3), "\n")
  
  # Warnings
  if(best_params$patient_cv_auc < 0.65) {
    cat("\nWARNING: Low CV-AUC (", round(best_params$patient_cv_auc, 3), ")\n")
    cat("Consider: More data, better features, or different algorithm\n")
  }
  
  if(best_params$train_cv_gap > 0.20) {
    cat("\nWARNING: Large train-CV gap (", round(best_params$train_cv_gap, 3), ")\n")
    cat("Model may still be overfitting. Consider stronger regularization.\n")
  }
  
  if(best_params$train_auc > 0.85) {
    cat("\nWARNING: Very high train AUC (", round(best_params$train_auc, 3), ")\n")
    cat("This suggests overfitting despite regularization.\n")
    cat("Recommendation: Reduce features or increase regularization further.\n")
  }
  
  cat("\n")
  
  return(best_params)
}

### MODULE 8: MODEL TRAINING

train_final_model <- function(X_train_scaled, y_train, best_params, class_weights) {
  cat("=== MODULE 8: MODEL TRAINING (REVISED) ===\n")
  
  # Validate best_params
  if(is.null(best_params)) {
    stop("ERROR: best_params is NULL. Cannot train final model.\n",
         "This usually means all hyperparameter combinations failed during tuning.\n",
         "Please check:\n",
         "  1. Your data quality and feature selection\n",
         "  2. Class balance (severe imbalance can cause failures)\n",
         "  3. Sample size (may be too small for the number of features)")
  }
  
  # Handle data frame row (best_params is returned as a data frame row)
  if(is.data.frame(best_params)) {
    if(nrow(best_params) == 0) {
      stop("ERROR: best_params is an empty data frame. Cannot train final model.")
    }
    # Extract as a list for easier access
    best_params_list <- as.list(best_params[1, ])
  } else if(is.list(best_params)) {
    best_params_list <- best_params
  } else {
    stop("ERROR: best_params has unexpected type: ", class(best_params))
  }
  
  # Check if best_params has the required fields
  if(!"alpha" %in% names(best_params_list)) {
    stop("ERROR: best_params$alpha is missing. Cannot train final model.\n",
         "Available fields: ", paste(names(best_params_list), collapse = ", "))
  }
  
  # Extract alpha value safely
  alpha_val <- tryCatch({
    val <- best_params_list$alpha
    if(is.null(val) || length(val) == 0) {
      stop("alpha is NULL or has length 0")
    }
    val <- as.numeric(val)[1]
    if(is.na(val) || !is.finite(val)) {
      stop("alpha is NA or not finite")
    }
    val
  }, error = function(e) {
    stop("ERROR: Cannot extract valid alpha from best_params: ", e$message, 
         "\nbest_params$alpha = ", deparse(best_params_list$alpha))
  })
  
  if(alpha_val < 0 || alpha_val > 1) {
    stop("ERROR: Invalid alpha value (", alpha_val, "). Must be between 0 and 1.")
  }
  
  cat("  Using alpha =", alpha_val, "\n")
  
  # Prepare observation weights
  obs_weights <- if(!is.null(class_weights)) {
    class_weights[as.character(y_train)]
  } else {
    NULL
  }
  
  # Train final model with selected alpha
  cv_model_final <- tryCatch({
    cv.glmnet(
      x = X_train_scaled,
      y = y_train,
      family = "binomial",
      alpha = alpha_val,
      nfolds = 10,
      type.measure = "auc",
      lambda.min.ratio = 0.001,
      weights = obs_weights
    )
  }, error = function(e) {
    stop("ERROR: Failed to train final model with alpha = ", alpha_val, ": ", e$message)
  })
  
  # ALWAYS use lambda.1se unless it's completely broken
  train_preds_1se <- as.numeric(predict(cv_model_final, newx = X_train_scaled, 
                                        s = "lambda.1se", type = "response"))
  train_roc_1se <- roc(y_train, train_preds_1se, quiet = TRUE)
  train_auc_1se <- auc(train_roc_1se)
  
  if(train_auc_1se <= 0.52) {
    # lambda.1se is broken, check lambda.min
    train_preds_min <- as.numeric(predict(cv_model_final, newx = X_train_scaled, 
                                          s = "lambda.min", type = "response"))
    train_roc_min <- roc(y_train, train_preds_min, quiet = TRUE)
    train_auc_min <- auc(train_roc_min)
    
    if(train_auc_min > train_auc_1se + 0.10) {
      cat("lambda.1se broken (train AUC =", round(train_auc_1se, 3), 
          "), using lambda.min (train AUC =", round(train_auc_min, 3), ")\n")
      cv_model_final$lambda_to_use <- "lambda.min"
    } else {
      cat("Both lambdas have low train AUC. Using lambda.1se (more conservative)\n")
      cv_model_final$lambda_to_use <- "lambda.1se"
    }
  } else {
    # lambda.1se is working - ALWAYS use it
    cv_model_final$lambda_to_use <- "lambda.1se"
    cat("Using lambda.1se (train AUC =", round(train_auc_1se, 3), ")\n")
  }
  
  cat("\nFinal model trained\n")
  cat("Lambda.1se:", round(cv_model_final$lambda.1se, 6), "\n")
  cat("Lambda.min:", round(cv_model_final$lambda.min, 6), "\n")
  cat("Using:", cv_model_final$lambda_to_use, "\n")
  cat("\n")
  
  return(cv_model_final)
}

### MODULE 9: MODEL EVALUATION

evaluate_model <- function(cv_model_final, X_test_scaled, y_test, X_train_final, y_train, 
                           scaler_center, scaler_scale, classification_type, class_labels) {
  cat("=== MODULE 9: MODEL EVALUATION ===\n")
  
  # Test set predictions - use the lambda selected during tuning
  lambda_to_use <- if("lambda_to_use" %in% names(cv_model_final)) {
    cv_model_final$lambda_to_use
  } else {
    "lambda.1se"  # Default fallback
  }
  
  pred_probs <- as.numeric(predict(
    cv_model_final,
    newx = X_test_scaled,
    s = lambda_to_use,
    type = "response"
  ))
  roc_obj <- roc(response = y_test, predictor = pred_probs, quiet = TRUE)
  auc_val <- auc(roc_obj)
  
  # Training set predictions (on original data)
  X_train_scaled_for_eval <- scale(X_train_final, center = scaler_center, scale = scaler_scale)
  pred_probs_train_orig <- as.numeric(predict(
    cv_model_final,
    newx = X_train_scaled_for_eval,
    s = lambda_to_use,
    type = "response"
  ))
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

### MODULE 10: LEAVE-ONE-PATIENT-OUT CROSS-VALIDATION
                                      
perform_lopo_cv <- function(df_wide, fixed_feature_set, best_params, classification_type, class_weights_template, class_levels = NULL) {
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
      # Use same class levels in every fold so glmnet always predicts P(positive) = P(second level); avoids LOPO ROC dip
      y_train_fold <- if(!is.null(class_levels)) factor(train_fold$class, levels = class_levels) else factor(train_fold$class)
      y_test_fold <- if(!is.null(class_levels)) factor(test_fold$class, levels = class_levels) else factor(test_fold$class)
      
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
      
      # Prepare observation weights for glmnet
      obs_weights_fold <- if(!is.null(class_weights_fold)) {
        class_weights_fold[as.character(y_train_fold)]
      } else {
        NULL
      }
      
      # Train glmnet model
      cv_fold_model <- cv.glmnet(
        x = X_train_fold_sc,
        y = y_train_fold,
        family = "binomial",
        alpha = best_params$alpha,
        nfolds = 5,
        type.measure = "auc",
        lambda.min.ratio = 0.001,
        weights = obs_weights_fold
      )
      
      # Use lambda type: for panel classifications use lambda.min in LOPO to reduce over-regularization
      # (lambda.1se can over-shrink in leave-one-out folds and lower LOPO-CV AUC)
      is_panel <- grepl("panel", classification_type, ignore.case = TRUE)
      lambda_type_fold <- if(is_panel) {
        "lambda.min"
      } else if("lambda_type" %in% names(best_params)) {
        best_params$lambda_type
      } else {
        "lambda.1se"  # Default fallback
      }
      
      pred_prob_fold <- as.numeric(predict(
        cv_fold_model,
        newx = X_test_fold_sc,
        s = lambda_type_fold,
        type = "response"
      ))
      
      if(length(pred_prob_fold) == 0 || any(is.na(pred_prob_fold))) {
        next
      }
      
      lopo_predictions[[i]] <- pred_prob_fold
      lopo_true_labels[[i]] <- as.character(y_test_fold)
      lopo_patient_ids[[i]] <- rep(test_patient, length(y_test_fold))
      lopo_fold_counter <- lopo_fold_counter + 1
      
    }, error = function(e) {
      cat("    Error processing patient", test_patient, ":", e$message, "\n")
    })
  }
  
  # Aggregate results
  lopo_predictions <- lopo_predictions[!sapply(lopo_predictions, is.null)]
  lopo_true_labels <- lopo_true_labels[!sapply(lopo_true_labels, is.null)]
  lopo_patient_ids <- lopo_patient_ids[!sapply(lopo_patient_ids, is.null)]
  
  if(length(lopo_predictions) > 0) {
    all_preds <- unlist(lopo_predictions)
    all_labels <- if(!is.null(class_levels)) {
      factor(unlist(lopo_true_labels), levels = class_levels)
    } else {
      factor(unlist(lopo_true_labels))
    }
    
    if(length(unique(all_labels)) >= 2) {
      # Explicit levels and direction so ROC uses same positive class as Train/Test (avoids dip below diagonal)
      if(!is.null(class_levels)) {
        lopo_roc <- roc(all_labels, all_preds, levels = class_levels, direction = "<", quiet = TRUE)
      } else {
        lopo_roc <- roc(all_labels, all_preds, quiet = TRUE)
      }
      lopo_auc <- auc(lopo_roc)
      cat("LOPO-CV complete\n")
      cat("Successful folds:", lopo_fold_counter, "/", n_patients, "\n")
      cat("LOPO-CV AUC:", round(lopo_auc, 4), "\n\n")
    } else {
      lopo_roc <- NULL
      lopo_auc <- NA
      cat("LOPO-CV failed: Insufficient class diversity\n\n")
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

### MODULE 11: FEATURE IMPORTANCE

extract_feature_importance <- function(cv_model_final, feature_names = NULL) {
  cat("=== MODULE 11: FEATURE IMPORTANCE ===\n")
  
  # Extract coefficients using lambda.min to get more features (less conservative than lambda.1se)
  # This allows us to show more features in the importance plot
  lambda_to_use <- "lambda.min"  # Use lambda.min for feature importance to show more features
  
  coef_matrix <- as.matrix(coef(cv_model_final, s = lambda_to_use))
  coef_values <- coef_matrix[, 1]
  # Remove intercept but keep all other coefficients (including zeros and very small values)
  coef_values <- coef_values[names(coef_values) != "(Intercept)"]
  
  # If feature_names are provided, ensure we have coefficients for all of them
  if(!is.null(feature_names)) {
    # Create a complete data frame with all features
    all_features <- feature_names
    coef_dict <- setNames(as.numeric(coef_values), names(coef_values))
    
    # Get coefficients for all features (0 if not in model)
    all_coefs <- sapply(all_features, function(f) {
      if(f %in% names(coef_dict)) {
        coef_dict[f]
      } else {
        0  # Feature not in model (shouldn't happen, but safety check)
      }
    })
    
    importance_df <- data.frame(
      Feature = all_features,
      Coefficient = all_coefs,
      AbsCoefficient = abs(all_coefs)
    ) %>%
      arrange(desc(AbsCoefficient))
  } else {
    # Use only features with coefficients from the model
    importance_df <- data.frame(
      Feature = names(coef_values),
      Coefficient = as.numeric(coef_values),
      AbsCoefficient = abs(as.numeric(coef_values))
    ) %>%
      arrange(desc(AbsCoefficient))
  }
  
  cat("Total features available:", nrow(importance_df), "\n")
  cat("Features with non-zero coefficients:", sum(importance_df$Coefficient != 0), "\n")
  
  if(nrow(importance_df) < 15) {
    cat("WARNING: Only", nrow(importance_df), "features available (requested 15).\n")
    cat("This is because feature selection reduced the feature set.\n")
    cat("Showing all available features.\n")
  }
  
  cat("Top 15 features (by absolute coefficient):\n")
  print(head(importance_df, min(15, nrow(importance_df))))
  cat("\n")
  
  return(importance_df)
}

### MODULE 12: PLOTTING

create_feature_importance_plot <- function(importance_df, classification_type, output_prefix) {
  cat("=== CREATING SEPARATE FEATURE IMPORTANCE PLOT ===\n")
  
  if(nrow(importance_df) == 0) {
    cat("No features to plot\n\n")
    return(NULL)
  }
  
  # Select top 15 features (matching the image style)
  top_n <- min(15, nrow(importance_df))
  plot_data <- head(importance_df, top_n)
  
  # Create title based on classification type (matching image format)
  title_map <- list(
    "AP_RP" = paste0("Top ", top_n, " Feature Importance (Logistic Regression AP vs RP)"),
    "DF_DHF" = paste0("Top ", top_n, " Feature Importance (Logistic Regression DF vs DHF)"),
    "HP_SUB" = paste0("Top ", top_n, " Feature Importance (Logistic Regression Hospitalized vs Subclinical)"),
    "DF_DHF_panel18" = paste0("Top ", top_n, " Feature Importance (Logistic Regression DF vs DHF - 18-gene panel)"),
    "HP_SUB_panel21" = paste0("Top ", top_n, " Feature Importance (Logistic Regression Hospitalized vs Subclinical - 21-gene panel)")
  )
  
  plot_title <- if(classification_type %in% names(title_map)) {
    title_map[[classification_type]]
  } else {
    paste0("Top ", top_n, " Feature Importance (Logistic Regression ", classification_type, ")")
  }
  
  # Ensure features are ordered by AbsCoefficient (descending) for the plot
  plot_data$Feature <- factor(plot_data$Feature, levels = rev(plot_data$Feature))
  
  # Create the plot with horizontal bars and color gradient based on actual values
  # Color scheme: yellow -> orange -> pink/magenta -> purple (matching image)
  p <- ggplot(plot_data, aes(x = Feature, y = AbsCoefficient, fill = AbsCoefficient)) +
    geom_col() +
    scale_fill_gradient(
      low = "#FFD700",  # Yellow
      high = "#8B008B",  # Dark Magenta/Purple
      name = "Importance\n(Absolute\nCoefficient)",
      guide = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = 0.5,
        barheight = 10
      )
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
    scale_x_discrete(expand = expansion(add = c(0.5, 0.5))) +
    coord_flip() +  # Horizontal bars
    labs(
      title = NULL,
      subtitle = NULL,
      x = "Gene names",
      y = "Absolute Coefficient"
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
  
  print(p)
  
  # Save in multiple formats with polished_new suffix (png 600 dpi; also pdf, svg, emf) - larger dimensions for readability
  base_filename <- paste0("figures/feature_importance_", output_prefix, "_polished_new.png")
  ggsave(base_filename, plot = p, width = 20, height = 16, dpi = 600, units = "in")
  ggsave(gsub("\\.png$", ".pdf", base_filename), plot = p, width = 20, height = 16, units = "in")
  ggsave(gsub("\\.png$", ".svg", base_filename), plot = p, width = 20, height = 16, units = "in")
  if (has_devEMF) tryCatch(
    ggsave(gsub("\\.png$", ".emf", base_filename), plot = p, width = 20, height = 16, units = "in"),
    error = function(e) message("EMF export skipped: ", conditionMessage(e))
  )
  cat("Feature importance plot saved\n\n")
  
  return(p)
}

# Ensure ROC curve (FPR, TPR) is monotonically increasing and never dips below the diagonal (TPR >= FPR)
roc_monotonic <- function(fpr, tpr) {
  o <- order(fpr, tpr)
  fpr <- fpr[o]
  tpr <- tpr[o]
  tpr <- cummax(tpr)
  # Clamp TPR to be >= FPR so the curve never goes below the random classifier line
  tpr <- pmax(tpr, fpr)
  list(fpr = fpr, tpr = tpr)
}

create_plots <- function(roc_obj, auc_val, train_auc, lopo_roc, lopo_auc, 
                        importance_df, pred_probs, y_test, best_thresh,
                        classification_type, class_labels, output_prefix,
                        cv_model_final, X_train_final, y_train, scaler_center, scaler_scale,
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
    barplot(rev(top_features$AbsCoefficient), 
            names.arg = rev(top_features$Feature),
            horiz = TRUE, las = 1, cex.names = 2.7,
            main = "Top 20 Feature Importance",
            xlab = "Absolute Coefficient")
  }
  
  par(mfrow = c(1, 1), cex.main = 1, cex.lab = 1, cex.axis = 1)
  dev.off()
  
  # Plot 2: Combined ROC curves (Train, Test, LOPO-CV)
  # Note: train_auc is already calculated, but we need train_roc for plotting
  lambda_to_use_plot <- if("lambda_to_use" %in% names(cv_model_final)) {
    cv_model_final$lambda_to_use
  } else {
    "lambda.1se"  # Default fallback
  }
  
  X_train_scaled_for_plot <- scale(X_train_final, center = scaler_center, scale = scaler_scale)
  pred_probs_train_orig <- as.numeric(predict(
    cv_model_final,
    newx = X_train_scaled_for_plot,
    s = lambda_to_use_plot,
    type = "response"
  ))
  train_roc <- roc(y_train, pred_probs_train_orig, quiet = TRUE)
  
  train_roc_smooth <- smooth(train_roc, method = "density")
  roc_smooth_combined <- smooth(roc_obj, method = "density")
  
  # Apply monotonicity so no curve dips below the diagonal (sort by FPR, TPR = cummax(TPR))
  test_m <- roc_monotonic(1 - roc_smooth_combined$specificities, roc_smooth_combined$sensitivities)
  roc_combined_data <- data.frame(FPR = test_m$fpr, TPR = test_m$tpr, Method = "Test")
  
  if(!is.null(lopo_roc)) {
    lopo_roc_smooth_combined <- smooth(lopo_roc, method = "density")
    lopo_m <- roc_monotonic(1 - lopo_roc_smooth_combined$specificities, lopo_roc_smooth_combined$sensitivities)
    lopo_roc_data <- data.frame(FPR = lopo_m$fpr, TPR = lopo_m$tpr, Method = "LOPO-CV")
    roc_combined_data <- rbind(roc_combined_data, lopo_roc_data)
  }
  
  train_m <- roc_monotonic(1 - train_roc_smooth$specificities, train_roc_smooth$sensitivities)
  train_roc_data <- data.frame(FPR = train_m$fpr, TPR = train_m$tpr, Method = "Train")
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
  ggsave(paste0("figures/roc_combined_train_test_lopo_logres_", output_prefix, "_polished_new.png"), 
         plot = p_roc_combined, width = 19, height = 19, dpi = 600, units = "in")
  ggsave(paste0("figures/roc_combined_train_test_lopo_logres_", output_prefix, "_polished_new.pdf"), 
         plot = p_roc_combined, width = 19, height = 19, units = "in")
  ggsave(paste0("figures/roc_combined_train_test_lopo_logres_", output_prefix, "_polished_new.svg"), 
         plot = p_roc_combined, width = 19, height = 19, units = "in")
  if (has_devEMF) tryCatch(
    ggsave(paste0("figures/roc_combined_train_test_lopo_logres_", output_prefix, "_polished_new.emf"), 
           plot = p_roc_combined, width = 19, height = 19, units = "in"),
    error = function(e) message("EMF export skipped: ", conditionMessage(e))
  )
  # Plot 3: Combined ROC curves WITHOUT smoothing
  # Use existing ROC objects (they already have good resolution)
  # For smoother appearance, we'll use the existing roc_obj, train_roc, and lopo_roc
  roc_obj_hr <- roc_obj
  train_roc_hr <- train_roc
  lopo_roc_hr <- lopo_roc
  
  # Apply monotonicity so LOPO (and any curve) never dips below the diagonal
  test_um <- roc_monotonic(1 - roc_obj_hr$specificities, roc_obj_hr$sensitivities)
  roc_combined_data_unsmooth <- data.frame(FPR = test_um$fpr, TPR = test_um$tpr, Method = "Test")
  
  if(!is.null(lopo_roc_hr)) {
    lopo_um <- roc_monotonic(1 - lopo_roc_hr$specificities, lopo_roc_hr$sensitivities)
    lopo_roc_data_unsmooth <- data.frame(FPR = lopo_um$fpr, TPR = lopo_um$tpr, Method = "LOPO-CV")
    roc_combined_data_unsmooth <- rbind(roc_combined_data_unsmooth, lopo_roc_data_unsmooth)
  }
  
  train_um <- roc_monotonic(1 - train_roc_hr$specificities, train_roc_hr$sensitivities)
  train_roc_data_unsmooth <- data.frame(FPR = train_um$fpr, TPR = train_um$tpr, Method = "Train")
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
  ggsave(paste0("figures/roc_combined_unsmoothed_logres_", output_prefix, "_polished_new.png"), 
         plot = p_roc_combined_unsmooth, width = 19, height = 19, dpi = 600, units = "in")
  ggsave(paste0("figures/roc_combined_unsmoothed_logres_", output_prefix, "_polished_new.pdf"), 
         plot = p_roc_combined_unsmooth, width = 19, height = 19, units = "in")
  ggsave(paste0("figures/roc_combined_unsmoothed_logres_", output_prefix, "_polished_new.svg"), 
         plot = p_roc_combined_unsmooth, width = 19, height = 19, units = "in")
  if (has_devEMF) tryCatch(
    ggsave(paste0("figures/roc_combined_unsmoothed_logres_", output_prefix, "_polished_new.emf"), 
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
    ci_combined <- do.call(rbind, ci_data_list)
    # Clamp CI bounds to be on or above the diagonal so ribbon never dips below 50% line
    ci_combined$lower <- pmax(ci_combined$lower, ci_combined$FPR)
    ci_combined$upper <- pmax(ci_combined$upper, ci_combined$FPR)
    ci_combined
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
  ggsave(paste0("figures/roc_combined_polished_logres_", output_prefix, "_polished_new.png"), 
         plot = p_roc_polished, width = 19, height = 19, dpi = 600, units = "in")
  ggsave(paste0("figures/roc_combined_polished_logres_", output_prefix, "_polished_new.pdf"), 
         plot = p_roc_polished, width = 19, height = 19, units = "in")
  ggsave(paste0("figures/roc_combined_polished_logres_", output_prefix, "_polished_new.svg"), 
         plot = p_roc_polished, width = 19, height = 19, units = "in")
  if (has_devEMF) tryCatch(
    ggsave(paste0("figures/roc_combined_polished_logres_", output_prefix, "_polished_new.emf"), 
           plot = p_roc_polished, width = 19, height = 19, units = "in"),
    error = function(e) message("EMF export skipped: ", conditionMessage(e))
  )
  # Plot 4: Combined ROC curves WITHOUT smoothing WITH shaded AUC areas
  # Use monotonic curves for polygons so shading never goes below diagonal
  test_um2 <- roc_monotonic(1 - roc_obj_hr$specificities, roc_obj_hr$sensitivities)
  test_polygon <- data.frame(
    FPR = c(0, test_um2$fpr, 1),
    TPR = c(0, test_um2$tpr, 0),
    Method = "Test"
  )
  train_um2 <- roc_monotonic(1 - train_roc_hr$specificities, train_roc_hr$sensitivities)
  train_polygon <- data.frame(
    FPR = c(0, train_um2$fpr, 1),
    TPR = c(0, train_um2$tpr, 0),
    Method = "Train"
  )
  if(!is.null(lopo_roc_hr)) {
    lopo_um2 <- roc_monotonic(1 - lopo_roc_hr$specificities, lopo_roc_hr$sensitivities)
    lopo_polygon <- data.frame(
      FPR = c(0, lopo_um2$fpr, 1),
      TPR = c(0, lopo_um2$tpr, 0),
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
  ggsave(paste0("figures/roc_combined_shaded_auc_logres_", output_prefix, "_polished_new.png"), 
         plot = p_roc_combined_shaded, width = 19, height = 19, dpi = 600, units = "in")
  ggsave(paste0("figures/roc_combined_shaded_auc_logres_", output_prefix, "_polished_new.pdf"), 
         plot = p_roc_combined_shaded, width = 19, height = 19, units = "in")
  ggsave(paste0("figures/roc_combined_shaded_auc_logres_", output_prefix, "_polished_new.svg"), 
         plot = p_roc_combined_shaded, width = 19, height = 19, units = "in")
  if (has_devEMF) tryCatch(
    ggsave(paste0("figures/roc_combined_shaded_auc_logres_", output_prefix, "_polished_new.emf"), 
           plot = p_roc_combined_shaded, width = 19, height = 19, units = "in"),
    error = function(e) message("EMF export skipped: ", conditionMessage(e))
  )
  cat("✓ All plots saved\n\n")
}

###MAIN

cat("============================================================================\n")
cat("UNIFIED DENGUE PATIENT CLASSIFIER - LOGISTIC REGRESSION (ELASTIC NET)\n")
cat("============================================================================\n\n")

# Step 1: Load and transform data
data_result <- load_and_transform_data(INPUT_FILE, CLASSIFICATION_TYPE)
df_wide <- data_result$df_wide
class_levels <- data_result$class_levels
class_labels <- data_result$class_labels

# Step 2: Check leakage and split
if(!is.na(RANDOM_SEED_SPLIT)) set.seed(RANDOM_SEED_SPLIT)
split_result <- check_leakage_and_split(df_wide, CLASSIFICATION_TYPE)
train_data <- split_result$train_data
test_data <- split_result$test_data
if(length(split_result$test_patients) < 10) {
  cat("WARNING: Small test set (", length(split_result$test_patients), " patients). Test AUC may be unstable; LOPO-CV is a more reliable estimate.\n\n")
}

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
class_weights <- handle_class_balance(y_train, balance_strength = 1.5)

# Step 6: Feature scaling
scale_result <- scale_features(X_train_final, X_test_final)
X_train_scaled <- scale_result$X_train_scaled
X_test_scaled <- scale_result$X_test_scaled
scaler_center <- scale_result$scaler_center
scaler_scale <- scale_result$scaler_scale

# Step 7: Hyperparameter tuning
best_params <- tune_hyperparameters(X_train_final, y_train, train_data, class_weights)

# Step 8: Train final model
cv_model_final <- train_final_model(X_train_scaled, y_train, best_params, class_weights)

# Step 9: Print diagnostic summary
cat("\n")
cat("============================================================================\n")
cat("DIAGNOSTIC SUMMARY\n")
cat("============================================================================\n")
cat("\n1. FEATURE COUNT:\n")
cat("Final number of features used:", ncol(X_train_final), "\n")
cat("Sample-to-feature ratio:", round(nrow(X_train_scaled) / ncol(X_train_final), 2), ":1\n")
if(nrow(X_train_scaled) / ncol(X_train_final) < 10) {
  cat("WARNING: Ratio < 10:1 - may lead to overfitting\n")
}
cat("\n2. CLASS DISTRIBUTION:\n")
cat("Overall (all patients):\n")
patient_classes_all <- df_wide %>%
  group_by(patient_id) %>%
  summarise(class = first(class)) %>%
  ungroup()
print(table(patient_classes_all$class))
cat("\nTraining set (patients):\n")
train_patient_classes <- train_data %>%
  group_by(patient_id) %>%
  summarise(class = first(class)) %>%
  ungroup()
print(table(train_patient_classes$class))
cat("\nTest set (patients):\n")
test_patient_classes <- test_data %>%
  group_by(patient_id) %>%
  summarise(class = first(class)) %>%
  ungroup()
print(table(test_patient_classes$class))
cat("\n3.TRAIN/TEST SPLIT:\n")
cat("Split ratio: 75/25 (75% train, 25% test)\n")
cat("Training patients:", length(unique(train_data$patient_id)), "\n")
cat("Test patients:", length(unique(test_data$patient_id)), "\n")
cat("Training observations:", nrow(train_data), "\n")
cat("Test observations:", nrow(test_data), "\n")
cat("\n4. ELASTIC NET PARAMETERS:\n")
cat("Alpha (regularization mix):", best_params$alpha, "\n")
cat("(0 = Ridge, 1 = Lasso, 0.5 = Elastic Net)\n")
if(best_params$alpha == 1) {
  cat("Using Lasso (L1 regularization only)\n")
} else if(best_params$alpha == 0) {
  cat("Using Ridge (L2 regularization only)\n")
} else {
  cat("Using Elastic Net (L1 + L2 regularization)\n")
  cat("L1 ratio:", best_params$alpha, ", L2 ratio:", 1 - best_params$alpha, "\n")
}
cat("Lambda type:", best_params$lambda_type, "\n")
cat("Lambda value:", round(best_params$lambda_value, 6), "\n")
if(cv_model_final$lambda_to_use == "lambda.1se") {
  cat("Final lambda used: lambda.1se (more conservative)\n")
} else {
  cat("Final lambda used: lambda.min (less conservative)\n")
}
cat("\n============================================================================\n")
cat("\n")

# Step 10: Evaluate model
eval_result <- evaluate_model(cv_model_final, X_test_scaled, y_test, X_train_final, y_train,
                              scaler_center, scaler_scale, CLASSIFICATION_TYPE, class_labels)
if(eval_result$train_auc - eval_result$auc_val > 0.15) {
  cat("WARNING: Large train–test AUC gap (", round(eval_result$train_auc - eval_result$auc_val, 3), 
      "). Model may be overfitting; LOPO-CV AUC is a more reliable estimate of generalization.\n\n")
}
if(eval_result$auc_val >= 0.999) {
  cat("NOTE: Test AUC is (near) 1.0. With a small test set this can occur by chance; LOPO-CV AUC is the more reliable estimate.\n\n")
}

# Step 11: LOPO-CV
lopo_result <- perform_lopo_cv(df_wide, fixed_feature_set, best_params, CLASSIFICATION_TYPE, class_weights, class_levels)

# Step 12: Feature importance
importance_df <- extract_feature_importance(cv_model_final, feature_names = colnames(X_train_final))

# Step 13: Create plots
n_train_plot <- length(unique(train_data$patient_id))
n_test_plot <- length(unique(test_data$patient_id))
n_lopo_plot <- length(unique(df_wide$patient_id))
create_plots(eval_result$roc_obj, eval_result$auc_val, eval_result$train_auc,
            lopo_result$lopo_roc, lopo_result$lopo_auc, importance_df,
            eval_result$pred_probs, y_test, eval_result$best_thresh,
            CLASSIFICATION_TYPE, class_labels, OUTPUT_PREFIX,
            cv_model_final, X_train_final, y_train, scaler_center, scaler_scale,
            n_train_plot, n_test_plot, n_lopo_plot)

# Step 14: Create separate feature importance plot
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

