# ============================================================================
# STRATIFIED CROSS-PLATFORM PAIR LEARNING
# Learn pairs from EBPP + 50% Xena, test on held-out 50% Xena
# ============================================================================

library(data.table)
library(xgboost)
library(ImmuneSubtypeClassifier)

# ============================================================================
# STEP 1: Load cleaned datasets
# ============================================================================

ebpp_clean <- fread('../data/formatted_full_L1000/EBpp_pancancer_matched.csv.gz')
xena_clean <- fread('../data/formatted_full_L1000/xena_rsem_tpm_matched.csv.gz')

cat("Loaded cleaned datasets:\n")
cat("EBPP:", nrow(ebpp_clean), "x", ncol(ebpp_clean), "\n")
cat("Xena:", nrow(xena_clean), "x", ncol(xena_clean), "\n\n")

# Verify columns match
if (!all(colnames(ebpp_clean) == colnames(xena_clean))) {
  stop("Column names don't match between EBPP and Xena!")
}

# ============================================================================
# STEP 1.5: LOG-TRANSFORM EBPP TO MATCH XENA PIPELINE
# Xena uses: log2(TPM + 0.001)
# ============================================================================

cat("=== Transforming EBPP to match Xena pipeline ===\n")

# Identify metadata columns (these will NOT be transformed)
metadata_cols <- c("SampleID", "Barcode", "Label", "label")
existing_metadata <- intersect(colnames(ebpp_clean), metadata_cols)

cat("Metadata columns (will be preserved):", paste(existing_metadata, collapse=", "), "\n")

# Gene columns are everything else
gene_cols <- setdiff(colnames(ebpp_clean), metadata_cols)
cat("Gene columns to transform:", length(gene_cols), "\n\n")

# Check current data ranges
cat("Before transformation:\n")
ebpp_sample_genes <- sample(gene_cols, min(5, length(gene_cols)))
cat("EBPP sample genes:\n")
for (g in ebpp_sample_genes) {
  cat(sprintf("  %s: [%.2f, %.2f]\n", g,
              min(ebpp_clean[[g]], na.rm=TRUE),
              max(ebpp_clean[[g]], na.rm=TRUE)))
}

xena_sample_genes <- sample(gene_cols, min(5, length(gene_cols)))
cat("\nXena sample genes:\n")
for (g in xena_sample_genes) {
  cat(sprintf("  %s: [%.2f, %.2f]\n", g,
              min(xena_clean[[g]], na.rm=TRUE),
              max(xena_clean[[g]], na.rm=TRUE)))
}

# Done in matching code
# # Log-transform EBPP gene columns ONLY
# # Xena pipeline: log2(TPM + 0.001)
# cat("\nApplying log2(TPM + 0.001) transformation to EBPP...\n")
# for (gene in gene_cols) {
#   ebpp_clean[[gene]] <- log2(ebpp_clean[[gene]] + 0.001)
# }
#
# cat("\nAfter log2(x + 0.001) transformation:\n")
# cat("EBPP sample genes:\n")
# for (g in ebpp_sample_genes) {
#   cat(sprintf("  %s: [%.2f, %.2f]\n", g,
#               min(ebpp_clean[[g]], na.rm=TRUE),
#               max(ebpp_clean[[g]], na.rm=TRUE)))
# }

# Verify metadata columns were NOT transformed
cat("\nMetadata verification:\n")
for (mc in existing_metadata) {
  if (mc == "Label") {
    cat(sprintf("  %s: %s (preserved)\n", mc,
                paste(head(unique(ebpp_clean[[mc]]), 3), collapse=", ")))
  } else if (mc %in% c("SampleID", "Barcode")) {
    cat(sprintf("  %s: %s... (preserved)\n", mc,
                paste(head(ebpp_clean[[mc]], 2), collapse=", ")))
  }
}

cat("\n✓ EBPP transformed to log2 scale, metadata preserved\n\n")

# ============================================================================
# STEP 2: Stratified split - Hold out 50% of Xena for testing
# ============================================================================

set.seed(412)

# Shuffle Xena
xena_shuffled <- xena_clean[sample(nrow(xena_clean)), ]

# Split 50/50
xena_train_size <- floor(0.5 * nrow(xena_shuffled))
xena_train <- xena_shuffled[1:xena_train_size, ]
xena_test <- xena_shuffled[(xena_train_size + 1):nrow(xena_shuffled), ]

cat("=== Stratified Split ===\n")
cat("EBPP (all for training, log-transformed):", nrow(ebpp_clean), "\n")
cat("Xena training (50%):", nrow(xena_train), "\n")
cat("Xena test holdout (50%):", nrow(xena_test), "\n\n")

# Verify metadata in splits
cat("Metadata columns in all datasets:\n")
cat("  EBPP:", paste(intersect(colnames(ebpp_clean), metadata_cols), collapse=", "), "\n")
cat("  Xena train:", paste(intersect(colnames(xena_train), metadata_cols), collapse=", "), "\n")
cat("  Xena test:", paste(intersect(colnames(xena_test), metadata_cols), collapse=", "), "\n\n")

# ============================================================================
# STEP 3: Create training set (ALL EBPP + 50% Xena)
# ============================================================================

# Add batch indicators
ebpp_clean$Batch <- "EBPP"
xena_train$Batch <- "Xena_train"
xena_test$Batch <- "Xena_test"

# Combine for training
train_combined <- rbind(ebpp_clean, xena_train)

# Shuffle training data
train_combined <- train_combined[sample(nrow(train_combined)), ]

cat("=== Training Data Composition ===\n")
cat("Total training samples:", nrow(train_combined), "\n")
cat("  EBPP (log2):", sum(train_combined$Batch == "EBPP"), "\n")
cat("  Xena:", sum(train_combined$Batch == "Xena_train"), "\n\n")

cat("Label distribution in training:\n")
print(table(train_combined$Label))
cat("\nBatch distribution in training:\n")
print(table(train_combined$Batch))
cat("\n")

# Final metadata check
final_metadata <- intersect(colnames(train_combined), c(metadata_cols, "Batch"))
cat("Final metadata columns in training set:\n")
cat("  ", paste(final_metadata, collapse=", "), "\n\n")

# Save training set, test set,
fwrite(train_combined, '../data/formatted_full_L1000/train_stratified.csv.gz')
fwrite(xena_test, '../data/formatted_full_L1000test_xena_holdout.csv.gz')

rm(ebpp_clean, xena_clean, xena_train, xena_shuffled)

# ============================================================================
# STEP 4: Select top variance genes from TRAINING data only
# ============================================================================

# Define metadata columns (including Batch now)
all_metadata_cols <- c("SampleID", "Barcode", "Label", "label", "Batch")
numeric_cols <- setdiff(colnames(train_combined), all_metadata_cols)

cat("Gene columns available:", length(numeric_cols), "\n")

# Calculate variance on TRAINING data only (no data leakage)
gene_variances <- apply(train_combined[, ..numeric_cols], 2, var, na.rm = TRUE)

# Select top variance genes
top_300_genes <- names(sort(gene_variances, decreasing = TRUE)[1:300])

cat("\nTop 10 genes by variance (from training data):\n")
print(head(top_300_genes, 10))
cat("\n")

# ============================================================================
# STEP 5: Create all pairs from TRAINING data only
# ============================================================================

create_all_pairs_from_genes <- function(data, gene_list) {
  gene_cols <- gene_list[gene_list %in% colnames(data)]
  if (length(gene_cols) < 2) {
    stop("Need at least 2 genes to create pairs")
  }

  pair_features <- list()

  for (i in 1:(length(gene_cols)-1)) {
    for (j in (i+1):length(gene_cols)) {
      gene1 <- gene_cols[i]
      gene2 <- gene_cols[j]
      pair_feat <- as.numeric(data[[gene1]] > data[[gene2]])
      pair_name <- paste0(gene1, "_X_", gene2)
      pair_features[[pair_name]] <- pair_feat
    }
  }

  pair_data <- as.data.table(pair_features)
  cat(sprintf("Created %d pairs from %d genes\n",
              ncol(pair_data), length(gene_cols)))
  return(pair_data)
}

# Create pairs from TRAINING data
all_pairs_train <- create_all_pairs_from_genes(train_combined, top_300_genes)

cat("Pair data dimensions:", nrow(all_pairs_train), "samples x",
    ncol(all_pairs_train), "pairs\n\n")

saveRDS(all_pairs_train, '../data/formatted_full_L1000/training/all_pairs_train.rds')

# ============================================================================
# STEP 6: Rank pairs by XGBoost importance on TRAINING data
# ============================================================================

rank_pairs_xgboost <- function(pair_data, labels, top_n = 200,
                               nrounds = 100, early_stopping = 10) {

  if (is.factor(labels)) {
    label_map <- levels(labels)
    labels_numeric <- as.numeric(labels) - 1
  } else {
    label_map <- sort(unique(labels))
    labels_numeric <- match(labels, label_map) - 1
  }
  num_classes <- length(unique(labels_numeric))

  cat(sprintf("Training XGBoost with %d pairs, %d samples, %d classes\n",
              ncol(pair_data), nrow(pair_data), num_classes))

  # Train/validation split within training data
  set.seed(42)
  n_samples <- nrow(pair_data)
  train_idx <- sample(1:n_samples, size = floor(0.8 * n_samples))
  val_idx <- setdiff(1:n_samples, train_idx)

  dtrain <- xgb.DMatrix(
    data = as.matrix(pair_data[train_idx, ]),
    label = labels_numeric[train_idx]
  )
  dval <- xgb.DMatrix(
    data = as.matrix(pair_data[val_idx, ]),
    label = labels_numeric[val_idx]
  )

  params <- list(
    objective = "multi:softmax",
    num_class = num_classes,
    max_depth = 4,
    eta = 0.3,
    subsample = 0.8,
    colsample_bytree = 0.8,
    eval_metric = "mlogloss"
  )

  watchlist <- list(train = dtrain, val = dval)
  xgb_model <- xgb.train(
    params = params,
    data = dtrain,
    nrounds = nrounds,
    watchlist = watchlist,
    early_stopping_rounds = early_stopping,
    verbose = 1
  )

  importance_matrix <- xgb.importance(
    feature_names = colnames(pair_data),
    model = xgb_model
  )

  importance_matrix$gene1 <- sapply(strsplit(importance_matrix$Feature, "_X_"), `[`, 1)
  importance_matrix$gene2 <- sapply(strsplit(importance_matrix$Feature, "_X_"), `[`, 2)

  top_pairs <- head(importance_matrix, top_n)

  cat(sprintf("\nTop %d pairs selected from %d total\n",
              nrow(top_pairs), nrow(importance_matrix)))
  cat(sprintf("Top pair: %s (importance: %.4f)\n",
              top_pairs$Feature[1], top_pairs$Gain[1]))

  return(top_pairs)
}

# Learn pairs from training data (log-EBPP + 50% Xena)
cat("=== Learning pairs from STRATIFIED TRAINING SET ===\n")
cat("Training on:", nrow(train_combined), "samples (log2-EBPP + 50% Xena)\n")
cat("Held out:", nrow(xena_test), "samples (50% Xena for testing)\n\n")

top_pairs <- rank_pairs_xgboost(
  pair_data = all_pairs_train,
  labels = train_combined$Label,
  top_n = 100,
  nrounds = 50,
  early_stopping = 5
)

# ============================================================================
# STEP 7: Convert to pair_list format
# ============================================================================

convert_to_pair_list <- function(ranked_pairs) {
  pair_list <- c()
  for (i in 1:nrow(ranked_pairs)) {
    pair_list <- c(pair_list, ranked_pairs$gene1[i], ranked_pairs$gene2[i])
  }
  return(pair_list)
}

stratified_pair_list <- convert_to_pair_list(top_pairs)

cat("\nStratified cross-platform pair list created:\n")
cat("Length:", length(stratified_pair_list), "(", length(stratified_pair_list)/2, "pairs )\n")
cat("First 10 genes:\n")
print(head(stratified_pair_list, 10))

# ============================================================================
# STEP 8: Verify genes exist and save everything
# ============================================================================

genes_in_pairs <- unique(stratified_pair_list)
common_genes_check <- setdiff(colnames(train_combined), all_metadata_cols)
missing_genes <- setdiff(genes_in_pairs, common_genes_check)

cat("\nValidation:\n")
cat("Unique genes in pair_list:", length(genes_in_pairs), "\n")
cat("Missing genes:", length(missing_genes), "\n")

if (length(missing_genes) > 0) {
  cat("ERROR: These genes are missing:\n")
  print(missing_genes)
  stop("Fix missing genes before proceeding")
} else {
  cat("✓ All genes in pair_list exist in datasets!\n")
}


saveRDS(stratified_pair_list, '../models/pair_list_stratified.rds')
saveRDS(top_pairs, '../models/top_pairs_stratified.rds')

cat("\n=== SAVED FILES ===\n")
cat("✓ Training data: ../data/formatted_full_L1000/train_stratified.csv.gz\n")
cat("   (", nrow(train_combined), "samples:", sum(train_combined$Batch == "EBPP"),
    "log2-EBPP +", sum(train_combined$Batch == "Xena_train"), "Xena )\n")
cat("   Columns:", paste(head(colnames(train_combined), 10), collapse=", "), "...\n")
cat("✓ Test data: ../data/formatted_full_L1000/test_xena_holdout.csv.gz\n")
cat("   (", nrow(xena_test), "Xena holdout samples )\n")
cat("   Columns:", paste(head(colnames(xena_test), 10), collapse=", "), "...\n")
cat("✓ Pair list: ../models/pair_list_stratified.rds\n")
cat("✓ Top pairs: ../models/top_pairs_stratified.rds\n")

cat("\n=== SUCCESS ===\n")
cat("EBPP was transformed using Xena's exact method: log2(TPM + 0.001)\n")
cat("Metadata preserved:", paste(final_metadata, collapse=", "), "\n")
cat("Pairs were learned from", nrow(train_combined), "training samples.\n")
cat("These include", sum(train_combined$Batch == "EBPP"), "log2-EBPP +",
    sum(train_combined$Batch == "Xena_train"), "Xena samples.\n")
cat("", nrow(xena_test), "Xena samples held out for unbiased testing.\n")
cat("Pairs should be robust across both pipelines!\n")

# ============================================================================
# STEP 9: RUN IT
# ============================================================================

# Conservative parameters for better C4/C6
conservative_params <- list(
  max_depth = 8,
  eta = 0.3,
  nrounds = 64,
  early_stopping_rounds = 4,
  gamma = 0.2,
  lambda = 1.8,
  alpha = 0.2,
  ensemble_size = 11,
  sample_prop = 0.8,
  feature_prop = 0.8,
  subsample = 0.8
)

stratified_pair_list <- readRDS('../models/pair_list_stratified.rds')
this_set <- editPairList(stratified_pair_list, 'IGJ')
#this_set <- stratified_pair_list[1:64]

# Now use the cleaned pair_list
result <- build_robencla_classifier(
  data_path='../data/formatted_full_L1000/train_stratified.csv.gz',
  test_path='../data/formatted_full_L1000test_xena_holdout.csv.gz',
  output_path = '../models/immune_optimized_pairs_64.rds',
  pair_list = this_set,
  sig_list = NULL,
  param_list = conservative_params,
  data_mode = c("namedpairs"),
  train_fraction = NULL,
  seed = 412,
  sample_id = "SampleID"  # Specify the sample ID column
)


# ============================================================================
# STEP 10: Train with only EBpp, but the combined pair_list
# ============================================================================

# Conservative parameters for better C4/C6
conservative_params <- list(
  max_depth = 8,
  eta = 0.3,
  nrounds = 64,
  early_stopping_rounds = 4,
  gamma = 0.2,
  lambda = 1.8,
  alpha = 0.2,
  ensemble_size = 11,
  sample_prop = 0.8,
  feature_prop = 0.8,
  subsample = 0.8
)

stratified_pair_list <- readRDS('../models/pair_list_stratified.rds')
this_set <- stratified_pair_list[1:64]

new_list <- editPairList(this_set, 'IGJ')

# Now use the cleaned pair_list
result <- build_robencla_classifier(
  data_path='../data/formatted_full_L1000/EBpp_pancancer_matched.csv.gz',
  test_path=NULL, #'../data/formatted_full_L1000/xena_rsem_tpm.csv.gz',
  output_path = '../models/immune_optimized_pairs_100.rds',
  pair_list = new_list,
  sig_list = NULL,
  param_list = conservative_params,
  data_mode = c("namedpairs"),
  train_fraction = 1.0,
  seed = 412,
  sample_id = "Barcode"  # Specify the sample ID column
)

library(ImmuneSubtypeClassifier)
Xs <- readr::read_csv('../data/formatted_full_L1000/xena_rsem_tpm.csv.gz')
result <- callSubtypes(Xs,
                       model = NULL,
                       model_path = '../models/immune_optimized_pairs_100.rds',
                       geneid = "symbol",
                       sampleid = 'Barcode',
                       labelid='Label')

library(caret)

# Create confusion matrix object
cm <- confusionMatrix(
  factor(result$BestCall, levels = 1:6),
  factor(result$Label, levels = 1:6)
)

# Get key metrics
cm$overall  # Overall accuracy, Kappa, etc.
cm$byClass  # Per-class metrics

# Manual calculations:
conf_table <- table(result$Label, result$BestCall)

# 1. Overall Accuracy
accuracy <- sum(diag(conf_table)) / sum(conf_table)
# (1638 + 2371 + 2022 + 856 + 329 + 61) / total ≈ 0.88

# 2. Per-class metrics
per_class <- sapply(1:6, function(i) {
  TP <- conf_table[i, i]
  FP <- sum(conf_table[, i]) - TP
  FN <- sum(conf_table[i, ]) - TP
  TN <- sum(conf_table) - TP - FP - FN

  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)  # sensitivity
  f1 <- 2 * (precision * recall) / (precision + recall)

  c(Precision = precision, Recall = recall, F1 = f1)
})
colnames(per_class) <- paste0("Class", 1:6)
t(per_class)

# 3. Macro-averaged F1 (unweighted mean across classes)
macro_f1 <- mean(per_class["F1", ])

# 4. Weighted F1 (weighted by class size)
class_sizes <- rowSums(conf_table)
weighted_f1 <- sum(per_class["F1", ] * class_sizes) / sum(class_sizes)

# 5. Cohen's Kappa (accounts for chance agreement)
kappa <- cm$overall["Kappa"]

# Quick summary
list(
  Accuracy = accuracy,
  Macro_F1 = macro_f1,
  Weighted_F1 = weighted_f1,
  Kappa = kappa
)
