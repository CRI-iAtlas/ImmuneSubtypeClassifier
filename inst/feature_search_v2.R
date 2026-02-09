

# ============================================================================
# REQUIRED LIBRARIES
# ============================================================================

# Core data manipulation
library(data.table)      # Fast data reading/writing and manipulation

# Machine learning
library(xgboost)         # For XGBoost models and feature importance ranking

# Your custom packages
library(ImmuneSubtypeClassifier)  # For genesetsymbols if using signatures

# Optional but recommended
library(dplyr)           # For data wrangling (if you use %>% pipes)
library(stringr)         # For string manipulation (str_detect, etc.)

# If using ROCit for signature pairs
library(ROCit)           # For AUC calculation in sigpairs mode

# R6 for object-oriented programming (should be auto-loaded with robencla)
library(R6)              # For R6 class definitions

# ============================================================================
# INSTALLATION COMMANDS (if needed)
# ============================================================================

# Core packages
# install.packages("data.table")
# install.packages("xgboost")
# install.packages("dplyr")
# install.packages("stringr")
# install.packages("R6")
# install.packages("ROCit")

# Your custom package
# devtools::install_github("Gibbsdavidl/ImmuneSubtypeClassifier")

# ============================================================================
# MINIMAL REQUIRED SET FOR YOUR CURRENT WORKFLOW
# ============================================================================

library(data.table)
library(xgboost)
library(ImmuneSubtypeClassifier)  # If using genesetsymbols

# ============================================================================
# CLEAN GENE NAME NORMALIZATION FOR BOTH DATASETS
# ============================================================================

# Load both datasets
ebpp <- data.table::fread('../data/formatted_full_L1000/EBpp_pancancer.csv.gz')
xena <- data.table::fread('../data/formatted_full_L1000/xena_rsem_tpm.csv.gz')

cat("Original dimensions:\n")
cat("EBPP:", nrow(ebpp), "x", ncol(ebpp), "\n")
cat("Xena:", nrow(xena), "x", ncol(xena), "\n\n")

# ============================================================================
# STEP 1: Define a consistent gene name normalization function
# ============================================================================

normalize_gene_names <- function(gene_names) {
  # Apply transformations in a consistent order
  gene_names <- gsub(" ", "_", gene_names)      # spaces to underscores
  gene_names <- gsub("-", "_", gene_names)      # hyphens to underscores
  gene_names <- gsub("\\.", "_", gene_names)    # dots to underscores
  return(gene_names)
}

# ============================================================================
# STEP 2: Normalize column names in both datasets
# ============================================================================

# Save original names for reference
ebpp_original_names <- colnames(ebpp)
xena_original_names <- colnames(xena)

# Apply normalization
colnames(ebpp) <- normalize_gene_names(colnames(ebpp))
colnames(xena) <- normalize_gene_names(colnames(xena))

cat("After normalization:\n")
cat("EBPP columns:", ncol(ebpp), "\n")
cat("Xena columns:", ncol(xena), "\n\n")

# ============================================================================
# STEP 3: Find common genes between datasets
# ============================================================================

# Identify special columns (non-gene columns)
special_cols_ebpp <- c("Label", "Barcode")
special_cols_xena <- c("Label", "Barcode")

# Gene columns only
ebpp_genes <- setdiff(colnames(ebpp), special_cols_ebpp)
xena_genes <- setdiff(colnames(xena), special_cols_xena)

cat("Gene columns:\n")
cat("EBPP:", length(ebpp_genes), "genes\n")
cat("Xena:", length(xena_genes), "genes\n\n")

# Find common genes
common_genes <- intersect(ebpp_genes, xena_genes)
cat("Common genes:", length(common_genes), "\n\n")

# Find genes unique to each dataset
ebpp_only <- setdiff(ebpp_genes, xena_genes)
xena_only <- setdiff(xena_genes, ebpp_genes)

cat("Genes only in EBPP:", length(ebpp_only), "\n")
if (length(ebpp_only) > 0) {
  cat("Examples:", head(ebpp_only, 10), "\n")
}
cat("\n")

cat("Genes only in Xena:", length(xena_only), "\n")
if (length(xena_only) > 0) {
  cat("Examples:", head(xena_only, 10), "\n")
}
cat("\n")

# ============================================================================
# STEP 4: Handle special cases (like IGJ -> IGJP1)
# ============================================================================

# Check for genes that might be related but named differently
potential_mappings <- list()

# Example: Check if IGJ is in EBPP but only IGJP1 is in Xena
if ("IGJ" %in% ebpp_genes && !("IGJ" %in% xena_genes) && "IGJP1" %in% xena_genes) {
  cat("Found mapping: IGJ (EBPP) -> IGJP1 (Xena)\n")
  potential_mappings[["IGJ"]] <- "IGJP1"

  # Decide: Do we want to rename IGJP1 to IGJ in Xena?
  # Option 1: Rename in Xena (treat as equivalent)
  colnames(xena)[colnames(xena) == "IGJP1"] <- "IGJ"
  common_genes <- c(common_genes, "IGJ")
  xena_only <- setdiff(xena_only, "IGJP1")
  cat("-> Renamed IGJP1 to IGJ in Xena\n\n")
}

# ============================================================================
# STEP 5: Subset both datasets to common genes only
# ============================================================================

# Keep label column if it exists
ebpp_label_col <- intersect(colnames(ebpp), c("Label", "label"))
xena_label_col <- intersect(colnames(xena), c("Label", "label"))

# Subset to common genes + label
if (length(ebpp_label_col) > 0) {
  ebpp_clean <- ebpp[, c(ebpp_label_col, common_genes), with = FALSE]
} else {
  ebpp_clean <- ebpp[, common_genes, with = FALSE]
}

if (length(xena_label_col) > 0) {
  xena_clean <- xena[, c(xena_label_col, common_genes), with = FALSE]
} else {
  xena_clean <- xena[, common_genes, with = FALSE]
}

cat("Cleaned dimensions:\n")
cat("EBPP:", nrow(ebpp_clean), "x", ncol(ebpp_clean), "\n")
cat("Xena:", nrow(xena_clean), "x", ncol(xena_clean), "\n\n")

# ============================================================================
# STEP 6: Verify column names match exactly
# ============================================================================

if (all(colnames(ebpp_clean) == colnames(xena_clean))) {
  cat("✓ Column names match perfectly!\n\n")
} else {
  cat("✗ Column names still don't match:\n")
  mismatches <- which(colnames(ebpp_clean) != colnames(xena_clean))
  cat("Mismatched positions:", mismatches, "\n")
  print(data.frame(
    Position = mismatches,
    EBPP = colnames(ebpp_clean)[mismatches],
    Xena = colnames(xena_clean)[mismatches]
  ))
}

# ============================================================================
# STEP 7: Save cleaned datasets
# ============================================================================

data.table::fwrite(ebpp_clean, '../data/formatted_full_L1000/EBpp_pancancer_cleaned.csv.gz')
data.table::fwrite(xena_clean, '../data/formatted_full_L1000/xena_rsem_tpm_cleaned.csv.gz')

cat("\nCleaned datasets saved:\n")
cat("- EBpp_pancancer_cleaned.csv.gz\n")
cat("- xena_rsem_tpm_cleaned.csv.gz\n\n")



# ============================================================================
# STEP 1: Load cleaned datasets
# ============================================================================

library(data.table)
library(xgboost)

ebpp_clean <- data.table::fread('../data/formatted_full_L1000/EBpp_pancancer_cleaned.csv.gz')

cat("Loaded cleaned EBPP data:", nrow(ebpp_clean), "x", ncol(ebpp_clean), "\n")

# ============================================================================
# STEP 2: Create pairs from top variance genes
# ============================================================================

# Get labels
train_labels <- ebpp_clean$Label

# Get numeric gene columns only (exclude Label column)
numeric_cols <- setdiff(colnames(ebpp_clean), c("Label", "label"))

# Calculate variance and select top genes
gene_variances <- apply(ebpp_clean[, ..numeric_cols], 2, var)
top_500_genes <- names(sort(gene_variances, decreasing = TRUE)[1:500])

cat("Top 10 genes by variance:\n")
print(head(top_500_genes, 10))

# ============================================================================
# STEP 3: Create all pairs from top genes
# ============================================================================

# Your create_all_pairs_from_genes function
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

all_pairs <- create_all_pairs_from_genes(ebpp_clean, top_500_genes)

# ============================================================================
# STEP 4: Rank pairs by XGBoost importance
# ============================================================================

rank_pairs_xgboost <- function(pair_data, labels, top_n = 200,
                               nrounds = 100, early_stopping = 10) {
  library(xgboost)

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

  # Train/validation split
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

top_pairs <- rank_pairs_xgboost(
  pair_data = all_pairs,
  labels = train_labels,
  top_n = 100,
  nrounds = 50,
  early_stopping = 5
)

# ============================================================================
# STEP 5: Convert to pair_list format
# ============================================================================

convert_to_pair_list <- function(ranked_pairs) {
  pair_list <- c()
  for (i in 1:nrow(ranked_pairs)) {
    pair_list <- c(pair_list, ranked_pairs$gene1[i], ranked_pairs$gene2[i])
  }
  return(pair_list)
}

new_pair_list <- convert_to_pair_list(top_pairs)

cat("\nPair list created:\n")
cat("Length:", length(new_pair_list), "(", length(new_pair_list)/2, "pairs )\n")
cat("First 10 genes:\n")
print(head(new_pair_list, 10))

# ============================================================================
# STEP 6: Verify all genes exist in cleaned data
# ============================================================================

genes_in_pairs <- unique(new_pair_list)
common_genes <- setdiff(colnames(ebpp_clean), c("Label", "label"))
missing_genes <- setdiff(genes_in_pairs, common_genes)

cat("\nValidation:\n")
cat("Unique genes in pair_list:", length(genes_in_pairs), "\n")
cat("Missing genes:", length(missing_genes), "\n")

if (length(missing_genes) > 0) {
  cat("ERROR: These genes are missing:\n")
  print(missing_genes)
  stop("Fix missing genes before proceeding")
} else {
  cat("✓ All genes in pair_list exist in cleaned datasets!\n")
}

# ============================================================================
# STEP 7: Save the pair_list
# ============================================================================

saveRDS(new_pair_list, '../models/pair_list_normalized.rds')
saveRDS(top_pairs, '../models/top_pairs_ranked.rds')

cat("\n✓ Pair list saved to: ../models/pair_list_normalized.rds\n")
cat("✓ Top pairs saved to: ../models/top_pairs_ranked.rds\n")
cat("\nYou can now use this pair_list for robencla training!\n")

# Quick fix - add SampleID to cleaned datasets
library(data.table)

ebpp_clean <- data.table::fread('../data/formatted_full_L1000/EBpp_pancancer_cleaned.csv.gz')
xena_clean <- data.table::fread('../data/formatted_full_L1000/xena_rsem_tpm_cleaned.csv.gz')

# Add SampleID as first column
ebpp_clean <- cbind(data.table(SampleID = paste0("EBPP_", 1:nrow(ebpp_clean))), ebpp_clean)
xena_clean <- cbind(data.table(SampleID = paste0("Xena_", 1:nrow(xena_clean))), xena_clean)

# Save
data.table::fwrite(ebpp_clean, '../data/formatted_full_L1000/EBpp_pancancer_cleaned.csv.gz')
data.table::fwrite(xena_clean, '../data/formatted_full_L1000/xena_rsem_tpm_cleaned.csv.gz')

cat("✓ SampleID added to both datasets\n")


# Conservative parameters for better C4/C6
conservative_params <- list(
  max_depth = 6,
  eta = 0.3,
  nrounds = 32,
  early_stopping_rounds = 2,
  gamma = 0.5,
  lambda = 2.0,
  alpha = 0.5,
  ensemble_size = 7,
  sample_prop = 0.7,
  feature_prop = 0.7,
  subsample = 0.7
)


# Now use the cleaned pair_list
result <- build_robencla_classifier(
  data_path='../data/formatted_full_L1000/EBpp_pancancer_cleaned.csv.gz',
  test_path='../data/formatted_full_L1000/xena_rsem_tpm_cleaned.csv.gz',
  output_path = '../models/immune_optimized_pairs.rds',
  pair_list = new_pair_list,
  sig_list = NULL,
  param_list = conservative_params,
  data_mode = c("namedpairs"),
  train_fraction = 0.8,
  seed = 412,
  sample_id = "SampleID"  # Specify the sample ID column
)



