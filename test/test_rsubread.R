
library(data.table)
library(ImmuneSubtypeClassifier)

rsub <- fread('../data/formatted_full_L1000/Rsubread_tpm.csv.gz')

model <- readRDS('../models/immune_optimized_pairs_64.rds')

pair_list <- model$pair_list

# Apply the same normalization used in training data
normalize_gene_names <- function(gene_names) {
  gene_names <- gsub(" ", "_", gene_names)
  gene_names <- gsub("-", "_", gene_names)
  gene_names <- gsub("\\.", "_", gene_names)  # Convert dots to underscores
  return(gene_names)
}

# Store original names for reference
rsub_original_names <- colnames(rsub)

# Normalize column names
colnames(rsub) <- normalize_gene_names(colnames(rsub))


# Verify the fix
stat1_after <- grep("^STAT1", colnames(rsub), value = TRUE)
cat("STAT1 variants in rsub after normalization:\n")
print(stat1_after)
# Should now show: "STAT1" "STAT1_1"

# Final check
missing_genes <- unique(pair_list[!(pair_list %in% colnames(rsub))])
cat("\nMissing genes after normalization:\n")
print(missing_genes)
cat("Count:", length(missing_genes), "\n")

if (length(missing_genes) == 0) {
  cat("\n✓ All genes in pair_list now exist in rsub!\n")
} else {
  cat("\n✗ Still missing", length(missing_genes), "genes\n")
  cat("Examples:", head(missing_genes, 10), "\n")
}

model$predict(
  data_frame = rsub,
  label_name = 'Label',
  sample_id = 'Barcode'
)

results <- model$results()

conf_matrix <- table(Predicted = results$BestCalls, Actual = model$test_label)

rownames(conf_matrix) <- c("C1", "C2", "C3", "C4", "C5", "C6")
colnames(conf_matrix) <- c("1", "2", "3", "4", "5", "6")

# Calculate per-class metrics
calculate_metrics <- function(conf_mat, class_idx) {
  tp <- conf_mat[class_idx, class_idx]
  fp <- sum(conf_mat[class_idx, ]) - tp
  fn <- sum(conf_mat[, class_idx]) - tp
  tn <- sum(conf_mat) - tp - fp - fn

  sensitivity <- tp / (tp + fn)
  precision <- tp / (tp + fp)
  f1 <- 2 * (precision * sensitivity) / (precision + sensitivity)

  return(c(Sensitivity=sensitivity, Precision=precision, F1=f1))
}

for (i in 1:6) {
  metrics <- calculate_metrics(conf_matrix, i)
  cat(sprintf("C%d: Sens=%.3f, Prec=%.3f, F1=%.3f\n",
              i, metrics[1], metrics[2], metrics[3]))
}

# Overall accuracy
overall_acc <- sum(diag(conf_matrix)) / sum(conf_matrix)
cat(sprintf("\nOverall Accuracy: %.3f\n", overall_acc))
