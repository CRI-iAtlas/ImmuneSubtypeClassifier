library(data.table)
library(readr)

# Load ORIGINAL data with Barcodes
xena_orig <- read_csv('../data/formatted_full_L1000/xena_rsem_tpm.csv.gz')
ebpp_orig <- read_csv('../data/formatted_full_L1000/EBpp_pancancer.csv.gz')

cat("Original dimensions:\n")
cat("EBPP:", nrow(ebpp_orig), "x", ncol(ebpp_orig), "\n")
cat("Xena:", nrow(xena_orig), "x", ncol(xena_orig), "\n\n")

# Find shared barcodes
shared_barcodes <- intersect(xena_orig$Barcode, ebpp_orig$Barcode)
cat("Shared barcodes:", length(shared_barcodes), "\n")

# Find shared genes (normalize names first)
normalize_gene_names <- function(gene_names) {
  gene_names <- gsub(" ", "_", gene_names)
  gene_names <- gsub("-", "_", gene_names)
  gene_names <- gsub("\\.", "_", gene_names)
  return(gene_names)
}

colnames(xena_orig) <- normalize_gene_names(colnames(xena_orig))
colnames(ebpp_orig) <- normalize_gene_names(colnames(ebpp_orig))

shared_genes <- intersect(colnames(xena_orig), colnames(ebpp_orig))
shared_genes <- setdiff(shared_genes, c("Barcode", "Label", "SampleID", "sample"))

cat("Shared genes:", length(shared_genes), "\n\n")

# # Handle IGJ -> IGJP1 mapping if needed
# if ("IGJP1" %in% colnames(xena_orig) && !("IGJ" %in% colnames(xena_orig)) && "IGJ" %in% colnames(ebpp_orig)) {
#   colnames(xena_orig)[colnames(xena_orig) == "IGJP1"] <- "IGJ"
#   shared_genes <- c(shared_genes, "IGJ")
#   cat("Renamed IGJP1 to IGJ in Xena\n")
# }

# Subset to shared barcodes and genes, PRESERVING ORDER BY BARCODE
xena_matched <- as.data.table(xena_orig[xena_orig$Barcode %in% shared_barcodes, ])
ebpp_matched <- as.data.table(ebpp_orig[ebpp_orig$Barcode %in% shared_barcodes, ])

# Sort both by Barcode
setkey(xena_matched, Barcode)
setkey(ebpp_matched, Barcode)

# Verify alignment
cat("Barcodes aligned:", all(xena_matched$Barcode == ebpp_matched$Barcode), "\n\n")

# Keep only shared genes + Barcode + Label
keep_cols <- c("Barcode", "Label", shared_genes)
keep_cols <- intersect(keep_cols, colnames(xena_matched))

xena_clean_matched <- xena_matched[, ..keep_cols]
ebpp_clean_matched <- ebpp_matched[, ..keep_cols]

# Add sequential SampleID (but keep Barcode for verification)
xena_clean_matched$SampleID <- paste0("Sample_", 1:nrow(xena_clean_matched))
ebpp_clean_matched$SampleID <- paste0("Sample_", 1:nrow(ebpp_clean_matched))

# Reorder columns
final_cols <- c("SampleID", "Barcode", "Label", shared_genes)
final_cols <- intersect(final_cols, colnames(xena_clean_matched))

xena_clean_matched <- xena_clean_matched[, ..final_cols]
ebpp_clean_matched <- ebpp_clean_matched[, ..final_cols]

cat("Final matched dimensions:\n")
cat("EBPP:", nrow(ebpp_clean_matched), "x", ncol(ebpp_clean_matched), "\n")
cat("Xena:", nrow(xena_clean_matched), "x", ncol(xena_clean_matched), "\n\n")

# Log-transform EBPP gene columns ONLY
# Xena pipeline: log2(TPM + 0.001)
cat("\nApplying log2(TPM + 0.001) transformation to EBPP...\n")
for (gene in shared_genes) {
  ebpp_clean_matched[[gene]] <- log2(ebpp_clean_matched[[gene]] + 0.001)
}

# Save properly matched datasets
fwrite(ebpp_clean_matched, '../data/formatted_full_L1000/EBpp_pancancer_matched.csv.gz')
fwrite(xena_clean_matched, '../data/formatted_full_L1000/xena_rsem_tpm_matched.csv.gz')

cat("✓ Saved properly matched datasets with Barcode alignment preserved!\n")
