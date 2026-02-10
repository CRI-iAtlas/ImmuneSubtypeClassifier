

library(ImmuneSubtypeClassifier)

# get the list of genes used for the model
# ...there's model_genes_list$model_genes and model_genes_list$gene_map
model_genes_list <- modelGenes(model_path='../models/immune_optimized_99_pairs.rds')

# confirm the gene_map contains all genes
length(model_genes$model_genes) == nrow(model_genes$gene_map)

# call the subtypes, can also pass in a data.frame / tibble etc.
res0 <- callSubtypes('../data/formatted_full_L1000/xena_rsem_tpm.csv.gz',
                     model = NULL,  # genes in columns and samples in rows.
                     model_path = '../models/immune_optimized_99_pairs.rds',
                     geneid =   "symbol",  # how are the gene IDs encoded
                     sampleid = 'Barcode', # column name with sample IDs
                     labelid=  'Label')    # column name with sample Labels (optional)

# get the trained robencla model.
model <- res0$Model

# get the results table, $BestCall is the predicted subtype
pred  <- res0$Pred

# confusion matrix
table(pred$BestCall, pred$Label)

# get class prediction metrics.
model$classification_metrics(labels = pred$Label, calls = pred$BestCall)


