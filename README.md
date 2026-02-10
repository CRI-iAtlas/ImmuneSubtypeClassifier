# ImmuneSubtypeClassifier

**ImmuneSubtypeClassifier** is an R package for robust immune subtype classification of tumor samples using gene expression data. It uses the "Robencla" (Robust Ensemble Classifier) framework (github.com/gibbsdavidl/Robencla), utilizing an ensemble of XGBoost models to assign samples to one of six immune subtypes (C1–C6).

## ⚠️ Critical Dependency Note

**This package requires `xgboost` version < 2.0.0.**

Due to major changes in the XGBoost model serialization format, models trained on version 1.x cannot be loaded by version 2.x or 3.x. The pre-trained models included in this package were built using XGBoost 1.7.x.

To ensure compatibility, this package enforces `xgboost (< 2.0.0)`. Using this package with renv is encouraged,
but if you are setting up your environment manually:

```r
# Install a compatible version of xgboost
remotes::install_version("xgboost", version = "1.7.8.1")

```

## Installation

You can install the development version of ImmuneSubtypeClassifier from GitHub:

```r
# Install devtools if you haven't already
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("gibbsdavidl/robencla")
devtools::install_github("gibbsdavidl/ImmuneSubtypeClassifier")

```

## Quick Start: Predicting Subtypes

The function `getFeaturesPairList()` returns a named list, C1-C6, where feature pairs
are defined.  Pairs are (1,2), (3,4), (5,6), etc. for each cluster label.

The function `getFeaturesGeneTable()` returns the gene ID table used with the 
TCGA PanCancer EB++ expression training data (hg19).  Use this function to subset your data with 
columns `Symbol`, `Entrez`, and `Ensembl`.

The main function `callSubtypes()` handles gene matching, data transformation, and prediction.

```r

library(ImmuneSubtypeClassifier)

# Get the list of feature-pair genes used for the model.
# ...there's model_genes_list$model_genes and model_genes_list$gene_map
model_genes_list <- modelGenes(model_path='../models/immune_optimized_99_pairs.rds')

# Confirm the gene_map contains all genes.
length(model_genes$model_genes) == nrow(model_genes$gene_map)

# Call the subtypes, can also pass in a data.frame / tibble etc.
result <- callSubtypes(X_or_path = '../data/gene_expression_rsem_tpm.csv.gz',
                     model = NULL,  # genes in columns and samples in rows.
                     model_path = '../models/immune_optimized_99_pairs.rds',
                     geneid =   "symbol",  # how are the gene IDs encoded
                     sampleid = 'Barcode', # column name with sample IDs
                     labelid=  'Label')    # column name with sample Labels (optional)

# get the trained robencla model.
model <- result$Model

# get the results table, $BestCall is the predicted subtype
pred  <- result$Pred

# confusion matrix
table(pred$BestCall, pred$Label)

# get class prediction metrics.
model$classification_metrics(labels = pred$Label, calls = pred$BestCall)


```

### Output Format

The output is a data frame containing:

* **SampleIDs**: The sample identifiers.
* **BestCall**: The predicted immune subtype (1–6).
* **1–6**: Probability scores for each subtype.

| SampleIDs | BestCall | 1 | 2 | 3 | 4 | 5 | 6 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Sample_01 | 3 | 0.02 | 0.10 | **0.85** | 0.01 | 0.01 | 0.01 |
| Sample_02 | 1 | **0.78** | 0.05 | 0.10 | 0.02 | 0.04 | 0.01 |

## Advanced: Retraining the Model

If you have a labeled training set and wish to rebuild the model (e.g., after updating XGBoost or changing feature pairs), you can use `build_robencla_classifier`.  I can also send you the PanCancer training set, already formatted.

### Training Data Requirements

* **CSV file** with samples as rows and genes as columns.
* Must contain a **Label** column (e.g., "C1", "C2") and a **Sample ID** column.

```r

library(ImmuneSubtypeClassifier)

# Conservative parameters, see robencla & xgboost for parameters 
conservative_params <- list(
  max_depth = 8,
  eta = 0.2,
  nrounds = 64,
  early_stopping_rounds = 4,
  gamma = 0.3,
  lambda = 1.9,
  alpha = 0.3,
  ensemble_size = 11,
  sample_prop = 0.8,
  feature_prop = 0.8,
  subsample = 0.8
)

# then load up a list of gene pairs...
pair_list <- readRDS('../models/pair_list_stratified.rds')
# if there's a gene you want to try removing
this_set <- editPairList(pair_list, 'IGJ')

# Now use the cleaned pair_list
result <- build_robencla_classifier(
  data_path='../data/training_expression_data.csv.gz',
  test_path='../data/testing_expression_data.csv.gz',
  output_path = '../models/immune_optimized_99_pairs.rds',
  pair_list = this_set,
  sig_list = NULL,
  param_list = conservative_params,
  data_mode = c("namedpairs"),
  train_fraction = NULL,  # used if test_path is NULL
  seed = 412,
  sample_id = "Barcode"  # Specify the sample ID column
)



```

## Immune Subtypes

The classifier assigns one of six immune subtypes:

* **C1 (Wound Healing):** High proliferation, angiogenic gene expression.
* **C2 (IFN-gamma Dominant):** High M1/M2 macrophage polarization, strong CD8 signal.
* **C3 (Inflammatory):** Elevated Th17 and Th1 genes, low tumor cell proliferation.
* **C4 (Lymphocyte Depleted):** Macrophage sequestration, Th2 shift.
* **C5 (Immunologically Quiet):** Low lymphocyte and macrophage responses.
* **C6 (TGF-beta Dominant):** High TGF-beta signature, lymphocytic infiltration.

## Troubleshooting

### Error: `xgboost.Booster object is corrupted`

This means you are trying to load the model with a newer version of XGBoost (>= 2.0) than was used to train it.
**Fix:** Downgrade XGBoost to 1.7.8.1 or re-train the model using `build_robencla_classifier`.

## License

This project is licensed under the terms of the MIT license.
![License: MIT](https://img.shields.io)

```
