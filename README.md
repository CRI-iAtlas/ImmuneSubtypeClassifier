---
title: "ImmuneSubtypeClassifier"
output:
  html_document: default
  'html_document:': default
---

# ImmuneSubtypeClassifier #
This is an R package for classification of PanCancer immune subtypes. Five gene signatures were used in the initial clustering of tumor samples as part of the Immune Landscape of Cancer manuscript, providing 485 genes that were used to create quartile and binary-valued gene-pair features. With these features, an ensemble of XGBoost classifiers was trained to predict subtype membership, where each member of the ensemble was trained on 70% of 9,129 samples.

```{r}
library(devtools)
install_github("Gibbsdavidl/ImmuneSubtypeClassifier")

library(ImmuneSubtypeClassifier)
```

To get a list of the genes needed:
```{r}
data(ebpp_gene)

head(ebpp_genes_sig)  ### 485 genes are needed
```

To make calls on new data, 

```{r}
Xtest <- as.matrix(X) # has gene IDs in rownames and sample IDs in column names

calls <- callEnsemble(X=Xtest, geneids='symbol')

# or in parallel .. **Not working** #
calls <- parCallEnsemble(X=Xtest, geneids='symbol', numCores=4)

```
Where gene IDs are 'symbol', 'entrez', or 'ensembl'.

If you want to be safe, map your own gene IDs to symbols, 
using your favorite method.

The resulting 'calls' will have 'best calls' in the first column, and probabilities
of belonging to each subtype after that.

* inst/how_the_model_was_fit.Rmd and inst/algorithm_details.txt 
Have information on how the model was built.

* inst/data/five_signature_mclust_ensemble_results.tsv.gz
Contains TCGA subtype membership from the manuscript, suggest using column 'ClusterModel1'.

* inst/important_genes_signature_model.tsv  
A list of the important features in each subtype/ensemble member.

Also see scripts in the test directory for more detailed instructions on
fitting one subtype model, a model per subtype and ensembles of models.
