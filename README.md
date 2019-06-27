---
title: "ImmuneSubtypeClassifier"
output:
  html_document: default
  'html_document:': default
---

# ImmuneSubtypeClassifier #
An R package for classification of immune subtypes, in cancer, using binned and binary gene-pairs expression data.

```{r}
library(devtools)
install_github("Gibbsdavidl/ImmuneSubtypeClassifier")

library(ImmuneSubtypeClassifier)
```

To get a list of the genes needed:
```{r}
data(ebpp_genes)

head(ebpp_genes_sig)  ### 485 genes are needed
```

To make calls on new data, 

```{r}
Xtest <- as.matrix(X) # has gene IDs in rownames and sample IDs in column names

calls <- callEnsemble(X=Xtest, geneids='symbol')
```
Where gene IDs are 'symbol', 'entrez', or 'ensembl' .

The resulting 'calls' will have 'best calls' in the first column, and probabilities
of belonging to each subtype after that.

See inst/how_the_model_was_fit.Rmd and inst/algorithm_details.txt for more information.

See inst/important_genes_signature_model.tsv  
for a list of the important features in each subtype/ensemble member.

Also see scripts in the test directory for more detailed instructions on
fitting one subtype model, a model per subtype and ensembles of models.
