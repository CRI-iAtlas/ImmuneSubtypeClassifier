---
title: "ImmuneSubtypeClassifier"
output:
  html_document: default
  'html_document:': default
---

# ImmuneSubtypeClassifier #
An R package for classification of immune subtypes, in cancer, using binned and binary gene-pairs expression data.

```{r}
# Install devtools from CRAN
install.packages("devtools")

# Or the development version from GitHub:
# install.packages("devtools")
devtools::install_github("Gibbsdavidl/ImmuneSubtypeClassifier")

library(ImmuneSubtypeClassifier)
```

To make calls on new data, 

```{r}
Xtest <- as.matrix(X) # has gene IDs in rownames and sample IDs in column names

calls <- callEnsemble(ens, Xtest)
```
Where gene IDs are Symbols, Entrez, or Ensemble IDs.

The resulting 'calls' will have 'best calls' in the first column, and probabilities
of belonging to each subtype after that.

See inst/how_the_model_was_fit.Rmd and inst/algorithm_details.txt for more information.

Also see scripts in the test directory for more detailed instructions on
fitting one subtype model, a model per subtype and ensembles of models.
