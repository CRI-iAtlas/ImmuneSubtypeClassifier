---
title: "ImmuneSubtypeClassifier"
output:
  html_document: default
  'html_document:': default
---

# ImmuneSubtypeClassifier #
An R package for classification of immune subtypes, in cancer, using binned gene expression data.

```{r}
# Install devtools from CRAN
install.packages("devtools")

# Or the development version from GitHub:
# install.packages("devtools")
devtools::install_github("Gibbsdavidl/ImmuneSubtypeClassifier")

library(ImmuneSubtypeClassifier)
library(readr)
```

To make calls on new data, with genes IDs in the first column and samples names after that.

```{r}
Xnew <- read_tsv('new_data.tsv')

calls <- callEnsemble(ens, Xtest)
```

The resulting 'calls' will have best calls in the first column, and probabilities
of belonging to each subtype after that.

See how_the_model_was_fit.Rmd for details on how the model was fit.

See scripts in the test directory for more detailed instructions on
fitting one subtype model, a model per subtype and ensembles of models.
