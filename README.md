---
title: "ImmuneSubtypeClassifier"
output:
  html_document: default
  'html_document:': default
---

# ImmuneSubtypeClassifier #
An R package for classification of immune subtypes, in cancer, using gene expression data.

Currently, building the classifier... instuctions on using the classifier coming soon!



```{r}
# Install devtools from CRAN
install.packages("devtools")

# Or the development version from GitHub:
# install.packages("devtools")
devtools::install_github("Gibbsdavidl/ImmuneSubtypeClassifier")

library(utils)
library(readr)
library(stringr)
```

First we read in the PanCancer expression matrix,
'EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv'.

```{r}
download_file(src='http://api.gdc.cancer.gov/data/3586c0da-64d0-4b74-a449-5ff4d9136611', dest='~/ebpp.tsv')
ebpp <- read_tsv('~/ebpp.tsv')
```
  
Then we get the reported TCGA immune subtypes.

```{r}
reportedScores <- read_tsv('data/five_signature_mclust_ensemble_results.tsv.gz') # in the package data dir
reportedScores <- as.data.frame(reportedScores)
rownames(reportedScores) <- str_replace_all(reportedScores$AliquotBarcode, pattern = '\\.', replacement = '-')
```

A little processing on the EB++ gene identifiers
```{r}
geneList <- str_split(ebpp$gene_id, pattern='\\|')
geneSymbols <- unlist( lapply(geneList, function(a) a[2]) )
```

To get the matrix and phenotypes in the same order:

```{r}
X <- ebpp[, rownames(reportedScores)]
rownames(X) <- geneSymbols
Xmat <- as.matrix(X)
Y <- reportedScores[,"ClusterModel1"]
```

Then to preprocess and filter the data:

```{r}
library(ImmuneSubtypeClassifier)
res0     <- trainDataProc(Xmat, Y, cluster='1')
testRes  <- res0$testRes   # scores used in filtering
breakVec <- res0$breakVec  # break points used in binning
dat      <- res0$dat
```

To fit one model:

```{r}
params <- list(max_depth = 2, eta = 0.5, nrounds = 33, nthread = 5)
mod1 <- fitOneModel(dat$Xbin, dat$Ybin, params, breakVec)
```

To use cross validation to select the number of training rounds (trees):

```{r}
params <- list(max_depth = 3, eta = 0.3, nrounds = 150, nthread = 5, nfold=5)
mod1 <- cvFitOneModel(dat$Xbin, dat$Ybin, params, breakVec)
```

And to see model performance:

```{r}
res1 <- modelPerf(mod1$bst, dat$Xbin, dat$Ybin)
res1$modelError
res1$plot
``` 

To fit the list of subtype models (one for each subtype):

```{r}
xbgParams <- list(max_depth = 5, eta = 0.5, nrounds = 150, nthread = 5, nfold=5)
mods <- fitSubtypeModel(Xmat, Y, xbgParams, breakVec)
```

Calling subtypes on new data (genes in rows, samples in columns):

```{r}
Xtest <- read_tsv('path_to_my_data.tsv')
calls <- callSubtypes(mods, Xtest)
```

Each subtype (C1-C6) gives a probability of belonging to that class.
A sample, then, can take the subtype with the highest probability.

And to see model performance:

```{r}
Ytest <- read_tsv('path_to_my_data.tsv')
res1 <- modelPerf2(calls, Ytest, '1')
res1$modelError
res1$plot
``` 
