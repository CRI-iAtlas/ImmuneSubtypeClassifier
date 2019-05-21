# ImmuneSubtypeClassifier #
An R package for classification of immune subtypes, in cancer, using gene expression data.

```{r}
# Install devtools from CRAN
install.packages("devtools")

# Or the development version from GitHub:
# install.packages("devtools")
devtools::install_github("Gibbsdavidl/ImmuneSubtypeClassifier")
```

First we read in the PanCancer expression matrix.

```{r}
pathToEbpp <- 'ebpp.rda'
load(pathToEbpp)
```
  
Then we get the reported TCGA immune subtypes.

```{r}
reportedScores <- read.table('five_signature_mclust_ensemble_results.tsv.gz', sep='\t', header=T, stringsAsFactors = F)
rownames(reportedScores) <- str_replace_all(reportedScores$AliquotBarcode, pattern = '\\.', replacement = '-')
```

To get the matrix and phenotypes in the same order:

```{r}
X <- ebpp[, rownames(reportedScores)]
Xmat <- as.matrix(X)
Y <- reportedScores[,"ClusterModel1"]
```

Then to preprocess the data:

```{r}
res0 <- trainDataProc(Xmat, Y, cluster='1')
dat <- res0$dat
testRes <- res0$testRes
```

To fit one model:

```{r}
params <- list(max_depth = 2, eta = 0.5, nrounds = 33, nthread = 5)
mod1 <- fitOneModel(dat$Xbin, dat$Ybin, params)
```

And to see model performance:

```{r}
res1 <- modelPerf(mod1, dat$Xbin, dat$Ybin)
res1$modelError
res1$plot
``` 

To use cross validation to select the number of training rounds (trees):

```{r}
params <- list(max_depth = 3, eta = 0.3, nrounds = 150, nthread = 5, nfold=5)
mod1 <- cvFitOneModel(dat$Xbin, dat$Ybin, params)
```





