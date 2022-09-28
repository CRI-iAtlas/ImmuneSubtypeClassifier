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
install_github("CRI-iAtlas/ImmuneSubtypeClassifier")

# Right Now, the newest version of xgboost is incompatible. 
# Please use a prior version, 1.0.0.1 or 1.0.0.2 should work.
devtools::install_version("xgboost", version = "1.0.0.1")

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

But, to see where and if gene ID matches have failed:

```{r}

calls <- geneMatchErrorReport(X=Xtest, geneid='symbol')

```
This returns the proportion of missing genes (from 485 total)
and a data.frame of missing gene IDs.

```
$matchError
[1] 0.03505155

$missingGenes
        Symbol Entrez         Ensembl
1794  C12orf24  29902 ENSG00000204856
1841  C13orf18  80183 ENSG00000102445
1844  C13orf27  93081 ENSG00000151287
```

The resulting 'calls' will have 'best calls' in the first column, and probabilities
of belonging to each subtype after that.

* inst/how_the_model_was_fit.Rmd and inst/algorithm_details.txt 
Have information on how the model was built.

* inst/data/five_signature_mclust_ensemble_results.tsv.gz
Contains TCGA subtype membership from the manuscript, suggest using column 'ClusterModel1'.

* inst/important_features_in_the_ensemble_model.tsv  
A list of the important features in each subtype/ensemble member.

Also see scripts in the test directory for more detailed instructions on
fitting one subtype model, a model per subtype and ensembles of models.


This following script should work.
```{r}
library(readr)
library(ImmuneSubtypeClassifier)

download.file(url = 'https://raw.githubusercontent.com/CRI-iAtlas/shiny-iatlas/develop/data/ebpp_test1_1to20.tsv', destfile = 'ebpp_test.tsv')
dat <- read_tsv('ebpp_test.tsv')

dat2 <- as.data.frame(dat[!duplicated(dat$GeneID),])
Xmat <- dat2[,-1]
rownames(Xmat) <- dat2[,1]

res0 <- ImmuneSubtypeClassifier::callEnsemble(X = Xmat, geneids = 'symbol')
res0

   SampleIDs BestCall            1            2            3            4            5            6
1        XY1        4 1.121173e-07 1.873900e-06 6.311234e-02 0.9027497470 5.132385e-02 5.121503e-05
2        XY2        3 4.243225e-02 2.309184e-06 5.495564e-01 0.0084770960 1.000665e-04 2.261875e-04
3        XY3        4 3.095071e-04 1.029907e-06 1.635270e-01 0.9063920975 5.731729e-03 1.741630e-04
4        XY4        4 1.154656e-04 3.823049e-07 2.888787e-02 0.8831390142 1.236814e-03 7.847964e-05
5        XY5        4 1.193054e-07 8.741189e-06 5.060996e-01 0.9260828793 1.232723e-02 5.358944e-04
6        XY6        4 2.479082e-04 5.528710e-05 1.183165e-04 0.9923758209 1.019272e-03 2.854635e-04
7        XY7        6 9.421768e-03 1.079501e-03 1.000559e-01 0.0084094275 3.649302e-05 6.113418e-01
8        XY8        3 6.143126e-04 2.594000e-06 8.872822e-01 0.1306741871 8.199657e-04 1.142911e-04
9        XY9        4 1.001003e-04 4.479314e-06 3.989228e-03 0.9793768525 2.182218e-03 1.554600e-04
10      XY10        2 4.624279e-06 9.888145e-01 7.359737e-06 0.0053910189 3.058250e-05 6.004089e-04
11      XY11        4 5.837736e-05 7.877929e-06 1.970483e-03 0.9895575941 4.063185e-03 6.252101e-04
12      XY12        3 1.944198e-06 2.005060e-06 3.616153e-01 0.5258071870 1.320550e-02 2.124778e-04
13      XY13        4 9.002434e-07 9.563100e-04 2.673559e-06 0.9931332767 1.996231e-04 1.273039e-04
14      XY14        4 1.235332e-05 3.907399e-05 2.785671e-03 0.9972651005 2.466791e-03 5.618507e-04
15      XY15        4 4.831095e-05 1.127145e-04 1.685999e-04 0.9935480356 4.465038e-03 2.503012e-04
16      XY16        3 2.731656e-03 1.722274e-05 9.742068e-01 0.0007705171 4.555504e-05 6.263918e-03
17      XY17        4 3.431090e-06 1.327718e-03 1.570120e-05 0.9971910715 3.616062e-05 1.134532e-04
18      XY18        3 3.816116e-07 5.000733e-03 8.685111e-01 0.0030849737 2.343850e-04 2.590704e-03
19      XY19        4 4.787456e-05 2.856209e-06 5.303099e-01 0.9858489633 2.194258e-02 9.931145e-05
20      XY20        4 9.208488e-04 1.361713e-04 3.989896e-04 0.9504880905 1.666186e-03 1.220569e-02
```

These results match what's found on cri-iatlas.org / tools.



In looking at feature importance: 
You will see that really important features for classification are based on doing the binary gene-gene comparison, but on a signature level.  It summarizes the question "are the genes in signature 1 (s1) expressed at a lower level than signature 2 (s2)?"  In short "s1s2".  

```
label signature_name
s1	  LIexpression_score	
s2	  CSF1_response	
s3	  Module3_IFN_score	
s4	  TGFB_score_21050467	
s5	  CHANG_CORE_SERUM_RESPONSE_UP	
```

