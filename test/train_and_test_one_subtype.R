

##### Testing 1 subtype classifier. #####

library(devtools)
library(readr)

setwd('~/Code/iatlas/ImmuneSubtypeClassifier/')
# using the package

ebpp <- read_tsv('~/Work/PanCancer_Data/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv')

reportedScores <- read_tsv('data/five_signature_mclust_ensemble_results.tsv.gz') # in the package data dir
reportedScores <- as.data.frame(reportedScores)
rownames(reportedScores) <- str_replace_all(reportedScores$AliquotBarcode, pattern = '\\.', replacement = '-')

geneList <- str_split(ebpp$gene_id, pattern='\\|')
geneSymbols <- unlist( lapply(geneList, function(a) a[1]) )
ddx <- which(duplicated(geneSymbols))


X <- as.data.frame(ebpp[, rownames(reportedScores)])
X <- X[-ddx,]
rownames(X) <- geneSymbols[-ddx]
Xmat <- as.matrix(X)
Y <- reportedScores[,"ClusterModel1"]

devtools::install_github("Gibbsdavidl/ImmuneSubtypeClassifier", force = T)
reload(pkgload::inst('ImmuneSubtypeClassifier'))

library(ImmuneSubtypeClassifier)


idx <- sample(1:ncol(X), size = 0.1 * ncol(X), replace=F)
jdx <- setdiff(1:ncol(X), idx)
Xtrain <- Xmat[,jdx]
Ytrain <- Y[jdx]
Xtest  <- Xmat[,idx]
Ytest  <- Y[idx]

################################3

#save(Xtest, Ytest, Ytrain, Xtrain, idx, jdx, file='testdata.rda')
load('testdata.rda')

cluster <- '4'

# first process training data
breakVec=c(0, 0.33, 0.66, 1.0)
res0 <- trainDataProc(Xtrain, Ytrain, testRes=NULL, cores=2, cluster=cluster, ptail=0.02, breakVec)
traindat <- res0$dat

# then fit a model
params <- list(max_depth = 5, eta = 0.5, nrounds = 150, nthread = 5, nfold=5)

mod1 <- cvFitOneModel(traindat$Xbin, traindat$Ybin, params, breakVec, traindat$Genes)

# then process the test data
testXbin <- dataProc(Xtest, mod1)
testYbin <- sapply(Ytest, function(yi) if (yi == cluster){1} else {0})

# then model performance
perf1 <- modelPerf(mod1$bst, testXbin, testYbin, cluster)

# ROC plot
perf1$plot


