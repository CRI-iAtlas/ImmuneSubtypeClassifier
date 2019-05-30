
## testing ##

library(devtools)
library(readr)
# setwd('~/Work/iAtlas/Subtypes/Subtype-Classifier/')
# using the package

devtools::install_github("Gibbsdavidl/ImmuneSubtypeClassifier", force = T)
#reload(pkgload::inst('ImmuneSubtypeClassifier'))

library(ImmuneSubtypeClassifier)

# PanCancer batch corrected expression matrix
ebpp <- read_tsv('EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv')

geneList <- str_split(ebpp$gene_id, pattern='\\|')
geneSymbols <- unlist( lapply(geneList, function(a) a[1]) )

# add to a data dir.
reportedScores <- read.table('five_signature_mclust_ensemble_results.tsv.gz', sep='\t', header=T, stringsAsFactors = F)
rownames(reportedScores) <- str_replace_all(reportedScores$AliquotBarcode, pattern = '\\.', replacement = '-')

# shared barcodes
bs <- intersect(rownames(reportedScores),colnames(ebpp))

ddx <- which(duplicated(geneSymbols))

X <- as.data.frame(ebpp[-ddx, rownames(reportedScores)])
rownames(X) <- geneSymbols[-ddx]
Xmat <- as.matrix(X)
Y <- reportedScores[,"ClusterModel1"]

# save memory
# sample our training and testing groups
idx <- sample(1:ncol(X), size = 0.4 * ncol(X), replace=F)
jdx <- setdiff(1:ncol(X), idx)

Xtrain <- X[,jdx]
Ytrain <- Y[jdx]
Xtest  <- X[,idx]
Ytest <- Y[idx]

rm(ebpp, X, Xmat)
gc()

#fitting all models
# first process training data
#breakVec=c(0, 0.25, 0.5, 0.75, 1.0)
breakVec=c(0, 0.33, 0.66, 1.0)
mods <- fitSubtypeModel(Xtrain,Ytrain, ptail = 0.02, params=list(max_depth = 5, eta = 0.5, nrounds = 100, nthread = 5, nfold=5), breakVec = breakVec)

# calling subtypes on the test set
calls <- callSubtypes(mods, Xtest)

save(idx, jdx, mods, calls, Ytest, file='results.rda')


