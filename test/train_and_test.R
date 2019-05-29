
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
ebpp <- as.data.frame(ebpp)
rownames(ebpp) <- ebpp[,1]
ebpp <- ebpp[,-1]

# add to a data dir.
reportedScores <- read.table('five_signature_mclust_ensemble_results.tsv.gz', sep='\t', header=T, stringsAsFactors = F)
rownames(reportedScores) <- str_replace_all(reportedScores$AliquotBarcode, pattern = '\\.', replacement = '-')

# shared barcodes
bs <- intersect(rownames(reportedScores),colnames(ebpp))

# the main data matrices
X <- ebpp[, bs]
Y <- reportedScores[bs,"ClusterModel1"]

# save memory
rm(ebpp)
gc()

# sample our training and testing groups
idx <- sample(1:ncol(X), size = 0.4 * ncol(X), replace=F)
jdx <- setdiff(1:ncol(X), idx)

Xtrain <- X[,jdx]
Ytrain <- Y[jdx]
Xtest  <- X[,idx]
Ytest <- Y[idx]

#fitting all models
mods <- fitSubtypeModel(Xtrain,Ytrain,params=params)

# calling subtypes on the test set
calls <- callSubtypes(mods, Xtest)


