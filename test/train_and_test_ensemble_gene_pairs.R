
## testing ##

library(devtools)
library(readr)
library(dplyr)
# setwd('~/Work/iAtlas/Subtypes/Subtype-Classifier/')
# using the package

devtools::install_github("Gibbsdavidl/ImmuneSubtypeClassifier", force = T)
reload(pkgload::inst('ImmuneSubtypeClassifier'))
library(ImmuneSubtypeClassifier)

load('/home/davidgibbs/Work/iAtlas/Subtypes/Subtype-Classifier/ensemble_classifier_impt_genes.rda')

# add to a data dir.
reportedScores <- read.table('~/Work/PanCancer_Data/five_signature_mclust_ensemble_results.tsv.gz', sep='\t', header=T, stringsAsFactors = F)
rownames(reportedScores) <- str_replace_all(reportedScores$AliquotBarcode, pattern = '\\.', replacement = '-')

# PanCancer batch corrected expression matrix
ebpp <- read_tsv('~/Work/PanCancer_Data/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv')
#ebpp <- read_table('/home/davidgibbs/Work/iAtlas/Subtypes/Cluster_Work/ebppSubset.tsv.bz2')
##load('/home/davidgibbs/Work/iAtlas/Subtypes/Cluster_Work/ebpp_subset.rda') ### TURNS out many of the important genes not here...

geneList <- str_split(ebpp$gene_id, pattern='\\|')
geneSymbols <- unlist( lapply(geneList, function(a) a[1]) )
# remove duplicate gene names (mostly '?'s)
ddx <- which(duplicated(geneSymbols))


# shared barcodes
bs <- intersect(rownames(reportedScores),colnames(datSubset))

selectGenes <-
allgenes %>%
  group_by(Subtype1, Gene) %>%
  summarise(GainSum = sum(Gain)) %>%
  arrange(desc(GainSum), .by_group=T) %>%
  top_n(n=10)

ebpp2 <- ebpp[geneSymbols %in% c(as.character(selectGenes$Gene), ebpp_subset_genes),]
ebpp2$gene_id <- as.character(geneSymbols[(geneSymbols %in% c(as.character(selectGenes$Gene), ebpp_subset_genes))])
ddx <- which(duplicated(ebpp2$gene_id))
ebpp2 <- ebpp2[-ddx,]
genestokeep <- ebpp2$gene_id
ebpp2 <- as.matrix(ebpp2[,-1])
rownames(ebpp2) <- genestokeep

ebpp <- ebpp2
rm(ebpp2, ddx, ebpp_subset_genes, genestokeep, gi, gidx, values, Y)
gc()


# main matrices
Xmat <- as.matrix(ebpp)
Y <- reportedScores[bs,"ClusterModel1"]

#faster to start from here#
#save(Xmat, Y, geneList, file='~/ebpp_with_subtypes.rda')
load('~/ebpp_with_subtypes.rda')
#load('~/ebpp_subset_genes.rda')
# sample our training and testing groups
idx <- sample(1:ncol(Xmat), size = 0.2 * ncol(Xmat), replace=F)
jdx <- setdiff(1:ncol(Xmat), idx)
Xtrain <- Xmat[,jdx]
Ytrain <- Y[jdx]
Xtest  <- Xmat[,idx]
Ytest <- Y[idx]

# save memory
rm(ebpp, X, Xmat)
gc()

#fitting all models
breakVec=c(0, 0.25, 0.5, 0.75, 1.0)
params=list(max_depth = 5, eta = 0.5, nrounds = 100, nthread = 5, nfold=5)

# list of models
ens <- fitEnsembleModel(Xtrain, Ytrain, n=10, sampSize=0.7, ptail=0.02, params=params, breakVec=breakVec, mtype='pairs')

# calling subtypes on the test set
calls <- callEnsemble(ens, Xtest)

# model performance plots
perfs <- subtypePerf(calls, Ytest)


library(gridExtra)
x <- grid.arrange(perfs[[1]]$plot,perfs[[2]]$plot,perfs[[3]]$plot,perfs[[4]]$plot,perfs[[5]]$plot,perfs[[6]]$plot, ncol=6, nrow=1 )
ggsave(x, file='roc_plot.png')
