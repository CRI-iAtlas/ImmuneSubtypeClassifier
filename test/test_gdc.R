
# Test on GDC data #

library(devtools)
library(readr)

devtools::install_github("Gibbsdavidl/ImmuneSubtypeClassifier", force = T)
reload(pkgload::inst('ImmuneSubtypeClassifier'))


reportedScores <- read_tsv('~/Work/PanCancer_Data/five_signature_mclust_ensemble_results.tsv.gz') # in the package data dir
reportedScores <- as.data.frame(reportedScores)
rownames(reportedScores) <- str_replace_all(reportedScores$AliquotBarcode, pattern = '\\.', replacement = '-')

library(ImmuneSubtypeClassifier)

gdc <- read.table('~/Work/iAtlas/Subtypes/Subtype-Classifier/GDC_Test/gdc_kirc_dat.csv', header=T, stringsAsFactors = F, sep=',')

