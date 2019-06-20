
# Test on GDC data #

library(devtools)
library(readr)
library(stringr)

devtools::install_github("Gibbsdavidl/ImmuneSubtypeClassifier", force = T)
reload(pkgload::inst('ImmuneSubtypeClassifier'))
library(ImmuneSubtypeClassifier)

#load('~/ebpp_with_subtypes.rda')

reportedScores <- read_tsv('~/Work/PanCancer_Data/five_signature_mclust_ensemble_results.tsv.gz') # in the package data dir
reportedScores <- as.data.frame(reportedScores)
rownames(reportedScores) <- str_replace_all(reportedScores$AliquotBarcode, pattern = '\\.', replacement = '-')


#gdc <- read.table('~/Work/iAtlas/Subtypes/Subtype-Classifier/GDC_Test/gdc_kirc_dat.csv', header=T, stringsAsFactors = F, sep=',')

xena <- read.table('~/xena_tcga_RSEM_gene_tpm.gz', header=T, stringsAsFactors = F, sep='\t')
colnames(xena) <- str_replace_all(colnames(xena), pattern = '\\.', replacement = '-')

#library("AnnotationDbi")
#library("org.Hs.eg.db")

#ensemble <- str_split(kirc$Ensembl_ID, pattern = '\\.')
#head(ensemble)
#ensemble <- unlist(lapply(ensemble, function(a) a[1]))
#symbols <- mapIds(org.Hs.eg.db, keys=ensemble, column="SYMBOL", keytype="ENSEMBL", multiVals="first")

#idx <- !duplicated(symbols)
#kirc2 <- kirc[idx,-1]
#symb2 <- symbols[idx]
#symb2[is.na(symb2)] <- "g123"
#rownames(kirc2) <- symb2


ids <- intersect(colnames(xena), str_sub(string = reportedScores$SampleBarcode, start = 1, end = 15))

xena2 <- xena[,ids]
rep2 <- reportedScores[reportedScores$SampleBarcode %in% ids,]

all(rep2$SampleBarcode == colnames(xena2))
#[1] TRUE

save(xena2, rep2, file='~/xena_test_data.rda')

#######################################################################################33

# restart R

devtools::install_github("Gibbsdavidl/ImmuneSubtypeClassifier", force = T)
reload(pkgload::inst('ImmuneSubtypeClassifier'))
library(ImmuneSubtypeClassifier)

load('~/xena_test_data.rda')

xcalls <- callEnsemble(X=xena2, geneids = 'ensembl')

subtypePerf(xcalls, Ytest = rep2$ClusterModel1)

table(T=rep2$ClusterModel1, Xena=xcalls$BestCall)

