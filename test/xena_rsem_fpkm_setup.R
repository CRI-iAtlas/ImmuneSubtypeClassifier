# xena rsem fpkm setup


reportedScores <- read.table('~/Work/PanCancer_Data/five_signature_mclust_ensemble_results.tsv.gz', sep='\t', header=T, stringsAsFactors = F)
rownames(reportedScores) <- str_replace_all(reportedScores$AliquotBarcode, pattern = '\\.', replacement = '-')
reportedScores$SampleBarcodeShort <- str_sub(reportedScores$SampleBarcode, 1, 15)

dat <- read.table('tcga_RSEM_gene_fpkm.gz', stringsAsFactors=F, sep='\t', header=T)

load('~/Code/iatlas/ImmuneSubtypeClassifier/data/ebpp_gene.rda')

library(stringr)
 
reportedScores <- read.table('~/Work/PanCancer_Data/five_signature_mclust_ensemble_results.tsv.gz', sep='\t', header=T, stringsAsFactors = F)
rownames(reportedScores) <- str_replace_all(reportedScores$AliquotBarcode, pattern = '\\.', replacement = '-')
reportedScores$SampleBarcodeShort <- str_sub(reportedScores$SampleBarcode, 1, 15)

genes <- str_split(dat$sample, pattern='\\.')
genes2 <- unlist(lapply(genes, function(a) a[1]))

dat2 <- dat[genes2 %in% ebpp_genes_sig$Ensembl,]

reportedScores$ShortBarcode <- str_sub(reportedScores$AliquotBarcode, start=1, end=15)

ids <- intersect(reportedScores$ShortBarcode, colnames(dat2))

dat3 <- dat2[,ids]

rep1 <- reportedScores[reportedScores$ShortBarcode %in% ids,]

rep2 <- rep2[colnames(dat3),]

rownames(dat3) <- dat2$sample

fdat <- dat3
frep <- rep2

all(colnames(fdat) == frep$ShortBarcode)

