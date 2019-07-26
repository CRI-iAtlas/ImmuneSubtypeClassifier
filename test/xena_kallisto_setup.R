

library(readr)
library(stringr)
library(dplyr)

xena <- read_tsv('tcga_Kallisto_tpm.gz')


reportedScores <- read.table('~/Work/PanCancer_Data/five_signature_mclust_ensemble_results.tsv.gz', sep='\t', header=T, stringsAsFactors = F)
rownames(reportedScores) <- str_replace_all(reportedScores$AliquotBarcode, pattern = '\\.', replacement = '-')
reportedScores$SampleBarcodeShort <- str_sub(reportedScores$SampleBarcode, 1, 15)

kdx <- which(colnames(xena) %in% reportedScores$SampleBarcodeShort)
xena <- xena[, c(1,kdx)]

#t2g <- read_tsv('Utils_transcripts_to_genes.txt', col_names=F)
t2g <- read_tsv('gencode.v23.annotation.transcript.probemap', col_names=T)


idx <- sample(1:ncol(xena), size=1000, replace=F)

xena2 <- xena[,c(1,idx)]

xena2 <- inner_join(x=t2g, y=xena2, by=c('X1'='sample')) 

xena3 <- xena2[,-c(1,3)]

xena4 <- xena3 %>% group_by(X2) %>% summarise_all(funs(sum))

jdx <- match(x=colnames(xena4), table=reportedScores$SampleBarcodeShort)

all(colnames(xena4) == reportedScores$SampleBarcodeShort[jdx], na.rm=T)
#[1] TRUE

rep2 <- reportedScores[jdx,]  # first row is matched to 'X2' in xena4

all(rep4$SampleBarcodeShort == colnames(xena4)[-1])
#[1] TRUE

Xmat <- as.matrix(xena4 %>% select(-X2))
rownames(Xmat) <- xena4$X2

save(Xmat, Y, file='xena_kallisto_tpm_1k.rda')


