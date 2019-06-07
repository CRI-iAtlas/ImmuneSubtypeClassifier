
# Test on GDC data #

library(devtools)
library(readr)

devtools::install_github("Gibbsdavidl/ImmuneSubtypeClassifier", force = T)
reload(pkgload::inst('ImmuneSubtypeClassifier'))

load('~/ebpp_with_subtypes.rda')

reportedScores <- read_tsv('~/Work/PanCancer_Data/five_signature_mclust_ensemble_results.tsv.gz') # in the package data dir
reportedScores <- as.data.frame(reportedScores)
rownames(reportedScores) <- str_replace_all(reportedScores$AliquotBarcode, pattern = '\\.', replacement = '-')

library(ImmuneSubtypeClassifier)

#gdc <- read.table('~/Work/iAtlas/Subtypes/Subtype-Classifier/GDC_Test/gdc_kirc_dat.csv', header=T, stringsAsFactors = F, sep=',')

kirc <- read.table('~/Work/iAtlas/Subtypes/Subtype-Classifier/GDC_Test/TCGA-KIRC.htseq_fpkm-uq.tsv.gz', header=T, stringsAsFactors = F, sep='\t')
colnames(kirc) <- str_replace_all(colnames(kirc), pattern = '\\.', replacement = '-')

library("AnnotationDbi")
library("org.Hs.eg.db")

ensemble <- str_split(kirc$Ensembl_ID, pattern = '\\.')
head(ensemble)
ensemble <- unlist(lapply(ensemble, function(a) a[1]))
symbols <- mapIds(org.Hs.eg.db, keys=ensemble, column="SYMBOL", keytype="ENSEMBL", multiVals="first")

idx <- !duplicated(symbols)
kirc2 <- kirc[idx,-1]
symb2 <- symbols[idx]
symb2[is.na(symb2)] <- "g123"
rownames(kirc2) <- symb2

ids <- intersect(colnames(kirc2), reportedScores$SampleBarcode)

kirc2 <- kirc2[,ids]
rep2 <- reportedScores[reportedScores$SampleBarcode %in% ids,]

all(rep2$SampleBarcode == colnames(kirc2))
#[1] TRUE

xmatSampleBarcode <- str_sub(colnames(Xmat), 1, 16)
xmatSampleBarcode[which(duplicated(xmatSampleBarcode))] <- sapply(xmatSampleBarcode[which(duplicated(xmatSampleBarcode))], function(a) paste0(a,'-d2'))
colnames(Xmat) <- xmatSampleBarcode

z <- Xmat[,ids]

all(colnames(z) == colnames(kirc2))
#[1] TRUE

gids <- intersect(rownames(kirc2), rownames(z))

kirc2 <- kirc2[gids,]
z2 <- z[gids,]

all(rownames(z2) == rownames(kirc2))
#[1] TRUE

plot(x=as.numeric(z2[1000,]), y=as.numeric(kirc2[1000,]))

calls <- callEnsemble(kirc2, geneids = 'symbol')

zcalls <- callEnsemble(z, geneids = 'symbol')


table(Z=zcalls$BestCall, Kirc2=calls$BestCall)
Kirc2
Z   1   2   3   4
1   1   0   0   4
2   0  12   0   0
3   0  86 368  17
4   0   2   0  16
5   0   0   0   2
6   0   4   1   0

table(rep2$ClusterModel1, Kirc2=calls$BestCall)
Kirc2
    1   2   3   4
1   1   2   0   4
2   0  17   2   0
3   0  70 361  13
4   0   6   2  19
5   0   0   0   3
6   0   9   4   0


> table(rep2$ClusterModel1, Z=zcalls$BestCall)
Z
    1   2   3   4   5   6
1   5   1   1   0   0   0
2   0   9  10   0   0   0
3   0   1 443   0   0   0
4   0   1   8  18   0   0
5   0   0   1   0   2   0
6   0   0   8   0   0   5

missed <- which(calls$BestCall == 2 & rep2$ClusterModel1 == 2)

dim(calls)

a1 <- 'TCGA-A3-3306-01A-01R-0864-07'
s1 <- 'TCGA-A3-3306-01A'

s1 <- rep2$SampleBarcode[ missed[2] ]

load('~/Work/iAtlas/Subtypes/Subtype-Classifier/GDC_Test/work.rda')

plot(x=log2(z2[,s1]), y=kirc2[,s1])

zbin <- dataProc(z, ens[[2]], ci=2)
kbin <- dataProc(kirc2, ens[[2]], ci=2)

table(Z=zbin[s1,mods[[2]]$genes], K=kbin[s1,mods[[2]]$genes])

Gene Bins
    K
Z   1   2   3   4
1  11   2   0   0
2   4 244  32   0
3   1  42 127  42
4   0   4  48 169

### Main difference comes from different gene bins ###

> quantile(log2(z[,1]))
0%       25%       50%       75%      100%
-Inf  2.649856  7.680360  9.700893 16.817034

> quantile(kirc2[,1])
0%      25%      50%      75%     100%
0.00000 11.77782 15.77501 17.62235 26.42332

> quantile(log2(scale(z)[,1]), na.rm = T)
0%         25%         50%         75%        100%
-14.5944319  -3.3676726  -1.9180384  -0.5921198   5.1463145

> quantile(log2(scale(kirc2)[,1]), na.rm = T)
0%         25%         50%         75%        100%
-12.5732142  -1.5289511  -0.8521289  -0.3766692   1.0927962



kirc3 <- kirc2[mods[[2]]$genes,]
z3    <- z2[mods[[2]]$genes,]

kirc3 <- kirc2[mods[[2]]$genes,]
dim(kirc3)
#[1] 822 513

z3 <- z[mods[[2]]$genes,]
dim(z3)
# 822

plot(x=log2(z3[,s1]), y=kirc3[,s1])

library(ggplot2)

df <- data.frame(EBPP=log2(z3[,s1]), XENA=kirc3[,s1], EBPPbin=zbin[s1,], XENAbin=kbin[s1,], BothBins=sapply(1:ncol(zbin), function(i) paste0(zbin[s1,i],'_',kbin[s1,i])))

qplot(data=df, x = EBPP, y = XENA, col=as.factor(EBPPbin))

qplot(data=df, x = EBPP, y = XENA, col=as.factor(XENAbin))

qplot(data=df, x = EBPP, y = XENA, col=as.factor(BothBins))

### SO only a couple genes were zeros in the kirc3... not many


