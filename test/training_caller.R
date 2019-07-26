library(stringr)
library(ggplot2)

load('~/Work/iAtlas/Subtype-Classifier/datasets/pancancer_ebpp/pancancer_ebpp_calls.rda')
erep2 <- Y
#ecalls
#erep2

load('~/Work/iAtlas/Subtype-Classifier/datasets/xena_kallisto_tpm/xena_kallisto_tpm_results.rda')
#krep2
#kcalls

load('~/Work/iAtlas/Subtype-Classifier/datasets/xena_rsem_fpkm/xena_rsem_fpkm_results.rda')
#fcalls
#frep

load('~/Work/iAtlas/Subtype-Classifier/datasets/xena_rsem_tpm/xena_RSEM_tpm_results.rda')
#rcalls
#rep2

# the calls are not in the same order.

rownames(rcalls) <- rcalls$SampleIDs

shortnames <- str_sub(ecalls$SampleIDs, start=1, end=15)
didx <- which(duplicated(shortnames))
ecalls <- ecalls[-didx,]
rownames(ecalls) <- shortnames[-didx]
erep2 <- erep2[-didx]
names(erep2) <- rownames(ecalls)

shortnames <- str_replace_all(string=fcalls$SampleIDs, pattern = '\\.', replacement ='-')
rownames(fcalls) <- shortnames
rownames(frep) <- str_replace_all(string=frep$ShortBarcode, pattern = '\\.', replacement ='-')

shortnames <- str_replace_all(string=kcalls$SampleIDs, pattern = '\\.', replacement ='-')
rownames(kcalls) <- shortnames
rownames(krep2) <- str_replace_all(string=krep2$ShortSampleBarcode, pattern = '\\.', replacement ='-')

rownames(rep2) <- rep2$ShortBarcode


sidx <- rownames(rcalls)
ex <- ecalls[sidx,]
ey <- erep2[sidx]
fx <- fcalls[sidx,]
fy <- frep[sidx,]
kx <- kcalls[sidx,]
ky <- krep2[sidx,]
rx <- rcalls[sidx,]
ry <- rep2[sidx,]

qplot(x=rx$`1`, y=kx$`1`, xlab='xena rsem tpm', ylab='xena kallisto tpm', main='C1')
qplot(x=rx$`1`, y=fx$`1`, xlab='xena rsem tpm', ylab='xena rsem fpkm', main='C1')
qplot(x=ex$`1`, y=kx$`1`, xlab='PanCancer EBpp', ylab='xena kallisto tpm', main='C1')
qplot(x=ex$`2`, y=fx$`2`, xlab='PanCancer EBpp', ylab='xena rsem fpkm', main='C2')
qplot(x=ex$`3`, y=fx$`3`, xlab='PanCancer EBpp', ylab='xena rsem fpkm', main='C3')
qplot(x=ex$`4`, y=fx$`4`, xlab='PanCancer EBpp', ylab='xena rsem fpkm', main='C4')
qplot(x=fx$`1`, y=kx$`1`, xlab='xena rsem fpkm', ylab='xena kallisto tpm', main='C1', col=fx$BestCall)


# xena RSEMs very similar
# xena kallisto more similar to xena RSEM
# xena RSEM very diff from eb++
# xena kallisto less diff from ebp++


scaller <- xgboost::xgb.cv(data = as.matrix(ecalls[,3:8]), label = (erep2-1),
                           nfold = 5, nrounds = 250, early_stopping_rounds = 2,
                            params = list(objective='multi:softmax', num_class=6))


scaller <- xgboost::xgboost(data = as.matrix(ecalls[,3:8]),
                           label = (erep2-1),
                           nrounds = 8,
                           params = list(objective='multi:softmax', num_class=6))


tmat <- rbind(as.matrix(ex[,3:8]), as.matrix(kx[,3:8]), as.matrix(rx[,3:8]))
tvec <- c(ey, ky$ClusterModel1, ry$ClusterModel1)



xgboost::xgb.cv(data = tmat, label = (tvec-1),
                       nfold = 5, nrounds = 150, early_stopping_rounds = 2,
                       params = list(objective='multi:softmax', num_class=6))


scaller <- xgboost::xgboost(data = tmat,
                            label = (tvec-1),
                            nrounds = 14,
                            params = list(objective='multi:softmax', num_class=6))

# calling on the left out rsem tpm
rnew <- predict(scaller, as.matrix(rx[,3:8]))
tab1 <- table(NEW=rnew, REP=ry$ClusterModel1)

diag(tab1) / colSums(tab1)
#1         2         3         4         5         6
#0.8431752 0.9490261 0.8598972 0.8077945 0.9444444 0.3611111


enew <- predict(scaller, as.matrix(ex[,3:8]))
etab1 <- table(NEW=enew, REP=ey)
diag(etab1) / colSums(etab1)

fnew <- predict(scaller, as.matrix(fx[,3:8]))
ftab1 <- table(NEW=fnew, REP=fy$ClusterModel1)
diag(ftab1) / colSums(ftab1)


