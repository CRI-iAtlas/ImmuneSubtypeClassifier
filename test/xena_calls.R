library(ImmuneSubtypeClassifier)

load('xena_test_data.rda')

xcalls <- callEnsemble(X=xena2, geneids = 'ensembl')

xperf <- subtypePerf(xcalls, Ytest = rep2$ClusterModel1)

xtable1 <- table(T=rep2$ClusterModel1, Xena=xcalls$BestCall)

save(xcalls, xperf, xtable1, file='xena_results.rda')

