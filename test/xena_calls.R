library(ImmuneSubtypeClassifier)

load('~/xena_test_data.rda')

idx <- sample(x=1:8500, size = 1000, replace = F)
xenasm <- xena2[,idx]
repsm <- rep2[idx,]

rm(xena2, rep2)
gc()

xcalls <- callEnsemble(X=xenasm, geneids = 'ensembl')

xperfs <- subtypePerf(xcalls, Ytest = repsm$ClusterModel1)

xtable1 <- table(T=repsm$ClusterModel1, Xena=xcalls$BestCall)

save(xcalls, xperf, xtable1, file='xenasm_results.rda')

library(ggplot2)
library(gridExtra)
x <- grid.arrange(xperfs[[1]]$plot,xperfs[[2]]$plot,xperfs[[3]]$plot,xperfs[[4]]$plot,xperfs[[5]]$plot,xperfs[[6]]$plot, ncol=6, nrow=1 )
plot(x)
