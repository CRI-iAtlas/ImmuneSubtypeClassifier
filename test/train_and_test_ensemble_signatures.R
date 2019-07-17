
## testing ##

library(devtools)
library(readr)
# setwd('~/Work/iAtlas/Subtypes/Subtype-Classifier/')
# using the package

#devtools::install_github("Gibbsdavidl/ImmuneSubtypeClassifier", force = T)
devtools::install_github("Gibbsdavidl/ImmuneSubtypeClassifier")

reload(pkgload::inst('ImmuneSubtypeClassifier'))
library(ImmuneSubtypeClassifier)

load('~/ebpp_signature_genes.rda')

#fitting all models
breakVec=c(0, 0.25, 0.5, 0.75, 1.0)
params=list(max_depth = 5, eta = 0.3, nrounds = 100, nthread = 5, nfold=5)

# list of models
ens <- fitEnsembleModel(Xmat, Y, n=10, sampSize=0.7, ptail=0.05, params=params, breakVec=breakVec)

save(ens, file='~/ensemble_model.rda')

# calling subtypes on the test set
calls <- callEnsemble(ens, Xmat)

# model performance plots
perfs <- subtypePerf(calls, Y)

library(gridExtra)
x <- grid.arrange(perfs[[1]]$plot,perfs[[2]]$plot,perfs[[3]]$plot,perfs[[4]]$plot,perfs[[5]]$plot,perfs[[6]]$plot, ncol=6, nrow=1 )

save(ens, file='~/ensemble_model_signatures_genes.rda')
