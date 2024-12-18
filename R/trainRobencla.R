
### training robencla ###

library(robencla)

trainRobencla <- function(data,
                          label_name='ClusterLabel',
                          sample_id = 'SampleBarcode',
                          pair_list
                          ) {

  obj_name <- Sys.time()

  # Our classifier object named Anne.  Hi Anne!
  mod <- Robencla$new(obj_name)

  # XGBoost parameters to pass to each sub-classifier in the ensembles
  params <- list(
    max_depth=12,    # "height" of the tree, 6 is actually default. I think about 12 seems better.  (xgboost parameter)
    eta=0.3,        # this is the learning rate. smaller values slow it down, more conservative   (xgboost parameter)
    nrounds=64,     # number of rounds of training, lower numbers less overfitting (potentially)  (xgboost parameter)
    early_stopping_rounds=2, # number of rounds without improvment stops the training (xgboost early_stopping_rounds)
    nthreads=4,     # parallel threads
    gamma=0.2,        # min loss required to again split a leaf node. higher ~ more conservative (xgboost parameter)
    lambda=1.2,     # L2 regularization term on weights, higher number ~ more conservative (xgboost parameter)
    alpha=0.2,      # L1 regularization term on weights. higher number ~ more conservative (xgboost parameter)
    size=11,        # Size of the ensemble, per binary prediction
    sample_prop=0.8, # The percentage of data used to train each ensemble member.
    feature_prop=0.8, # The percentage of data used to train each ensemble member.
    subsample=0.8,    # the xgboost machines subsample at this rate.
    combine_function='median',  # How the ensemble should be combined. Only median currently.
    verbose=0)
  ###More on the xgboost parameters: https://xgboost.readthedocs.io/en/latest/parameter.html

  # First we use the training data
  mod$train (data_frame=data,
             label_name=label_name,
             sample_id=sample_id,
             data_mode=c('namedpairs'), # allpairs,namedpairs,sigpairs,quartiles,tertiles,binarize,ranks,original #
             signatures=NULL,
             pair_list=pair_list,  # subset to these genes.
             params=params)

  return(mod)
}

# define the pair list; top 10 scoring pairs in each class
pair_list <- list(
  C1=c("B2M","COL3A1","B2M","COL1A2","COL1A2","HLA-B","COL3A1","HLA-B","APOE","COL6A3","APOE","SDC1","APOE","COL6A1","APOE","MMP14","APOE","MMP2","COL1A2","SPARC"),
  C2=c('IFI27','RHOB','IFI6','RHOB','IFI27','LRP1','IFI27','NPC2','IFI6','NPC2','IFI27','TAGLN','RHOB','STAT1','IFI6','LRP1','IFI27','IGFBP3','IFI27','IGFBP4'),
  C3=c('COL3A1','SPARC','IGFBP4','JUP','CD59','JUP','NPC2','SLC25A5','IFI27','NPC2','IGFBP4','PFN1','HNRNPA2B1','IGFBP4','CCT5','NPC2','IFI27','MET','IGFBP4','TPI1'),
  C4=c('APOE','COL3A1','APOE','COL1A2','COL3A1','SPARC','APOC1','DSP','COL1A2','SPARC','APOC1','JUP','APOC1','COL6A3','APOC1','LYZ','APOC1','SDC1','DSP','RHOB'),
  C5=c('APOE','FN1','FN1','SPARC','APOE','COL3A1','APOE','COL1A2','B2M','SPARC','APOE','B2M','FN1','SLC1A3','COL3A1','SLC1A3','APOE','HLA-B','COL1A2','SLC1A3'),
  C6=c('B2M','COL3A1','B2M','COL1A2','COL6A3','ENO1','BSG','COL6A3','COL6A3','HNRNPA2B1','COL6A3','MYL6','BSG','COL6A1','COL6A3','DSP','BSG','MMP2','JUP','MMP2')
  )

# read in formatted data
ebpp <- readr::read_csv('/data/robencla_work/data/formatted_iatlas_485/EBpp_pancancer.csv')
#ebpp <- as.data.table(ebpp)
ebpp$Label <- sapply(ebpp$Label, function(a) paste0('C',a))


emod <- trainRobencla(ebpp,
                      label_name = 'Label',
                      sample_id  = 'Barcode',
                      pair_list)


emod$predict(data_frame = ebpp,
             label_name = 'Label',
             sample_id  = 'Barcode')


res0 <- emod$results(include_label = T)

table(res0$BestCalls, emod$test_label)

emod$classification_metrics()

#     C1   C2   C3   C4   C5   C6
#C1 2362   27   45   16    0    8
#C2   16 2557    3   13    0    1
#C3   29    4 2320   20    4   12
#C4    4    2   26 1084   20    1
#C5    0    0    1   26  361    0
#C6    5    2    2    0    0  158
