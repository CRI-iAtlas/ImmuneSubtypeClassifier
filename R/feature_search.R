
library(readr)

# read in formatted data
ebpp <- readr::read_csv('/data/robencla_work/data/formatted_iatlas_485/EBpp_pancancer.csv')
ebpp$Label <- sapply(ebpp$Label, function(a) paste0('C',a))


# Expected that data will be
# The label will be one of the columns, and have column name: "label_name".
# "feature_list" will be some list of column names.
# "this_label" will be the base label to compare to others
# "this_thresh" will not report feature-pairs with scores below thresh

feature_search <- function(data,  # samples in rows, features in cols
                           label_col_name, # column name containing labels
                           sample_col_name,  # column name containing sample barcodes / sample names
                           this_label,  # the base label to compare to others
                           this_thresh # will not report feature-pairs with prop_diffs below thresh
                          ) {


  # recorded diff
  diff_list <- c()

  # the pair count/index
  k <- 1

  # list of results
  resl <- list()

  # the unique labels in the data
  labels <- as.vector(data[label_col_name])[[1]]
  labels <- sort(unique(labels))

  #  index into this cluster
  idx <- data[,label_col_name] == this_label
  jdx <- data[,label_col_name] != this_label

  cidx <- !(colnames(data) %in% c(label_col_name,sample_col_name))
  cjdx <- !(colnames(data) %in% c(label_col_name,sample_col_name))

  # data labeled as THIS cluster
  labdat <- data[idx, cidx]
  # data labels as NOT THIS cluster
  notdat <- data[jdx,cjdx]

  # number of samples in each category
  m <- sum(idx)
  n <- sum(jdx)

  # number of genes in the table
  ngenes <- ncol(labdat)

  # for each pair of genes
  for (i in 1:(ngenes-1)){

    for (j in (i+1):(ngenes)) {

      # proportion of samples showing pattern gene1 > gene2
      prop1 <- (sum(labdat[,i] > labdat[,j])/m) # proportion i>j in this cluster
      prop2 <- (sum(notdat[,i] > notdat[,j])/n) # proportion i>j in others

      # if the difference in proportions, between cluster and not-cluster
      # is greater than 75%, then record the gene pair.
      if (abs(prop1-prop2) > this_thresh) {

        # next line: distnance between the two genes, within sample grouping.
        dist1 <- as.vector(labdat[,i] - labdat[,j])[[1]]
        dist2 <- as.vector(notdat[,i] - notdat[,j])[[1]]

        t_dist1   <- t.test(dist1, na.rm = T)$statistic
        t_dist2   <- t.test(dist2, na.rm = T)$statistic
        med_dist1 <- median(dist1, na.rm = T)
        med_dist2 <- median(dist2, na.rm=T)
        avg_dist1 <- mean(dist1, na.rm = T)
        avg_dist2 <- mean(dist2, na.rm=T)
        n1        <- length(dist1)
        n2        <- length(dist2)

        score <- abs(prop1-prop2)*(-1 * mean(dist1) * mean(dist2))

        g1 <- colnames(labdat)[i]
        g2 <- colnames(labdat)[j]

        if (score > 0) {
          #capture this result
          resl[[k]] <- c(this_label,
                         k,i,j,
                         g1,g2,
                         prop1,prop2,
                         (prop1-prop2),
                         score,
                         t_dist1, t_dist2,
                         med_dist1, med_dist2,
                         avg_dist1, avg_dist2,
                         n1, n2) # want diff to be 1 or -1
          k <- k+1
        }
      } else {
        diff_list <- c(diff_list, abs(prop1-prop2))
      }
    }
  } # end gene pair loop

  if (length(resl) > 0) {
    df <- do.call('rbind', resl)
    colnames(df) <- c('Label','k','i','j','Feature1','Feature2',
                      'Prop1','Prop2','PropDiff','Score','T_dist1','T_dist2',
                      'MedDist1','MedDist2','AvgDist1','AvgDist2','N1','N2')

    return(df)
  } else {
    return(diff_list)
  }
}


df1 <- feature_search(data = ebpp,
                      label_col_name = 'Label',
                      sample_col_name = 'Barcode',
                      this_label = 'C1',
                      this_thresh = 0.33)

df2 <- feature_search(data = ebpp,
                      label_col_name = 'Label',
                      sample_col_name = 'Barcode',
                      this_label = 'C2',
                      this_thresh = 0.33)

df3 <- feature_search(data = ebpp,
                      label_col_name = 'Label',
                      sample_col_name = 'Barcode',
                      this_label = 'C3',
                      this_thresh = 0.33)

df4 <- feature_search(data = ebpp,
                      label_col_name = 'Label',
                      sample_col_name = 'Barcode',
                      this_label = 'C4',
                      this_thresh = 0.33)

df5 <- feature_search(data = ebpp,
                      label_col_name = 'Label',
                      sample_col_name = 'Barcode',
                      this_label = 'C5',
                      this_thresh = 0.33)

df6 <- feature_search(data = ebpp,
                      label_col_name = 'Label',
                      sample_col_name = 'Barcode',
                      this_label = 'C6',
                      this_thresh = 0.33)


df <- rbind(df1,df2,df3,df4,df5,df6)


write.csv(df, file = '/data/robencla_work/ebpp_487_features.csv')


df <- readr::read_csv('/data/robencla_work/ebpp_487_features.csv')



# pull out the feature pair lists
label <- 'C6'
df1 <- df[df$Label == label,]
df1 <- as.data.frame( df1[order(df1$Score, decreasing = T),] )
x <- lapply(1:10, function(a) c( as.character(df1[a, 'Feature1']),
                                 as.character(df1[a, 'Feature2']) ))
paste(unlist(x), collapse = "','")


