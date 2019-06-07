

# Examining genes in the model #

load('~/Code/iatlas/ImmuneSubtypeClassifier/data/ensemble_model.rda')

res0 <- list()
for (j in 1:6) {  #-- length 6

  res1 <- list()
  for (i in 1:5) {   #-- length 5

    print(paste0(j, '  ', i))
    ei <- ens[[i]]  # get the ensemble member out, list of 6
    m <- ei[[j]]$bst    # get the subtype classifier out
    g <- xgboost::xgb.importance(model=m)
    print(head(g$Feature))
    res1[[i]]  <- g
  }

  res0[[j]] <- res1
}



jacIdx <- function(a,b, topn = 50) {
  length(intersect( a$Feature[1:topn], b$Feature[1:topn] ) ) / length( union ( a$Feature[1:topn], b$Feature[1:topn] ) )
}



mat <- matrix(data=0, ncol=30, nrow=30)
annot <- data.frame()

for (si in 1:6) {
  for (ti in 1:6) {

    for (i in 1:5) {
      for (j in 1:5) {

        x <- res0[[si]][[i]]
        y <- res0[[ti]][[j]]

        mi <- (si * 5) - (5 - i)
        mj <- (ti * 5) - (5 - j)

        val <- (jacIdx(x,y, 100))

        print(paste0(mi, '  ', mj))

        annot <- rbind(annot, data.frame(Subtype1=si, Subtype2=ti, Ens1=i, Ens2=j, JaccardIdx=val))

        mat[mi,mj]  <- val
      }
    }
  }
}

dfcol <- unique(annot[,c(1,3)])
rownames(dfcol) <- 1:30
rownames(mat)   <- 1:30
colnames(mat)   <- 1:30

library(pheatmap)
pheatmap(mat, scale = 'none', cluster_rows = F, cluster_cols = F, annotation_col = dfcol)


allgenes <- data.frame()

for (si in 1:6) {
    for (i in 1:5) {
        x <- res0[[si]][[i]]
        for (xi in 1:nrow(x)) {
          allgenes <- rbind(allgenes, data.frame(Subtype1=si, EnsembleMember=i, GeneNum=xi, Gene=x$Feature[xi], Gain=x$Gain[xi]))
        }
    }
}


length(unique(allgenes$Gene))
#[1] 2813

length(unique(allgenes$Gene[allgenes$Gain > 0.01])) #***
#[1] 200

length(unique(allgenes$Gene[allgenes$Gain > 0.05]))
#[1] 30

g <- ggplot(data=allgenes, aes(x=GeneNum, y=Gain, col=as.factor(Subtype1)))
g + geom_point() + xlim(c(0,100)) + geom_hline(yintercept = 0)





