

library(readr)
library(dplyr)
library(stringr)

path <- '~/Work/iAtlas/Subtypes/Subtype-Classifier/GDC_Test'

setwd(path)

dirs <- list.files('dat_dirs/')

datlist <- list()

for (di in dirs) {

  print(di)
  files <- list.files(paste0(path,'/dat_dirs/',di))

  if (str_detect(files[1], 'htseq.counts')) {
    dat <- read_tsv(paste0(path,'/dat_dirs/',di,'/',files[1]), col_names = F)
  }
  else if (str_detect(files[2], 'htseq.counts')) {
    dat <- read_tsv(paste0(path,'/dat_dirs/',di,'/',files[2]), col_names = F)
  }

  datlist[[di]] <- dat

  }

datmat <- datlist[[1]]
colnames(datmat) <- c('GeneID', names(datlist)[1])

for (i in 2:length(datlist)) {
  print(i)
  dn <- names(datlist)[i]
  datn <- datlist[[dn]]
  colnames(datn) <- c('GeneID', dn)
  datmat <- inner_join(datmat, datn)

}

dim(datmat)
#[1] 60488   612

library("AnnotationDbi")
library("org.Hs.eg.db")

ensemble <- str_split(datmat$GeneID, pattern = '\\.')
head(ensemble)
ensemble <- unlist(lapply(ensemble, function(a) a[1]))

symbols <- mapIds(org.Hs.eg.db, keys=ensemble, column="SYMBOL", keytype="ENSEMBL", multiVals="first")

length(symbols)

idx <- !duplicated(symbols)

subsyms <- symbols[idx]
subsyms[is.na(subsyms)] <- "ISNA123"

datmat2 <- as.data.frame(datmat[idx,-1])
rownames(datmat2) <- subsyms

# Load the package required to read JSON files.
library("rjson")

# Give the input file name to the function.
result <- fromJSON(file = "metadata.cart.2019-06-07.json")

length(result)
#[1] 611

tcgaBarcodes <- unlist(lapply(result, function(a) a$associated_entities[[1]]$entity_submitter_id))

write.table(datmat2, file='gdc_kirc_dat.csv', sep=',', quote = F)


write.table(data.frame(AliquotBarcode=tcgaBarcodes, ID=colnames(datmat2)), file='gdc_kirc_barcodes.csv', sep=',', quote=F, row.names = F)


