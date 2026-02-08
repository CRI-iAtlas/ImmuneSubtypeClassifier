
feature_search <- function() {


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

          # next line: distance between the two genes, within sample grouping.
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

}

#######################################################33



compute_t_score <- function(genepairs) {
  # Combines: proportion difference × statistical strength of separation in both groups
  TScore <- abs(genepairs$Prop1 - genepairs$Prop2) *
    abs(genepairs$T_dist1) *
    abs(genepairs$T_dist2)

  return(TScore)
}

genepairs$T_score <- compute_t_score(genepairs)


select_top_pairs <- function(genepairs, N = 50, metric = "Score", filtergenes) {

  top_pairs <- genepairs %>%
    group_by(Label) %>%
    arrange(desc(.data[[metric]])) %>%
    mutate(rank = row_number()) %>%
    ungroup() %>%
    arrange(Label, rank)

  # Greedily select N pairs per label, avoiding feature reuse
  selected <- list()

  for (label in unique(top_pairs$Label)) {
    label_pairs <- top_pairs %>% filter(Label == label)
    used_features <- character(0)
    label_selected <- tibble()

    for (i in 1:nrow(label_pairs)) {
      pair <- label_pairs[i, ]

      # Check if either feature already used
      if (!pair$Feature1 %in% used_features &&
          !pair$Feature2 %in% used_features) {
        label_selected <- bind_rows(label_selected, pair)
        used_features <- c(used_features, pair$Feature1, pair$Feature2)

        if (nrow(label_selected) >= N) break
      }
    }

    selected[[label]] <- label_selected
  }

  bind_rows(selected)
}


genepairs <- readr::read_csv('../data/PanCancer_TCGA/ebpp_485_features.csv.gz')
genepairs2 <- genepairs %>% filter( (Feature1 %in% colnames(test_data))  & (Feature2 %in% colnames(test_data)) )

# Use it
top_pairs <- select_top_pairs(genepairs2, N = 50, metric = "T_score")

# Verify
top_pairs %>%
  group_by(Label) %>%
  summarise(n_pairs = n(),
            n_unique_features = n_distinct(c(Feature1, Feature2)))


print_pairs_list <- function(top_pairs, N = 10) {
  top_pairs_subset <- top_pairs %>%
    group_by(Label) %>%
    slice_head(n = N) %>%
    ungroup()

  cat("list(\n")

  labels <- unique(top_pairs_subset$Label)

  for (i in seq_along(labels)) {
    label <- labels[i]
    label_pairs <- top_pairs_subset %>%
      filter(Label == label) %>%
      arrange(desc(T_score))  # or whatever ordering

    # Interleave Feature1 and Feature2
    genes <- c(rbind(label_pairs$Feature1, label_pairs$Feature2))

    # Format as quoted strings
    gene_str <- paste0('"', genes, '"', collapse = ", ")

    cat("  ", label, " = c(", gene_str, ")", sep = "")

    # Add comma unless it's the last one
    if (i < length(labels)) {
      cat(",\n")
    } else {
      cat("\n")
    }
  }

  cat(")\n")
}

# Usage
print_pairs_list(top_pairs, N = 40)

# Usage
#gene_list <- pairs_to_list(top_pairs, N = 60)

gene_list_60 <- list(
  C1 = c("HERC6", "ZWILCH", "IFIT1", "NUP85", "LMNB2", "MX1", "IFIT3", "SMC2", "CDCA4", "IFI44", "IFI44L", "PNO1", "IFIT2", "RFC3", "RASSF4", "RBM14", "HERC5", "WDHD1", "IFIH1", "NUDT15", "ACTL6A", "OAS2", "DDX60", "MYCBP", "RTP4", "SKA1", "GIMAP4", "LYAR", "CCT5", "STAT1", "CHEK1", "RSAD2", "IFI6", "SNRPB", "HLA-E", "TPM3", "DDX58", "SNRPA1", "BRCA2", "MNDA", "CD86", "PLK4", "BAK1", "IFI35", "LY86", "POLE2", "CPEB4", "PSMD12", "EPHB3", "HLA-DMA", "ISG15", "MYBL2", "CSF1R", "PLAUR", "BRIP1", "BTK", "CENPO", "LST1", "HLA-DRB1", "MMP11", "MMP14", "MSN", "PA2G4", "PSME1", "FNBP1", "IPO4", "APOE", "COL6A3", "CENPW", "OASL", "EPHB4", "PSMB8", "GBP1", "LOXL2", "B2M", "COL1A2", "IL10RA", "TEAD4", "GGH", "OAS1", "MMP3", "MT3", "COL3A1", "HLA-B", "CDK2", "SAMHD1", "HADHB", "MCM7", "EVI2A", "PMAIP1", "MMP1", "SLC1A3", "APOC1", "EIF4EBP1", "CD84", "HAS2"),
  C2 = c("OAS2", "WSB2", "OAS3", "SF3A1", "HADHB", "MX1", "CPEB4", "IFI44", "CTNNA1", "STAT1", "IFIH1", "RNF41", "IFI27", "RHOB", "ID2", "IFIT3", "TAP1", "VDAC1", "MPP1", "OASL", "NMI", "TCF7L2", "HERC6", "RMND5B", "EPS15", "OAS1", "CXCL10", "SH2B3", "GNG11", "RSAD2", "DDX60", "NUPL1", "IFI44L", "RASSF4", "IFI6", "RHOC", "GBP1", "NEO1", "ISG15", "SNX17", "IRF1", "UBE2J1", "HERC5", "MERTK", "ALKBH7", "IFI35", "DDX58", "SH3BP5L", "IFITM1", "VAT1", "ARHGAP1", "IFI16", "BST2", "EIF2AK1", "IFIT1", "RHOQ", "ID3", "SAMD9", "USPL1", "ZWILCH", "CD3E", "FLI1", "DLEU1", "POLE2", "BBC3", "RTP4", "GSTCD", "PLK4", "IL15RA", "TBXAS1", "C1orf54", "CD3D", "COPS6", "PSMB8", "CDCA4", "SMO", "FAM167A", "SKA1", "CCL5", "RASSF2", "MYBL2", "SLC1A3", "CENPW", "MAP3K10", "IFI30", "IGF2R", "MKKS", "SNRPD1", "CD8A", "ITGB3", "BAK1", "SAR1B", "BRCA2", "SLC25A40", "H2AFZ", "SRP9", "LY86", "SP140", "ACTL6A", "HADH"),
  C3 = c("MYBL2", "RMND5B", "ARHGAP15", "SKA1", "EPS15", "LMNB2", "ARHGAP1", "H2AFZ", "PLK4", "SLC25A40", "CDCA4", "FLI1", "CCT5", "IGF2R", "ACTL6A", "RNF41", "C1orf54", "CENPW", "MCM7", "SNX17", "APITD1", "ZWILCH", "MNAT1", "WDHD1", "GIMAP4", "SMC2", "OAS3", "UBE2J1", "ID2", "MCM3", "SH2B3", "SNRPD1", "RFC3", "USPL1", "HADHB", "PSMA7", "BAK1", "MPP1", "BRCA2", "PIK3CG", "NME1", "RHOQ", "CHEK1", "VEGFC", "EIF4EBP1", "GNG11", "CPVL", "RUVBL1", "CD33", "POLE2", "CPEB4", "PSMD14", "PXN", "SNRPB", "FAS", "NMI", "ITGB3", "PHF19", "IFI27", "MET", "FNBP1", "OAS2", "BBC3", "NUDT1", "PAICS", "SF3A1", "CEBPD", "ISG15", "CCBL2", "SNRPA1", "ALKBH7", "MRTO4", "LYAR", "MARVELD2", "CENPO", "DLEU1", "COL8A1", "PLAUR", "IGFBP4", "PFN1", "ADAMTS1", "SLC16A1", "NPC2", "SLC25A5", "OAS1", "WIPF1", "BTG1", "IFITM1", "NOTCH2", "PARP1", "FGR", "OASL", "CD59", "JUP", "BRIP1", "CD84", "CTNNA1", "TPM3", "PSMD2", "VAT1"),
  C4 = c("CD3E", "LY86", "CCL5", "MPP1", "DLEU1", "LCK", "CD33", "CD3D", "C1orf54", "CD8A", "MRPS28", "NUDT15", "HSPG2", "VAT1", "PTPRC", "VSIG4", "IRF1", "MKKS", "C3AR1", "LCP2", "CD48", "FCGR1A", "HADHB", "ITGA6", "CPEB4", "EPHA2", "ARRB2", "SYK", "CD79A", "RPP40", "SAR1B", "TNFAIP3", "APOO", "CDCA4", "ARHGDIA", "JUP", "APOE", "COL3A1", "CD14", "LCP1", "ID2", "MMP11", "COL8A1", "MERTK", "APOC1", "DSP", "IKZF1", "RNASE6", "IL7R", "TNFSF13B", "CD3G", "TAS2R5", "LST1", "LTB", "AIMP2", "MYCBP", "ALKBH7", "SEMA3F", "CD52", "COQ2", "COL1A2", "SPARC", "COPS6", "HN1L", "MAP3K10", "RUNX3", "COL6A3", "RHOB", "FCER1G", "PLAUR", "ATP6V0B", "THBS1", "MMP2", "SDCBP", "GPLD1", "MS4A1", "IFIT3", "OAS2", "ITGB8", "SLC1A3", "MMP3", "MT3"),
  C5 = c("COL6A3", "SAR1A", "CBX1", "HN1L", "COL3A1", "PSMD2", "MAPRE1", "PDIA4", "PFN1", "SLC1A3", "IGFBP4", "SF3A1", "RASSF2", "RPN1", "EMP2", "SKA2", "EPHB4", "RBM14", "CASP8", "CEP78", "MRPS28", "MYCBP", "SDC1", "SNX17", "CPEB4", "HSPG2", "MMP11", "PNO1", "FAM167A", "IRF1", "PSME1", "RASSF4", "LRRC40", "NMI", "IFIT2", "OAS2", "ARRB2", "PSMB8", "CDCA4", "MNAT1", "NEO1", "THBS1", "TNFRSF1A", "WSB2", "FGD1", "MYBL2", "IFRD2", "RMND5B", "HADHB", "JUP", "DCK", "OAS1", "MARVELD2", "SLC25A40", "COPS6", "SNRPB", "PLAUR", "USPL1", "IPO4", "PHF19", "SRP9", "TPM3", "ID2", "TIMP1", "COL1A2", "LRP1", "FARSA", "MCM3", "IFIT3", "OAS3", "EPHA2", "PITPNC1", "CCL5", "GSTCD", "ITGA5", "RNF41", "APOO", "BAK1", "CORO1C", "STAT1", "MSN", "PHLDA1", "FNBP1", "PML", "MARCKSL1", "MYL6", "CD33", "CD3E", "CSF1R", "IFI30", "MERTK", "TNFRSF12A", "LCK", "SP140", "EVI2A", "LYN", "MAP3K10", "PDLIM7", "DVL1", "IFI27"),
  C6 = c("LAIR1", "STRA13", "ITGB2", "PA2G4", "FCGR2A", "MRPL12", "FARSA", "LHFPL2", "MNDA", "MRPS28", "SELPLG", "SNRPD1", "LAPTM5", "XRCC6", "CD86", "LYAR", "NCKAP1L", "WDR77", "CD53", "RUVBL1", "CD33", "RPP40", "AIMP2", "C3AR1", "DOCK2", "UMPS", "SLC25A5", "THBS2", "NCF2", "NUDT15", "CYBB", "LSM4", "MRPS16", "WIPF1", "FCER1G", "MRTO4", "FCGR2B", "RFC3", "LCP2", "SNRPA1", "EVI2B", "RNF138", "HCLS1", "TOMM40", "PPIH", "RNASE6", "IMP4", "MYO1F", "BTK", "CHEK1", "FCGR3A", "PAICS", "RHOG", "SNRPA", "FMNL1", "NUP85", "ALOX5AP", "PNO1", "APOO", "SAMSN1", "C19orf48", "SH2B3", "IL10RA", "LSM3", "CSF1R", "NME1", "CYCS", "LAMA4", "SNRPB", "TNFRSF1A", "APITD1", "LY86", "CENPO", "PIK3CG", "FGR", "NOP16", "CDCA4", "HCK", "COL6A3", "HNRNPA2B1", "CD163", "COPS6", "CCT5", "ITGB5", "CENPN", "TBXAS1", "PLEK", "SKA2", "PLAUR", "SNRPC", "MRPL37", "TNFRSF1B", "BCCIP", "COL8A1", "LST1", "NUP35", "MMP2", "TPI1", "MKKS", "SLC7A7")
)


gene_list_20 <- list(
  C1 = c("HERC6", "ZWILCH", "IFIT1", "NUP85", "LMNB2", "MX1", "IFIT3", "SMC2", "CDCA4", "IFI44", "IFI44L", "PNO1", "IFIT2", "RFC3", "RASSF4", "RBM14", "HERC5", "WDHD1", "IFIH1", "NUDT15", "ACTL6A", "OAS2", "DDX60", "MYCBP", "RTP4", "SKA1", "GIMAP4", "LYAR", "CCT5", "STAT1", "CHEK1", "RSAD2", "IFI6", "SNRPB", "HLA-E", "TPM3", "DDX58", "SNRPA1", "BRCA2", "MNDA"),
  C2 = c("OAS2", "WSB2", "OAS3", "SF3A1", "HADHB", "MX1", "CPEB4", "IFI44", "CTNNA1", "STAT1", "IFIH1", "RNF41", "IFI27", "RHOB", "ID2", "IFIT3", "TAP1", "VDAC1", "MPP1", "OASL", "NMI", "TCF7L2", "HERC6", "RMND5B", "EPS15", "OAS1", "CXCL10", "SH2B3", "GNG11", "RSAD2", "DDX60", "NUPL1", "IFI44L", "RASSF4", "IFI6", "RHOC", "GBP1", "NEO1", "ISG15", "SNX17"),
  C3 = c("MYBL2", "RMND5B", "ARHGAP15", "SKA1", "EPS15", "LMNB2", "ARHGAP1", "H2AFZ", "PLK4", "SLC25A40", "CDCA4", "FLI1", "CCT5", "IGF2R", "ACTL6A", "RNF41", "C1orf54", "CENPW", "MCM7", "SNX17", "APITD1", "ZWILCH", "MNAT1", "WDHD1", "GIMAP4", "SMC2", "OAS3", "UBE2J1", "ID2", "MCM3", "SH2B3", "SNRPD1", "RFC3", "USPL1", "HADHB", "PSMA7", "BAK1", "MPP1", "BRCA2", "PIK3CG"),
  C4 = c("CD3E", "LY86", "CCL5", "MPP1", "DLEU1", "LCK", "CD33", "CD3D", "C1orf54", "CD8A", "MRPS28", "NUDT15", "HSPG2", "VAT1", "PTPRC", "VSIG4", "IRF1", "MKKS", "C3AR1", "LCP2", "CD48", "FCGR1A", "HADHB", "ITGA6", "CPEB4", "EPHA2", "ARRB2", "SYK", "CD79A", "RPP40", "SAR1B", "TNFAIP3", "APOO", "CDCA4", "ARHGDIA", "JUP", "APOE", "COL3A1", "CD14", "LCP1"),
  C5 = c("COL6A3", "SAR1A", "CBX1", "HN1L", "COL3A1", "PSMD2", "MAPRE1", "PDIA4", "PFN1", "SLC1A3", "IGFBP4", "SF3A1", "RASSF2", "RPN1", "EMP2", "SKA2", "EPHB4", "RBM14", "CASP8", "CEP78", "MRPS28", "MYCBP", "SDC1", "SNX17", "CPEB4", "HSPG2", "MMP11", "PNO1", "FAM167A", "IRF1", "PSME1", "RASSF4", "LRRC40", "NMI", "IFIT2", "OAS2", "ARRB2", "PSMB8", "CDCA4", "MNAT1"),
  C6 = c("LAIR1", "STRA13", "ITGB2", "PA2G4", "FCGR2A", "MRPL12", "FARSA", "LHFPL2", "MNDA", "MRPS28", "SELPLG", "SNRPD1", "LAPTM5", "XRCC6", "CD86", "LYAR", "NCKAP1L", "WDR77", "CD53", "RUVBL1", "CD33", "RPP40", "AIMP2", "C3AR1", "DOCK2", "UMPS", "SLC25A5", "THBS2", "NCF2", "NUDT15", "CYBB", "LSM4", "MRPS16", "WIPF1", "FCER1G", "MRTO4", "FCGR2B", "RFC3", "LCP2", "SNRPA1")
)

load('inst/extdata/geneSetSymbols.rda')

filter_sigs <- function(sigs, data) {
  # Get column names from data (exclude non-gene columns)
  available_genes <- colnames(data)

  # Filter each signature to only include available genes
  filtered_sigs <- lapply(sigs, function(sig_genes) {
    # Keep only genes that exist in the data
    filtered <- sig_genes[sig_genes %in% available_genes]

    # Warn if signature lost too many genes
    n_original <- length(sig_genes)
    n_kept <- length(filtered)
    if (n_kept == 0) {
      warning(sprintf("Signature has no genes in data! Original: %d genes", n_original))
    } else if (n_kept < n_original * 0.5) {
      warning(sprintf("Signature lost >50%% of genes: kept %d/%d", n_kept, n_original))
    }

    return(filtered)
  })

  # Remove empty signatures
  filtered_sigs <- filtered_sigs[sapply(filtered_sigs, length) > 0]

  return(filtered_sigs)
}

ebpp <- readr::read_csv('../data/PanCancer_TCGA/EBpp_pancancer_476.csv.gz')
xena <- readr::read_csv('../test_run_imm/data/xena_tpm/xena_rsem_tpm.csv')
genesetsymbols <- filter_sigs(genesetsymbols, ebpp)


# Conservative parameters for better C4/C6
conservative_params <- list(
  max_depth = 6,
  eta = 0.3,
  nrounds = 32,
  early_stopping_rounds = 2,
  gamma = 0.5,
  lambda = 2.0,
  alpha = 0.5,
  ensemble_size = 7,
  sample_prop = 0.7,
  feature_prop = 0.7,
  subsample = 0.7
)

library(ImmuneSubtypeClassifier)

result <- build_robencla_classifier(
  data_path='../data/PanCancer_TCGA/EBpp_pancancer_476.csv.gz',
  test_path='../test_run_imm/data/xena_tpm/xena_rsem_tpm.csv', #NULL
  output_path = '../models/immune_conservative.rds',
  pair_list = gene_list_20,
  sig_list = genesetsymbols,
  param_list = conservative_params,
  data_mode = c("sigpairs"), # , "tertiles" namedpairs
  train_fraction = 0.8,
  seed = 412
)



