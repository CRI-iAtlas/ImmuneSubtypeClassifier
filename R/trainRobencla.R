




#' Edit a pair_list, remove pairs containing identified items.
#'
#' Return the cleaned pair_list.
#'
#' @param pair_list Vector or named list of vectors. Pairs ordered as (1,2), (3,4), (5,6), etc.
#' @param items_to_remove Vector of items; pairs containing any of these will be removed.
#'
#' @return new_pair_list Same structure as input (vector or list) with pairs removed.
#'
#' @keywords internal
editPairList <- function(pair_list, items_to_remove) {

  # Handle vector input
  if (is.vector(pair_list) && !is.list(pair_list)) {
    # Check even length
    if (length(pair_list) %% 2 != 0) {
      stop("pair_list vector must have even length")
    }

    # Create pair indices
    pair_idx <- seq(1, length(pair_list), by = 2)

    # Check which pairs contain items to remove
    keep <- vapply(pair_idx, function(i) {
      !any(pair_list[i] %in% items_to_remove) &&
        !any(pair_list[i + 1] %in% items_to_remove)
    }, logical(1))

    # Extract kept pairs
    kept_indices <- sort(c(pair_idx[keep], pair_idx[keep] + 1))
    return(pair_list[kept_indices])
  }

  # Handle list input
  if (is.list(pair_list)) {
    new_list <- lapply(pair_list, function(vec) {
      if (length(vec) %% 2 != 0) {
        warning("Skipping list element with odd length")
        return(vec)
      }

      pair_idx <- seq(1, length(vec), by = 2)

      keep <- vapply(pair_idx, function(i) {
        !any(vec[i] %in% items_to_remove) &&
          !any(vec[i + 1] %in% items_to_remove)
      }, logical(1))

      kept_indices <- sort(c(pair_idx[keep], pair_idx[keep] + 1))
      vec[kept_indices]
    })

    return(new_list)
  }

  stop("pair_list must be a vector or list of vectors")
}


#' Default Feature Pair List for Immune Subtype Classification
#'
#' Returns the default named pair list used for immune subtype classification.
#' Each class (C1-C6) has 10 feature pairs (20 genes) selected based on
#' discriminative power.
#'
#' @return A named list with elements C1 through C6, each containing a
#'   character vector of 20 gene names representing 10 feature pairs.
#'
#' @examples
#' pair_list <- getFeaturesPairList()
#' names(pair_list)
#'
#' @export
getFeaturesPairList <- function() {
  data(gene_list, envir = environment())
}


#' Default Gene List for Immune Subtype Classification
#'
#' Returns the default named genes used for immune subtype classification.
#'
#' @return A named list with the 36 used genes.
#'
#' @examples
#' gene_list <- getFeaturesGeneTable()
#'
#' @export
getFeaturesGeneTable <- function() {
  gs <- unique(unlist(getFeaturesPairList()))
  data(ebpp_gene, envir = environment())
  df <- ebpp_genes_full[ebpp_genes_full$Symbol %in% gs,]
  # 2. Convert factors to characters (Modifying 'df')
  cols_to_fix <- c("Symbol", "Entrez", "Ensembl")
  df[cols_to_fix] <- lapply(df[cols_to_fix], as.character)
  return(
    df
  )
}


#' Train a Robencla Classifier for Immune Subtype Classification
#'
#' Trains a robencla ensemble classifier using named feature pairs and
#' XGBoost parameters tuned for immune subtype prediction.
#'
#' @param data A data.frame containing the training data with features,
#'   labels, and sample identifiers.
#' @param label_name Character string specifying the column name containing
#'   class labels. Default is "ClusterLabel".
#' @param sample_id Character string specifying the column name containing
#'   sample identifiers. Default is "SampleBarcode".
#' @param pair_list A named list where each element corresponds to a class
#'   and contains a character vector of feature names. Features are paired
#'   sequentially (1-2, 3-4, 5-6, ...) for named pair transformation.
#' @param sig_list List of signatures, named list of genes that belong to a
#'    signature.  The signatures will be used to compute P(sig_a__g_i > sig_b__g_j)
#' @param param_list A named list of training parameters. If NULL, uses defaults.
#'   Supported parameters:
#'   \itemize{
#'     \item max_depth: Integer, maximum tree depth (default: 12)
#'     \item eta: Numeric, learning rate (default: 0.3)
#'     \item nrounds: Integer, maximum boosting rounds (default: 64)
#'     \item early_stopping_rounds: Integer, early stopping patience (default: 2)
#'     \item nthreads: Integer, parallel threads (default: 4)
#'     \item gamma: Numeric, minimum loss reduction (default: 0.2)
#'     \item lambda: Numeric, L2 regularization (default: 1.2)
#'     \item alpha: Numeric, L1 regularization (default: 0.2)
#'     \item ensemble_size: Integer, number of ensemble models (default: 11)
#'     \item sample_prop: Numeric, sample proportion per model (default: 0.8)
#'     \item feature_prop: Numeric, feature proportion per model (default: 0.8)
#'     \item subsample: Numeric, XGBoost subsample ratio (default: 0.8)
#'     \item combine_function: Character, ensemble combination method (default: "median")
#'   }
#' @param verbose Integer, verbosity level. 0 for silent. Default is 0.
#' @param trim_model Logical, whether to call trim() on the model after
#'   training to reduce object size. Default is TRUE.
#'
#' @return A trained Robencla model object.
#'
#' @details
#' The classifier uses named feature pairs (namedpairs mode) where features
#' are paired sequentially within each class's pair_list. XGBoost parameters
#' are tuned for immune subtype classification tasks.
#'
#' For more information on XGBoost parameters, see:
#' \url{https://xgboost.readthedocs.io/en/latest/parameter.html}
#'
#' @examples
#' \dontrun{
#' # Define feature pairs for each class
#' pair_list <- list(
#'   C1 = c("GeneA", "GeneB", "GeneC", "GeneD"),
#'   C2 = c("GeneE", "GeneF", "GeneG", "GeneH")
#' )
#'
#' # Custom parameters
#' params <- list(
#'   max_depth = 6,
#'   eta = 0.1,
#'   nrounds = 200,
#'   ensemble_size = 21
#' )
#'
#' # Train the classifier
#' model <- trainRobencla(
#'   data = training_data,
#'   label_name = "Label",
#'   sample_id = "Barcode",
#'   pair_list = pair_list,
#'   param_list = params
#' )
#'
#' # Make predictions
#' model$predict(data_frame = test_data,
#'               label_name = "Label",
#'               sample_id = "SampleID")
#' results <- model$results()
#' }
#'
#' @seealso \code{\link[robencla]{Robencla}}
#'
#' @export
trainRobencla <- function(data,
                          label_name = "ClusterLabel",
                          sample_id = "SampleBarcode",
                          pair_list,
                          sig_list = NULL,
                          param_list = NULL,
                          data_mode=c("namedpairs"),
                          verbose = 0,
                          trim_model = TRUE) {

  # Default parameters
  default_params <- list(
    max_depth = 12,
    eta = 0.3,
    nrounds = 64,
    early_stopping_rounds = 2,
    nthreads = 4,
    gamma = 0.2,
    lambda = 1.2,
    alpha = 0.2,
    ensemble_size = 11,
    sample_prop = 0.8,
    feature_prop = 0.8,
    subsample = 0.8,
    combine_function = "median"
  )

  # Merge user params with defaults (user params override)
  if (!is.null(param_list)) {
    params <- modifyList(default_params, param_list)
  } else {
    params <- default_params
  }

  # Add verbose to params
  params$verbose <- verbose

  # Rename ensemble_size to size for robencla
  params$size <- params$ensemble_size
  params$ensemble_size <- NULL

  obj_name <- as.character(Sys.time())
  mod <- robencla::Robencla$new(obj_name)
  print(paste0("Robencla version: ", mod$version()))

  mod$train(
    data_frame = data,
    label_name = label_name,
    sample_id = sample_id,
    data_mode = data_mode,
    signatures = sig_list,
    pair_list = pair_list,
    params = params
  )

  if (trim_model) {
    mod$trim()
  }

  return(mod)
}


#' Build and Save Robencla Classifier
#'
#' Convenience function to build a robencla classifier from a data file,
#' evaluate performance, and save the trained model.
#'
#' @import data.table
#' @param data_path Character string, path to the training data CSV file.
#' @param test_path Character string, path to the testing data CSV file, can be NULL
#' @param output_path Character string, path where the trained model will
#'   be saved as an .rds file.
#' @param label_name Character string, column name containing class labels.
#'   Default is "Label".
#' @param sample_id Character string, column name containing sample IDs.
#'   Default is "Barcode".
#' @param pair_list Named list of feature pairs. If NULL, uses
#'   \code{getFeaturesPairList()}. Default is NULL.
#' @param sig_list List of signatures, named list of genes that belong to a
#'    signature.  The signatures will be used to compute P(sig_a__g_i > sig_b__g_j)
#' @param param_list Named list of training parameters. See \code{trainRobencla}
#'   for details. Default is NULL (uses default parameters).
#' @param data_mode List of data modalities to use. see gibbsdavidl/robencla
#' @param label_prefix Character string to prepend to numeric labels.
#'   Set to NULL to skip label transformation. Default is "C".
#' @param train_fraction Numeric, proportion of data to use for training
#'   when splitting. Use 1.0 for no split. Default is 1.0.
#' @param seed Integer, random seed for reproducible train/test splits.
#'   Default is 42.
#' @param evaluate Logical, whether to run evaluation and print metrics.
#'   Default is TRUE.
#'
#' @return The trained and trimmed Robencla model object (invisibly), along
#'   with test_data if a split was performed.
#'
#' @examples
#' \dontrun{
#' # With custom parameters
#' params <- list(
#'   max_depth = 6,
#'   eta = 0.1,
#'   nrounds = 200,
#'   ensemble_size = 21
#' )
#'
#' model <- build_robencla_classifier(
#'   data_path = "data/EBpp_pancancer.csv.gz",
#'   output_path = "models/robencla_trained_model.rds",
#'   label_name = "Label",
#'   sample_id = "Barcode",
#'   param_list = params,
#'   train_fraction = 0.8
#' )
#' }
#'
#' @export
build_robencla_classifier <- function(data_path,
                                      test_path = NULL,
                                      output_path = 'model.rds',
                                      label_name = "Label",
                                      sample_id = "Barcode",
                                      pair_list = NULL,
                                      sig_list = NULL,
                                      param_list = NULL,
                                      data_mode = c("namedpairs"),
                                      label_prefix = "C",
                                      train_fraction = 1.0,
                                      seed = 42,
                                      evaluate = TRUE) {

  if (is.null(pair_list)) {
    pair_list <- getFeaturesPairList()
  }

  message("Reading data from: ", data_path)
  data <- readr::read_csv(data_path, show_col_types = FALSE)

  if (!is.null(label_prefix)) {
    data[[label_name]] <- paste0(label_prefix, data[[label_name]])
  }

  # If there's a test set, use that, else split the training.
  if (!is.null(test_path)) {
    message("Reading test data from: ", test_path)
    test_data <- readr::read_csv(test_path, show_col_types = FALSE)
    if (!is.null(label_prefix)) {
      test_data[[label_name]] <- paste0(label_prefix, test_data[[label_name]])
    }
    train_data <- data
  } else if (train_fraction < 1.0) {
    # Split data
    set.seed(seed)
    train_idx <- data %>%
      group_by(.data[[label_name]]) %>%
      sample_frac(train_fraction) %>%
      pull({{sample_id}})

    train_data <- data %>% filter(.data[[sample_id]] %in% train_idx)
    test_data <- data %>% filter(!(.data[[sample_id]] %in% train_idx))

    message(sprintf("Split: %d train, %d test samples",
                    nrow(train_data), nrow(test_data)))
  } else {
    message("Using all data for training (no test split)")
    train_data <- data
    test_data <- NULL
  }

  message("Training robencla classifier...")
  model <- trainRobencla(
    data = train_data,
    label_name = label_name,
    sample_id = sample_id,
    pair_list = pair_list,
    sig_list = sig_list,
    param_list = param_list,
    data_mode = data_mode
  )

  if (evaluate) {
    # Evaluate on test data if available, otherwise training data
    eval_data <- if (!is.null(test_data)) test_data else train_data
    eval_label <- if (!is.null(test_data)) "test" else "training"

    message(sprintf("Evaluating on %s data...", eval_label))
    model$predict(
      data_frame = eval_data,
      label_name = label_name,
      sample_id = sample_id
    )
    results <- model$results()

    cat(sprintf("\nConfusion Matrix (%s data):\n", eval_label))
    print(table(Predicted = results$BestCalls, Actual = model$test_label))

    cat("\nClassification Metrics:\n")
    print(model$classification_metrics())
  }

  if (!is.null(output_path)) {
    message("Saving model to: ", output_path)
    saveRDS(model, file = output_path)
  }

  #invisible(list(model = model, test_data = test_data))
  return(list(model=model, results=results))
}

