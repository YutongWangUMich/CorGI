#' Compute the batch separation metric
#'
#' @param emb the embedding of the cells where each row is a cell, and each column is a feature (e.g. principal component)
#' @param batch the batch labels of the cells
#' @return Cohen's kappa of the the predicted batch labels with the true batch label on the entire dataset (i.e., training accuracy)
#' @export
batch_separation <- function(emb,batch){
  dat <- data.frame(emb, batch)
  model <- e1071::svm(batch ~ ., data = dat)
  return(caret::confusionMatrix(model$fitted, dat$batch)$overall[["Kappa"]])
}


#' Wrapper around the scmapCluster function from scmap
#' @param query the query dataset (SingleCellExperiment object)
#' @param ref the reference dataset (SingleCellExperiment object)
#' @param gene_set a list of characters representing the genes
#' @param threshold the threshold parameter used in scmapCluster, between 0 and 1
#' @return confusion matrix of the true labels in `ref` and the `scmapCluster` predicted labels
#' @export
run_scmap <- function(query, ref, gene_set, threshold, ...){
  # Specify the gene set used for scmapCluster
  rowData(ref)$scmap_features <-
    rowData(ref)$feature_symbol %in% gene_set

  ref <- scmap::indexCluster(ref)

  scmapCluster_results <- scmap::scmapCluster(
    projection = query,
    index_list = list(Reference = metadata(ref)$scmap_cluster_index),
    threshold = threshold
  )

  true_labels <- colData(query)$cell_type1
  pred_labels <-
    factor(scmapCluster_results$scmap_cluster_labs[, "Reference"])

  # make the factor levels uniform
  shared_levels <- union(levels(true_labels),
                         levels(pred_labels))
  true_labels <- forcats::fct_expand(true_labels, shared_levels)
  pred_labels <- forcats::fct_expand(pred_labels, shared_levels)

  caret::confusionMatrix(data = pred_labels,
                         reference = true_labels,
                         ...)
}
