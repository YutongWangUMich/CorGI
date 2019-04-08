#' Compute the batch separation metric
#'
#' @param emb the embedding of the cells where each row is a cell, and each column is a feature (e.g. principal component)
#' @param batch the batch labels of the cells
#' @return Cohen's kappa of the the predicted batch labels with the true batch label on the entire dataset (i.e., training accuracy)
#' @examples
#' # suppose that iris[["Species"]] are the batch labels
#' batch_separation(emb = iris[,1:2],batch = iris$Species) # lower batch separation present in dim = 1:2
#' plot(iris[,1],iris[,2],col = iris[["Species"]])
#'
#' batch_separation(emb = iris[,3:4],batch = iris$Species) # higher batch separation in dim = 3:4
#' plot(iris[,3],iris[,4],col = iris[["Species"]])
#' @export
batch_separation <- function(emb,batch){
  dat <- data.frame(emb, batch)
  model <- e1071::svm(batch ~ ., data = dat)
  return(caret::confusionMatrix(model$fitted, dat$batch)$overall[["Kappa"]])
}


#' Wrapper for \code{scmap::scmapCluster}
#' @param query the query dataset (SingleCellExperiment object)
#' @param ref the reference dataset (SingleCellExperiment object)
#' @param gene_set a list of characters representing the genes
#' @param threshold the threshold parameter used in scmapCluster, between 0 and 1
#' @param ... additional parameters passed on to \code{caret::confusionMatrix}. For example \code{dnn = c("Prediction","Truth")}.
#' @return confusion matrix of the true labels in `ref` and the `scmapCluster` predicted labels
#' @export
run_scmapCluster <- function(query, ref, gene_set, param, ...){

  # scmap requires the "feature_symbol" slot to be filled
  rowData(query)$feature_symbol <- rownames(query)
  rowData(ref)$feature_symbol <- rownames(ref)

  # Specify the gene set used for scmapCluster
  rowData(ref)$scmap_features <-
    rowData(ref)$feature_symbol %in% gene_set

  ref <- scmap::indexCluster(ref)

  scmapCluster_results <- scmap::scmapCluster(
    projection = query,
    index_list = list(Reference = metadata(ref)$scmap_cluster_index),
    threshold = param
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


#' Wrapper for \code{scmap::scmapCell}
#' @param query the query dataset (SingleCellExperiment object)
#' @param ref the reference dataset (SingleCellExperiment object)
#' @param gene_set a list of characters representing the genes
#' @param param the number of nearest neighbors in scmapCell, an integer
#' @param ... additional parameters passed on to \code{caret::confusionMatrix}. For example \code{dnn = c("Prediction","Truth")}.
#' @return confusion matrix of the true labels in `ref` and the `scmapCluster` predicted labels
#' @export
run_scmapCell <- function(query, ref, gene_set, param, ...) {

  # scmap requires the "feature_symbol" slot to be filled
  rowData(query)$feature_symbol <- rownames(query)
  rowData(ref)$feature_symbol <- rownames(ref)

  # Specify the gene set used for scmapCluster
  rowData(ref)$scmap_features <-
    rowData(ref)$feature_symbol %in% gene_set

  ref <- indexCell(ref)

  scmapCell_results <- scmapCell(
    projection = query,
    index_list = list(Reference = metadata(ref)$scmap_cell_index),
    w = param)

  scmapCell_clusters <- scmapCell2Cluster(scmapCell_results,
                                          list(as.character(colData(ref)$cell_type1)))


  true_labels <- colData(query)$cell_type1
  pred_labels <-
    factor(scmapCell_clusters$scmap_cluster_labs[, "Reference"])

  # make the factor levels uniform
  shared_levels <- union(levels(true_labels),
                         levels(pred_labels))
  true_labels <- forcats::fct_expand(true_labels, shared_levels)
  pred_labels <- forcats::fct_expand(pred_labels, shared_levels)
  caret::confusionMatrix(data = pred_labels,
                         reference = true_labels,
                         ...)
}




#' Run scmap over a range of parameters, over multiple gene sets
#' @param query a SingleCellExperiment object for the query dataset
#' @param reference a SingleCellExperiment object for the reference dataset
#' @param gene_sets a named list of gene sets
#' @param method_use a string denoting whether to use scmapCluster or scmapCell
#' @export
run_mapping_accuracy_comparison <-
  function(query,
           reference,
           gene_sets,
           method_use = "scmapCluster") {

  f <- NULL
  if(method_use == "scmapCluster"){
    params <- 0.1*(1:9)
    f <- corgi::run_scmapCluster
  }
  if(method_use == "scmapCell"){
    params <- 5*(1:6)
    f <- corgi::run_scmapCell
  }
  lapply(
    X = params,
    FUN = function(param) {
      gene_sets %>%
        lapply(
          FUN = function(gene_set) {
            f(
              query = query,
              ref = reference,
              gene_set = gene_set,
              param = param
            )
          }
        ) %>%
        lapply(
          FUN = function(confusionMat) {
            confusionMat$overall
          }
        ) %>%
        Reduce(f = rbind) %>%
        data.frame ->
        results

      results$Gene_set <- names(gene_sets)
      results$Param <- param
      rownames(results) <- NULL
      results
    }) %>%
    Reduce(f = rbind) -> results

  results$Gene_set <-
    factor(results$Gene_set,
           levels = names(gene_sets))
  return(results)
}
