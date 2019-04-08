#' Subsample cells random
#'
#' @param object an object that is column-subsettable, e.g., a matrix or a SingleCellExperiment object
#' @param n how many cells do you want
#' @export
sample_cells <- function(object, n){
  return(object[,sample(1:ncol(object),n)])
}

#' Combine several SingleCellExperiment objects
#'
#' @param sce_list a list of SingleCellExperiment objects
#' @param cell_type_use the name of the colData to be used for cell_type, must be available for each member of the list
#' @param ... addition parameters passed to the factors
#' @return a list with three objects: \code{combined} the combined count data, \code{cell_type} concatenated cell types, \code{batch} the batch labels
#' @export
combine_sces <- function(sce_list, cell_type_use = "cell_type1", ...){
  library(SingleCellExperiment)
  do.call(forcats::fct_c,
          lapply(sce_list,
                 function(sce) factor(colData(sce)[[cell_type_use]]))
  ) ->
    cell_type
  cell_type <- factor(cell_type, ...)
  do.call(c,
          lapply(names(sce_list), function(sce_name){
            rep(sce_name, ncol(sce_list[[sce_name]]))
          })
  ) ->
    batch

  batch <- factor(batch, levels = names(sce_list)) # order the levels by the order they were supplied

  shared_genes <- Reduce(intersect, lapply(sce_list, rownames))
  do.call(cbind, lapply(sce_list,
                        function(sce) assay(sce)[shared_genes, ]
  )
  ) ->
    combined
  sce <- SingleCellExperiment(assays = list(counts = combined))
  sce$cell_type <- cell_type
  sce$batch <- batch
  return(sce)
}


#' Convert a list of sets to a membership vector
#'
#' @param set_list a named list of sets
#' @export
get_membership_labels_from_set_list <- function(universe, set_list){
  labels <- rep(NA, length(universe))
  names(labels) <- universe
  for(set_name in names(set_list)){
    labels[universe %in% set_list[[set_name]]] <- set_name
  }
  return(factor(labels))
}
