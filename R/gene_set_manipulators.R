#' Highly dropped-out gene selection
#'
#' @param sce SingleCellExperiment object to be normalized
#' @param do_normalize boolean value of whether \code{sce} should be normalized prior to feature selection. Consider setting this to false if you prefer a different normalization technique.
#' @return ranked list such that most informative genes are at the top of the list
#' @export
HDG_ranking <- function(sce, do_normalize = T){
  library(SingleCellExperiment)
  if(do_normalize){
    sce <- scater::normalize(sce)
  }
  rowData(sce)$feature_symbol <- rownames(sce)
  sce <- scmap::selectFeatures(sce,suppress_plot = F)
  return(rownames(sce)[order(rowData(sce)[["scmap_scores"]],decreasing = T,na.last = T)])
}

#' Get compared gene sets
#'
#' Gene sets compared against CorGI in the manuscript
#' @param batch1_top_genes ordered list of genes (e.g., most highly variable to least)
#' @param batch1_name name of batch1 (e.g., if \code{batch1_top_genes} was computed using highly variable genes (HVG), then consider \code{batch1_name <- HVG(Batch1)})
#' @param batch2_top_genes see \code{batch1_top_genes}
#' @param batch2_name see \code{batch1_name}
#' @param desired_size how big the output gene_sets should be
#' @param marker_genes genes that are included in every gene set of the output. Defaults to the empty set, i.e., no markers provided.
#' @return a list of gene sets, where each item is a gene set of size \code{desired_size}. The list is ordered as follows:
#' \itemize{
#'   \item \code{batch1_top_genes[1:desired_size]}
#'   \item \code{batch2_top_genes[1:desired_size]}
#'   \item \code{union(batch1_top_genes[1:x], batch2_top_genes[1:y])} where \code{x} and \code{y} satisfies \code{|x-y| <= 1}
#'   \item \code{intersect(batch1_top_genes[1:x],batch2_top_genes[1:y])} where \code{x} and \code{y} satisfies \code{|x-y| <= 1}
#' }
#'   If \code{marker_genes} is nonempty, then all four gene sets above include \code{marker_genes}.
#' @examples
#' HVG_batch1 <- c(paste0("g",1:4),paste0("g",7:10))
#' HVG_batch2 <- paste0("g",3:9)
#' marker_genes <- "g0"
#' get_compared_gene_sets(HVG_batch1, "HVG(batch1)", HVG_batch2, "HVG(batch2)", 6, marker_genes)
#' @export
get_compared_gene_sets <- function(batch1_top_genes,
                                   batch1_name = "Batch1",
                                   batch2_top_genes,
                                   batch2_name = "Batch2",
                                   desired_size,
                                   marker_genes = c()){

  gene_sets <- list()
  n_genes <- min(length(batch1_top_genes),length(batch2_top_genes))

  # remove marker_genes from the ranked lists
  batch1_top_genes <- setdiff(batch1_top_genes,marker_genes)
  batch2_top_genes <- setdiff(batch2_top_genes,marker_genes)

  gene_sets[[batch1_name]] <- c(batch1_top_genes[1:(desired_size - length(marker_genes))], marker_genes)
  gene_sets[[batch2_name]] <- c(batch2_top_genes[1:(desired_size - length(marker_genes))], marker_genes)

  union_gene_set <- function(x,y){c(union(batch1_top_genes[1:x], batch2_top_genes[1:y]), marker_genes)}


  # gs_func is a function that takes 2 integer inputs, x and y, and output a gene set
  # Requirements:
  # 0 <= gs_func(x+1,y) - gs_func(x,y) <= 1
  # 0 <= gs_func(x,y+1) - gs_func(x,y) <= 1

  get_gs_from_gs_func <- function(gs_func){
    z <- 1
    while(length(gs_func(z,z)) < desired_size){
      z <- z+1
    }
    if(length(gs_func(z,z)) == desired_size){
      return(gs_func(z,z))
    }

    x <- z - 1
    y <- z - 1
    i <- 1
    while(length(gs_func(x,y)) < desired_size){
      if((i%%2) == 1){
        x <- x + 1
      }else{
        y <- y + 1
      }
      i <- i + 1
      if(x > length(batch1_top_genes) & y > length(batch2_top_genes)){
        stop("Not enough genes in batch1_top_genes or batch2_top_genes")
      }
    }
    print(list(x = x, y = y))
    return(gs_func(x,y))
  }

  gs_func_list <- list()
  gs_func_list[["Union"]] <-
    function(x, y){c(union(batch1_top_genes[1:x], batch2_top_genes[1:y]), marker_genes)}
  gs_func_list[["Intersection"]] <-
    function(x, y){c(intersect(batch1_top_genes[1:x], batch2_top_genes[1:y]), marker_genes)}

  for(gs_name in names(gs_func_list)){
    gene_sets[[gs_name]] <- get_gs_from_gs_func(gs_func_list[[gs_name]])
  }

  return(gene_sets)
}


#' @param x named numeric vector where the names are the gene names and the numbers are gene scores
#' @param n how many genes to return
#' @export
top_n_genes <- function(x,n) {
  return(
    head(names(sort(x,decreasing = T)),n)
  )
}
