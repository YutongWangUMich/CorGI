#' Get gene set of size \code{n} from arbitrary filter \code{gs_func}
#'
#' Binary search algorithm for finding gene sets of a desired size
#' @param gs_func monotonic nondecreasing function in one (integer-valued) input such that \code{gs_func(x)} is a gene set (represented by a vector)
#' @param n number of genes desired
#' @param LB lower bound for the input of \code{gs_func}
#' @param UB upper bound for the input of \code{gs_func}
#' @section Requirement:
#' Monotonic nondecreasing means that whenever \code{x1 >= x2}, \code{length(gs_func(x1)) >= length(gs_func(x2))}.
#' \code{length(gs_func(LB)) <= n <= length(gs_func(UB))} must hold
#' @return a gene set (a vector of strings)
#' @export
get_gene_set_of_size_n <- function(gs_func,n,LB,UB){
  f <- function(x) length(gs_func(x))
  x_test <- floor((UB-LB)/2)
  n_test <- f(x_test)
  while((n_test != n) & !(x_test %in% c(LB,UB))){
    if(n_test < n){
      LB <- x_test
      x_test <- ceiling((UB+x_test)/2)
    }else{
      UB <- x_test
      x_test <- floor((x_test+LB)/2)
    }
    n_test <- f(x_test)
  }
  return(gs_func(x_test))
}


#' @param batch1_top_genes ordered list of genes (e.g., most highly variable to least)
#' @param batch1_name name of batch1 (e.g., if \code{batch1_top_genes} was computed using highly variable genes (HVG), then consider \code{batch1_name <- HVG(Batch1)})
#' @param batch2_top_genes see above
#' @param batch2_name see above
#' @param desired_size how big the output gene_sets should be
#' @param marker_genes genes that are included in every gene set of the output. Defaults to the empty set, i.e., no markers provided.
#' @return a list of gene sets
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

  # gs_func_list <- list()
  # gs_func_list[["Union"]] <-
  #   function(x){c(union(batch1_top_genes[1:x], batch2_top_genes[1:x]), marker_genes)}
  # gs_func_list[["Intersection"]] <-
  #   function(x){c(intersect(batch1_top_genes[1:x], batch2_top_genes[1:x]), marker_genes)}
  #
  # for(gs_name in names(gs_func_list)){
  #   gene_sets[[gs_name]] <- get_gene_set_of_size_n(
  #     gs_func = gs_func_list[[gs_name]],
  #     n = desired_size,
  #     LB = 1,
  #     UB = n_genes)
  # }
  return(gene_sets)
}
