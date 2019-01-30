
upper_quantile_genes <- function(res,rng,q){
  n.sv <- ncol(res) - 1
  um <- list()
  # union.upper.medians <- list()
  for(j in rng){
    res[,j]/res[,n.sv+1] -> scores
    names(which(scores >= quantile(scores,q,na.rm=T))) -> um[[j]]
  }
  return(um)
  # return(res)
}



#' Corgi feature selection
#'
#' @param X gene-by-cell expression matrix for batch 1
#' @param Y gene-by-cell expression matrix for batch 2
#' @param n_phases how many phases of successive halving to do
#' @param time_per_phase how long each phase is (in seconds)
#' @param must_have_genes the set of marker genes that is used in every sampled gene set
#' @param n_singular_gaps how many singular value gaps to examine
#' @param n_genes_sample what size gene set to use
#' @param n_cells_X_sample how many cells to sample from X
#' @param n_cells_Y_sample how many cells to sample from Y
#' @return a list of list of scores
#' @export
run_corgi <- function(
  X,
  Y,
  n_phases = 6,
  time_per_phase = 20*60,
  must_have_genes = c(),
  n_singular_gaps = 3,
  n_genes_sample = NULL,
  n_cells_X_sample = NULL,
  n_cells_Y_sample = NULL
  ){
  if(is.null(rownames(X)) | is.null(rownames(Y))){
    stop("Rows of X and Y must have names")
  }

  if(is.null(colnames(X))|is.null(colnames(Y))){
    stop("Columns of X, Y must have names")
  }

  q <- 1/2 # successive halving

  res_list <- list()
  for(i in 1:n_phases){
    start_time <- Sys.time()

    print(start_time)
    res <- run_par(
      X,
      Y,
      run_time = time_per_phase,
      n_singular_gaps = n_singular_gaps,
      must_have_genes = must_have_genes,
      n_genes_sample = n_genes_sample,
      n_cells_X_sample = n_cells_X_sample,
      n_cells_Y_sample = n_cells_Y_sample
    )
    print(nrow(res))
    upper_quantile_genes(res,2:n_singular_gaps,q) -> outliers
    Reduce(union,outliers) -> genes_use
    union(genes_use,must_have_genes) -> genes_use
    res_list[[i]] <- res


    X <- X[genes_use,]
    Y <- Y[genes_use,]

    end_time <- Sys.time()
    print(end_time - start_time)
  }
  return(res_list)
}
