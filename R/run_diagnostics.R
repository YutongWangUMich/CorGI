#' The inner loop of corgi but with more time
#'
#' @param X gene-by-cell expression matrix for batch 1
#' @param Y gene-by-cell expression matrix for batch 2
#' @param run_time how long to run for
#' @param record_time_interval how long between each record
#' @param gene_use the set of genes used for sampling
#' @param must_have_genes the set of marker genes that is used in every sampled gene set
#' @param n_singular_gaps how many singular value gaps to examine
#' @param n_genes_sample what size gene set to use
#' @param n_cells_X_sample how many cells to sample from X
#' @param n_cells_Y_sample how many cells to sample from Y
#' @return a list of scores
run_diagnostics <- function(
  X,
  Y,
  run_time,
  record_time_interval,
  genes_use = NULL,
  must_have_genes = c(),
  n_singular_gaps = 3,
  n_genes_sample = NULL,
  n_cells_X_sample = NULL,
  n_cells_Y_sample = NULL
){

  # make sure all genes are the same
  stopifnot(all(rownames(X)==rownames(Y)))

  # make sure X and Y are non-negative
  X <- X - min(X)
  Y <- Y - min(Y)

  cells_X <- colnames(X)
  cells_Y <- colnames(Y)

  if(is.null(genes_use)){
    genes_use <- rownames(X)
  }

  if(is.null(n_cells_X_sample)){
    n_cells_X_sample <- min(300, ncol(X))
  }
  do_sample_X <- (n_cells_X_sample < ncol(X))

  if(is.null(n_cells_Y_sample)){
    n_cells_Y_sample <- min(300, ncol(Y))
  }
  do_sample_Y <- (n_cells_Y_sample < ncol(Y))

  if(is.null(n_genes_sample)){
    n_genes_sample <- min(round(length(genes_use)*0.1),200)
  }

  gene_sample_pool <- setdiff(genes_use,must_have_genes)

  Spearman_rho_singular_values <- function(X,Y,K){
    rho <- (1+cor(X,Y,method = 'spearman'))/2
    # W <- diag(rowSums(rho)^(-1/2)) %*% rho %*% diag(colSums(rho)^(-1/2))
    W <- t(colSums(rho)^(-1/2) * t(rowSums(rho)^(-1/2) * rho) )
    return(-diff(RSpectra::svds(W,k=(K+1),nu = 0, nv = 0)[["d"]]))
  }


  res_list <- list()

  res_temp <- matrix(0, n_singular_gaps + 1, length(gene_sample_pool))

  colnames(res_temp) <- gene_sample_pool
  rownames(res_temp) <-
    c(paste0("SVG", 1:n_singular_gaps), "num.times.sampled")

  start_time <- Sys.time()
  end_time <- start_time
  time_since_last_record <- start_time
  i <- 1
  while (as.numeric(difftime(
    time1 = end_time,
    time2 = start_time,
    units = "secs"
  )) < run_time) {
    gset <- sample(gene_sample_pool, n_genes_sample)
    gset_use <- union(gset, must_have_genes)

    if (do_sample_X) {
      X_use <- X[gset_use, sample(cells_X, n_cells_X_sample)]
    } else{
      X_use <- X[gset_use, ]
    }

    if (do_sample_Y) {
      Y_use <- Y[gset_use, sample(cells_Y, n_cells_Y_sample)]
    } else{
      Y_use <- Y[gset_use, ]
    }
    # make sure there is no sample with all zeros
    if (min(apply(X_use, 2, max), apply(Y_use, 2, max)) > 0) {
      res_temp[1:n_singular_gaps, gset] <-
        res_temp[1:n_singular_gaps, gset] +
        Spearman_rho_singular_values(X_use, Y_use, n_singular_gaps)

      res_temp["num.times.sampled", gset] <-
        res_temp["num.times.sampled", gset] + 1
    }
    end_time <- Sys.time()
    if (as.numeric(difftime(
      time1 = end_time,
      time2 = time_since_last_record,
      units = "secs"
    )) > record_time_interval) {
      print(paste0("Record number: ",i))
      time_since_last_record <- end_time
      res_list[[i]] <- t(res_temp)
      i <- i+1
    }

  }


  # return(res.temp)


  return(t(res_list))
}

