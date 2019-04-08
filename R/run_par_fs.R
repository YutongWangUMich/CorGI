#' @export
normalized_cor <- function(X, Y, normalization = T){
  rho <- (1 + cor(X, Y, method = 'spearman')) / 2
  if(normalization){
    rho <- t(colSums(rho)^(-1/2) * t(rowSums(rho)^(-1/2) * rho) )
  }
  return(rho)
}

#' The embedding as defined in Dhillon, Inderjit S. "Co-clustering documents and words using bipartite spectral graph partitioning." Proceedings of the seventh ACM SIGKDD international conference on Knowledge discovery and data mining. ACM, 2001.
#' @export
dhillon_emb <- function(X, Y, k){
  R <- corgi::normalized_cor(X, Y, normalization = F)
  left_degrees <- rowSums(R)
  right_degrees <- colSums(R)

  out <- RSpectra::svds(R,k = k)
  emb <- rbind(left_degrees^(-1/2)*out$u[,2:k],
               right_degrees^(-1/2)*out$v[,2:k])
  return(emb)
}



#' The user friendly version of corgi
#'
#' @param X gene-by-cell expression matrix for batch 1
#' @param Y gene-by-cell expression matrix for batch 2
#' @param run_time how long to run for
#' @param must_have_genes the set of marker genes that is used in every sampled gene set
#' @return a list of scores
#' @export

corgi <- function(X,
                  Y,
                  run_time,
                  must_have_genes = c(),
                  k = 3){
  X <- as.matrix(X)
  Y <- as.matrix(Y)


  fs <- list(
    principal_curve = function(X, Y){
      emb <- corgi::dhillon_emb(X, Y, k)
      princurve::principal_curve(emb)$dist
    },
    kmmd_statistic = function(X, Y){
      rk <- function(x) apply(x, 2, rank)
      capture.output(kmmd_out <- kernlab::kmmd(t(rk(X)), t(rk(Y))))
      kmmd_out@mmdstats[1]
    }
  )

  corgi::run_par_fs(
    X = X,
    Y = Y,
    run_time = run_time,
    fs = fs,
    must_have_genes = must_have_genes
  ) -> corgi_out
  return(corgi_out)
}

#' The inner loop of corgi
#'
#' @param X gene-by-cell expression matrix for batch 1
#' @param Y gene-by-cell expression matrix for batch 2
#' @param run_time how long to run for
#' @param fs a list of scoring functions f(X[F, ], Y[F, ]) -> real number
#' @param fs_var global variables, to save time
#' @param gene_use the set of genes used for sampling
#' @param must_have_genes the set of marker genes that is used in every sampled gene set
#' @param n_genes_sample what size gene set to use
#' @param n_cells_X_sample how many cells to sample from X
#' @param n_cells_Y_sample how many cells to sample from Y
#' @return a list of scores
#' @export
run_par_fs <- function(
  X,
  Y,
  run_time,
  fs,
  genes_use = NULL,
  must_have_genes = c(),
  n_genes_sample = NULL,
  n_cells_X_sample = NULL,
  n_cells_Y_sample = NULL,
  n_cores_use = NULL
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

  # setup parallel backend to use many processors
  if(is.null(n_cores_use)){
    cores = parallel::detectCores()
    n_cores_use = cores[1]-1
  }
  cl <- parallel::makeCluster(n_cores_use) # not to overload your computer
  print("number of cores used")
  print(n_cores_use)
  doParallel::registerDoParallel(cl)

  `%dopar%` <- foreach::`%dopar%`

  res <- foreach::foreach(
    i = 1:(n_cores_use),
    .combine = '+'
  ) %dopar% {
    res_temp <- matrix(0, length(fs) + 2, length(gene_sample_pool))

    colnames(res_temp) <- gene_sample_pool
    rownames(res_temp) <- c(names(fs), c("num.times.sampled", "failed"))

    start_time <- Sys.time()
    end_time <- start_time


    while(as.numeric(difftime(time1 = end_time,
                              time2 = start_time,
                              units = "secs")) < run_time){
      gset <- sample(gene_sample_pool,n_genes_sample)
      gset_use <- union(gset,must_have_genes)

      if(do_sample_X){
        X_use <- X[gset_use, sample(cells_X,n_cells_X_sample)]
      }else{
        X_use <- X[gset_use,]
      }

      if(do_sample_Y){
        Y_use <- Y[gset_use,sample(cells_Y,n_cells_Y_sample)]
      }else{
        Y_use <- Y[gset_use,]
      }
      # compute shared global variables
      if(min(apply(X_use,2,max),apply(Y_use,2,max))>0){
        for(fname in names(fs)){
          f <- fs[[fname]]
          res_temp[fname, gset]<-
            res_temp[fname, gset] +
            f(X_use, Y_use)
        }
        res_temp["num.times.sampled",gset] <- res_temp["num.times.sampled",gset] + 1
      }else{
        res_temp["failed", gset] <- res_temp["failed", gset] + 1
      }


      end_time <- Sys.time()
    }

    res_temp #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
  }
  #
  parallel::stopCluster(cl)
  return(t(res))
}

