#' The inner loop of corgi
#'
#' @param X gene-by-cell expression matrix for batch 1
#' @param Y gene-by-cell expression matrix for batch 2
#' @param run_time how long to run for
#' @param gene_use the set of genes used for sampling
#' @param must_have_genes the set of marker genes that is used in every sampled gene set
#' @param n_singular_gaps how many singular value gaps to examine
#' @param n_genes_sample what size gene set to use
#' @param n_cells_X_sample how many cells to sample from X
#' @param n_cells_Y_sample how many cells to sample from Y
#' @return a list of scores
#' @import doParallel foreach parallel
#' @export
run_par <- function(
  X,
  Y,
  run_time,
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

  if(is.null(n_cells_Y_sample)){
    n_cells_Y_sample <- min(300, ncol(Y))
  }

  if(is.null(n_genes_sample)){
    n_genes_sample <- min(round(length(genes_use)*0.1),200)
  }

  gene_sample_pool <- setdiff(genes_use,must_have_genes)

  Spearman_rho_Fiedler <- function(X,Y,K){
    rho <- (1+cor(X,Y,method = 'spearman'))/2
    W <- diag(rowSums(rho)^(-1/2)) %*% rho %*% diag(colSums(rho)^(-1/2))
    return(abs(diff(RSpectra::svds(W,k=(K+1))$d)))
  }




  # setup parallel backend to use many processors
  cores = detectCores()
  cl <- makeCluster(cores[1] - 1) # not to overload your computer
  registerDoParallel(cl)
  res <- foreach(
    i = 1:(cores[1] - 1),
    .combine = '+'
  ) %dopar% {
    res_temp <- matrix(0,length(gene_sample_pool),n_singular_gaps+1)

    rownames(res_temp) <- gene_sample_pool
    colnames(res_temp)<-c(paste0("SVG",1:n_singular_gaps), "num.times.sampled")

    start_time <- Sys.time()
    end_time <- start_time

    while(as.numeric(difftime(time1 = end_time,
                              time2 = start_time,
                              units = "secs")) < run_time){
      # print("hi")
      # print(genes.use)
      gset <- sample(gene_sample_pool,n_genes_sample)
      gset_must_have <- union(gset,must_have_genes)
      cset_X <- sample(cells_X,n_cells_X_sample)
      cset_Y <- sample(cells_Y,n_cells_Y_sample)
      # make sure there is no sample with all zeros
      if(min(apply(X[gset_must_have,cset_X],2,max),apply(Y[gset_must_have,cset_Y],2,max))>0){
        res_temp[gset,1:n_singular_gaps]<-
          res_temp[gset,1:n_singular_gaps] +
          pracma::repmat(
            Spearman_rho_Fiedler(X[gset_must_have,cset_X],
                                 Y[gset_must_have,cset_Y],
                                 n_singular_gaps),
                 n_genes_sample,1)

        res_temp[gset,"num.times.sampled"]<-
          res_temp[gset,"num.times.sampled"]+1
      }
      end_time <- Sys.time()
    }

    res_temp #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
    # return(res.temp)
  }
  #
  stopCluster(cl)
  return(res)
}



upper.quantile.genes <- function(res,rng,q){
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
#' @import doParallel foreach
#' @export
corgi <- function(
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

  q <- 1-0.5

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
    upper.quantile.genes(res,2:n_singular_gaps,q) -> outliers
    Reduce(union,outliers) -> genes.use
    res_list[[i]] <- res


    X <- X[genes.use,]
    Y <- Y[genes.use,]

    end_time <- Sys.time()
    print(end_time - start_time)
  }
  return(res_list)
}
