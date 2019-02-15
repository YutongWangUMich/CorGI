#' @export
spearman_rho_mds <- function(X){
  D <- sqrt((1-cor(X,method = "spearman"))/2) # Spearman rho distance matrix
  mds <- cmdscale(D,k=2) # Multidimensional-scaling
  return(mds)
}
