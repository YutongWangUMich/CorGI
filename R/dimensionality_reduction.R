#' Multidimensional scaling on Spearman rho correlation distance
#'
#' @param X the input matrix, where each column is a sample (cell) and each row is a feature (gene)
#' @param k number of latent dimensions, default to 2
#' @examples
#' mds <- spearman_rho_mds(t(as.matrix(mtcars)))
#' plot(mds)
#' text(mds, row.names(mtcars), cex=0.6, pos=4, col="red")
#' @export
spearman_rho_mds <- function(X, k = 2){
  D <- sqrt((1-cor(X,method = "spearman"))/2) # Spearman rho distance matrix
  mds <- cmdscale(D,k = k) # Multidimensional-scaling
  return(mds)
}
