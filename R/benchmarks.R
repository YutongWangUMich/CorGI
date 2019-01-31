#' Compute the batch mixing metric
#' @param emb the embedding of the cells where each row is a cell
#' @param Batch the batch labels of the cells
#' @export
batch_mixing <- function(emb,Batch){
  dat <- data.frame(emb,Batch)
  m <- e1071::svm(Batch~.,data=dat)
  tab <- table(m$fitted,Batch)
  tot <- sum(tab)
  error <- (tot - sum(diag(tab)))/tot
  return(error)
}




