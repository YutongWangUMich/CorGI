#' Compute the batch mixing metric
#' @param emb the embedding of the cells where each row is a cell
#' @param batch the batch labels of the cells
#' @export
batch_mixing <- function(emb,batch){
  dat <- data.frame(emb,batch)
  model <- e1071::svm(batch~.,data=dat)
  return(1- caret::confusionMatrix(model$fitted,dat$batch)$overall[["Accuracy"]])
}



#' @export
cluster_coherence <- function(emb, cell_type, cell_type_pred, train, test, knn_params = 5*(1:6)){
  results <- data.frame(
    knn_param = numeric(0),
    TPR = numeric(0),
    FPR = numeric(0),
    stringsAsFactors = F
  )
  i <- 1
  for (k in knn_params) {
    results[i, "knn_param"] <- k
    class::knn(
      train = emb[train,],
      test = emb[test,],
      cl = cell_type[train],
      k = k
    ) -> cell_type_knn

    caret::confusionMatrix(data = cell_type_knn,
                           reference = cell_type[test]) -> output

    results[i, "TPR"] <- output$byClass[paste("Class:", cell_type_pred), "Sensitivity"]
    results[i, "FPR"] <- 1 - output$byClass[paste("Class:", cell_type_pred), "Specificity"]
    i <- i + 1
  }
  return(results)
}

