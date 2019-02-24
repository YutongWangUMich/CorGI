#' @export
Seurat_preprocess <- function(mat1,mat2,batch_name1 = "b1", batch_name2 = "b2") {
  seurat_objects <- list()
  mats <- list(mat1, mat2)
  batch_names <- list(batch_name1, batch_name2)
  for (i in 1:2) {
    seu <- CreateSeuratObject(raw.data = mats[[i]])
    seu <- NormalizeData(object = seu)
    seu <- ScaleData(object = seu, display.progress = F)
    seu@meta.data[, "batch"] <- batch_names[[i]]
    seurat_objects[[i]] <- seu
  }
  names(seurat_objects) <- batch_names
  return(seurat_objects)
}

#' @export
Seurat_CCA_GeneLoading <- function(mat1, mat2, k = 3){
  seurat_objects <- Seurat_preprocess(mat1, mat2)
  RunCCA(object = seurat_objects[[1]],
         object2 = seurat_objects[[2]],
         genes.use = rownames(mat1),
         num.cc = k) -> combined
  return(combined@dr$cca@gene.loadings)
}


#' @export
CanonCorSpearman <- function(mat1, mat2, k = 3){
  mat3 <- cor(x = mat1, y = mat2, method = "spearman")
  RSpectra::svds(mat3, k = k)
}

#' @export
CCSpearman_GeneLoading <- function(mat1, mat2, k = 3){
  cca_output <- CanonCorSpearman(mat1, mat2, k = k)
  combined <- cbind(mat1, mat2)

  cca_combined <- rbind(cca_output$u, cca_output$v)
  return(apply(X = combined, MARGIN = 2, FUN = rank) %*% cca_combined)
}

#' @export
top_n_genes <- function(x,n) head(names(sort(x,decreasing = T)),n)
