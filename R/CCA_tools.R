#' @export
Seurat_preprocess <- function(mat1,mat2,batch_name1 = "b1", batch_name2 = "b2") {
  seurat_objects <- list()
  mats <- list(mat1, mat2)
  batch_names <- list(batch_name1, batch_name2)
  for (i in 1:2) {
    seu <- Seurat::CreateSeuratObject(raw.data = mats[[i]])
    seu <- Seurat::NormalizeData(object = seu)
    seu <- Seurat::ScaleData(object = seu, display.progress = F)
    seu@meta.data[, "batch"] <- batch_names[[i]]
    seurat_objects[[i]] <- seu
  }
  names(seurat_objects) <- batch_names
  return(seurat_objects)
}

#' Wrapper around Seurat's CCA to extract gene loadings
#'
#' @param mat1 batch 1 expression matrix, should be raw
#' @param mat2 same as above for batch 2
#' @param k how many CCA components to return, defaults to k
#' @return gene-by-CCA# matrix where # ranges from 1 to k
#' @export
Seurat_CCA_GeneLoading <- function(mat1, mat2, k = 3){
  seurat_objects <- Seurat_preprocess(mat1, mat2)
  Seurat::RunCCA(
    object = seurat_objects[[1]],
    object2 = seurat_objects[[2]],
    genes.use = rownames(mat1),
    num.cc = k
  ) -> combined
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
  gene_expr_combined <- apply(X = cbind(mat1, mat2), MARGIN = 2, FUN = rank)

  cca_combined <- rbind(cca_output$u, cca_output$v)
  return(gene_expr_combined %*% cca_combined)
}
