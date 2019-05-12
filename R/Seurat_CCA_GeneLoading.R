#' Wrapper around Seurat's CCA to extract gene loadings
#'
#' @param mat1 batch 1 expression matrix, should be raw
#' @param mat2 same as above for batch 2
#' @param k how many CCA components to return, defaults to k
#' @return gene-by-CCA# matrix where # ranges from 1 to k
#' @export
Seurat_CCA_GeneLoading <- function(mat1, mat2){

  Seurat_preprocess <- function(mat1,mat2) {
    seurat_objects <- list()
    mats <- list(mat1, mat2)
    names(mats) <- c("x","y")
    for (i in c("x","y")) {
      seu <- Seurat::CreateSeuratObject(raw.data = mats[[i]])
      seu <- Seurat::NormalizeData(object = seu)
      seu <- Seurat::ScaleData(object = seu, display.progress = F)
      seu <- Seurat::FindVariableGenes(object = seu, do.plot = FALSE)
      seu@meta.data[, "batch"] <- i
      seurat_objects[[i]] <- seu
    }

    return(seurat_objects)
  }
  seurat_objects <- Seurat_preprocess(mat1, mat2)

  get_hvg <- function(seu){
    rownames(x = head(x = seu@hvg.info, n = 2000))
  }

  hvg.union <- do.call(union,lapply(seurat_objects, get_hvg))

  Seurat::RunCCA(
    object = seurat_objects[[1]],
    object2 = seurat_objects[[2]],
    genes.use = hvg.union
  ) -> combined
  return(combined@dr$cca@gene.loadings)
}

