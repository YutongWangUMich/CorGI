#' @export
plot_dimensionality_reduction <- function(emb,batch,cell_type){
  qplot(emb[, 1], emb[, 2], color = cell_type, shape = batch) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none"
    )
}

#' @export
get_shape_legend <- function(batch,my_shape_palette){
  n_cells <- length(batch)
  cowplot::get_legend(
    ggplot2::qplot(1:n_cells, 1:n_cells, shape = batch) +
      scale_shape_manual(values = my_shape_palette) +
      guides(shape=guide_legend(title="Batch")) +
      theme(legend.title.align=0.5)
  )
}

#' @export
get_color_legend <- function(cell_type,my_color_palette,ncol = NULL,legend.position = "right",legend.title = "Cell type"){
  n_cells <- length(cell_type)
  cowplot::get_legend(
    ggplot2::qplot(1:n_cells, 1:n_cells, color = cell_type) +
      scale_color_manual(values = my_color_palette) +
      guides(col = guide_legend(title= legend.title, ncol = ncol)) +
      theme(legend.position = legend.position)
    )
}

#' Returns an empty plot with just the axes and the axes labels
#' @export
get_axes_legend <- function(emb_name){
  qplot(iris[, 1], iris[, 2], asp = 1, shape = NA, na.rm = TRUE) +
    xlab(paste0(emb_name, 1)) +
    ylab(paste0(emb_name, 2)) +
    theme(
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      text = element_text(size = 10)
    )
}

#' @export
get_scatterplots <- function(embeddings, batch, cell_type){
  lapply(X = names(gene_sets),
         FUN = function(gs_name) {
           emb <- embeddings[[gs_name]]
           plot_dimensionality_reduction(emb, batch, cell_type) +
             ggtitle(paste0(gs_name, ", ", round(corgi::batch_mixing(emb, batch), 2)))
         })
}

#' @export
get_AUC <- function(embeddings, cell_type, cell_type_pred, train, test){
  knn_results <-
    Reduce(rbind,lapply(X = names(embeddings),
                        FUN = function(gs_name){
                          emb <- embeddings[[gs_name]]
                          results <- corgi::cluster_coherence(emb, cell_type, cell_type_pred, train, test)
                          results$Gene_set <- gs_name
                          return(results)
                        }))
  return(knn_results)
}


#' @export
run_cluster_coherence_comparison <- function(query, reference){
  thresholds <- 0.1*(1:9)
  lapply(
    X = thresholds,
    FUN = function(threshold) {
      gene_sets %>%
        lapply(
          FUN = function(gene_set) {
            corgi::run_scmap(
              query = query,
              ref = reference,
              gene_set = gene_set,
              threshold = threshold
            )
          }
        ) %>%
        lapply(
          FUN = function(confusionMat) {
            confusionMat$overall
          }
        ) %>%
        Reduce(f = rbind) %>%
        data.frame ->
        results

      results$Gene_set <- names(gene_sets)
      results$Threshold <- threshold
      rownames(results) <- NULL
      results
    }) %>%
    Reduce(f = rbind) -> results

  results$Gene_set <-
    factor(results$Gene_set,
           levels = names(gene_sets))
  return(results)
}

#' @export
plot_cluster_coherence_comparison <- function(results){
  comparison_legend_options <- guide_legend(keywidth = 2, keyheight = 1, title = "Gene set")
  ggplot(results, aes(x=Threshold, y=Kappa, group=Gene_set)) +
    geom_line(aes(linetype = Gene_set)) +
    geom_point(aes(shape = Gene_set))+
    guides(linetype = comparison_legend_options,
           shape = comparison_legend_options) +
    theme_bw() +
    scale_x_continuous(breaks = 0.1*c(1,3,5,7,9)) +
    xlab("scmapCluster threshold") +
    ylab("Cohen's Kappa")
}
