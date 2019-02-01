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

