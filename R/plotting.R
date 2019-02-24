#' Wrapper around qplot for make cell scatter plots
#' @param emb n-by-2 matrix of cell coordinates, where n is the number of cells
#' @param batch factor or vector of length n
#' @param cell_type factor or vector of length n
#' @export
plot_dimensionality_reduction <- function(emb, batch, cell_type) {
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


#' Return a list of scatter plots
#'
#' @param embeddings a list of cell embeddings, i.e., \code{embeddings[[i]]} is an n-by-2 matrix of coordinates where n is the number of cells
#' @param batch factor or vector of length n
#' @param cell_type factor or vector of length n
#' @return a list of scatter plots corresponding to each embedding in \code{embeddings}
#' @export
get_scatterplots <- function(embeddings, batch, cell_type) {
  lapply(
    X = names(embeddings),
    FUN = function(gs_name) {
      emb <- embeddings[[gs_name]]
      plot_dimensionality_reduction(emb, batch, cell_type) +
        ggtitle(paste0(gs_name, ", ", round(
          corgi::batch_separation(emb, batch), 2
        )))
    }
  )
}


#' Legend for shapes in a scatter plot
#'
#' Useful for creating common shape legend for multiple scatter plots
#' @param batch factor or vector of the batch labels
#' @param my_shape_palette a vector of numbers corresponding to the shape palette http://www.sthda.com/english/wiki/ggplot2-point-shapes
#' @export
get_shape_legend <- function(batch, my_shape_palette) {
  n_cells <- length(batch)
  cowplot::get_legend(
    ggplot2::qplot(1:n_cells, 1:n_cells, shape = batch) +
      scale_shape_manual(values = my_shape_palette) +
      guides(shape = guide_legend(title = "Batch")) +
      theme(legend.title.align = 0.5)
  )
}

#' Legend for colors in a scatter plot
#'
#' Useful for creating common color legend for multiple scatter plots
#' @param cell_type factor or vector of the batch labels
#' @param my_color_palette a vector of colors e.g., \code{c("red","blue","green")} or hex color codes
#' @export
get_color_legend <-
  function(cell_type,
           my_color_palette,
           ncol = NULL,
           legend.position = "right",
           legend.title = "Cell type",
           ...) {
    n_cells <- length(cell_type)
    df <-
      data.frame(x = 1:n_cells,
                 y = 1:n_cells,
                 cell_type = cell_type)

    ggplot(df, aes(x = x, y = y)) + geom_point(aes(color = cell_type), ...) +
      scale_color_manual(values = my_color_palette) +
      guides(col = guide_legend(title = legend.title, ncol = ncol)) +
      theme(legend.position = legend.position) -> plt

    return(cowplot::get_legend(plt))
  }

#' Empty plot with just the axes and labels
#'
#' Useful for creating axes legend for multiple scatter plots with the same axes names
#' @param emb_name name of the embedding, for instance \code{emb_name = "PC"}
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
plot_mapping_accuracy_comparison <- function(results){
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


#' Plot scatterplots for each pair of columns in df1 and df2
#'
#' Code by Ben Bolker, see this \href{https://stackoverflow.com/a/40455040/636276}{stackoverflow answer}
#' @param df1 the first dataframe
#' @param df2 the second dataframe, rownames(df1) == rownames(df2) must hold
#' @author Ben Bolker
#' @examples
#' plot_two_dfs(iris[,1:2],iris[,3:4])
#' @export
plot_two_dfs <- function(df1,df2){
  library(dplyr)
  mfun <- function(x) {
    x %>%
      mutate(obs=seq(n())) %>%    ## add obs numbers
      gather(key=var,value=value,-obs)  ## reshape
  }
  ## combine
  df12 <- mfun(df1) %>% full_join(mfun(df2),by="obs")


  library(ggplot2);
  ggplot(df12,aes(value.x,value.y)) +
    geom_point(size=0.1)+
    facet_grid(var.y~var.x, scales = "free")+
    theme_bw() +
    theme(panel.margin=grid::unit(0,"lines")) ## squash panels together
}
