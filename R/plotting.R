#' Wrapper around qplot for making cell scatter plots
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
  library(ggplot2)
  comparison_legend_options <- guide_legend(keywidth = 2, keyheight = 1, title = "Gene set")
  ggplot(results, aes(x=Param, y=Kappa, group=Gene_set)) +
    geom_line(aes(linetype = Gene_set)) +
    geom_point(aes(shape = Gene_set))+
    guides(linetype = comparison_legend_options,
           shape = comparison_legend_options) +
    theme_bw() +
    scale_x_continuous(breaks = unique(results$Param)) +
    theme(panel.grid.minor.x = element_blank()) +
    xlab("scmap Parameter") +
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


#' Plot a venn diagram
#'
#' Code obtained from
#' https://scriptsandstatistics.wordpress.com/2018/04/26/how-to-plot-venn-diagrams-using-r-ggplot2-and-ggforce/
#' @param my_sets a named list of three sets
#' @examples plot_venn_diagram(list(A = 1:5, B = 3:7, C = 5:9))
#' @export
plot_venn_diagram <- function(my_sets) {
  library(ggforce)
  library(limma)
  library(dplyr)
  df.venn <- data.frame(
    x = c(0, 0.866,-0.866),
    y = c(1,-0.5,-0.5),
    labels = names(my_sets)
  )

  all_genes <- Reduce(union, my_sets)
  names(my_sets)[1:3]
  mydata <-
    data.frame(all_genes %in% my_sets[[1]],
               all_genes %in% my_sets[[2]],
               all_genes %in% my_sets[[3]])
  colnames(mydata) <- names(my_sets)

  vdc <- vennCounts(mydata)
  class(vdc) <- 'matrix'
  df.vdc <- as.data.frame(vdc)[-1, ] %>%
    mutate(
      x = c(0, 1.2, 0.8,-1.2,-0.8, 0, 0),
      y = c(1.2,-0.6, 0.5,-0.6, 0.5,-1, 0)
    )

  df.gene_set_names <-
    df.vdc[c(1, 2, 4), 4:6]
  df.gene_set_names

  df.gene_set_names$Name <- rev(names(my_sets))
  df.gene_set_names$y[1] <- df.gene_set_names$y[1] + 0.5
  df.gene_set_names$y[2:3] <- df.gene_set_names$y[2:3] - 0.5
  df.gene_set_names$x[2:3] <- df.gene_set_names$x[2:3] * 1.25


  ggplot(df.venn) +
    geom_circle(
      aes(
        x0 = x,
        y0 = y,
        r = 1.5,
        fill = labels
      ),
      alpha = .3,
      size = 1,
      colour = 'grey'
    ) +
    coord_fixed() +
    theme_void() +
    theme(legend.position = 'none') +
    scale_fill_manual(values = c('cornflowerblue', 'firebrick',  'gold')) +
    scale_colour_manual(values = c('cornflowerblue', 'firebrick', 'gold'),
                        guide = FALSE) +
    labs(fill = NULL) +
    annotate(
      "text",
      x = df.vdc$x,
      y = df.vdc$y,
      label = df.vdc$Counts,
      size = 5
    ) +
    annotate("text",
             x = df.gene_set_names$x,
             y = df.gene_set_names$y,
             label = df.gene_set_names$Name) ->
    venn_diagram

  venn_diagram
}


#' Plot gene-by-cell heatmap of a matrix
#'
#'
#' @param mat gene (row) by cell (column) matrix
#' @param annotation_col a data frame of meta data for the cells
#' @param order_cells_by which column of annotation_col to use
#' @param annotation_row a data frame of meta data for the genes
#' @param cluster_rows whether or not to cluster the rows, see pheatmap
#' @param cutree_rows number of clusters to be returned by the hiearchical clustering
#' @param annotate_gene_clusters boolean value of whether the gene cluster labels should be shown in the row annotation. This is because you might want to cluster the genes but not show the colors.
#' @param annotation_colors a list of color palettes (color hex codes), you don't have to supply all or even any colors. By default, the rainbow palette is applied to discrete color scales. For continuous data, black/white gradient is applied.
#' @export
plot_gene_by_cell_heatmap <- function(
  mat,
  annotation_col,
  order_cells_by,
  annotation_row = NULL,
  cluster_rows = T,
  cutree_rows = 7,
  annotate_gene_clusters = T,
  annotation_colors = list(),
  cluster_cols = F,
  ...
){
  cell_type <- annotation_col[[order_cells_by]]

  if(is.null(annotation_row) & cluster_rows & annotate_gene_clusters){
    annotation_row <- data.frame(row.names = rownames(mat))
  }
  if(cluster_rows & annotate_gene_clusters){

    pheatmap::pheatmap(mat[,order(cell_type)], silent = T) -> result

    gene_cluster <- factor(cutree(tree = result$tree_row,k = cutree_rows))
    annotation_row[["gene_cluster"]] <- gene_cluster
  }

  # convert logical meta data into factors so that it can be plotted
  for(md_name in names(annotation_col)){
    md <- annotation_col[[md_name]]
    if(is.logical(md)){
      annotation_col[[md_name]] <- factor(md)
      md <- annotation_col[[md_name]]
      if(!(md_name %in% names(annotation_colors))){
        annotation_colors[[md_name]] <- c("white", "black")
        names(annotation_colors[[md_name]]) <- c("FALSE", "TRUE")
      }
    }
  }

  for(md_name in names(annotation_row)){
    md <- annotation_row[[md_name]]
    if(is.logical(md)){
      annotation_row[[md_name]] <- factor(md)
      md <- annotation_row[[md_name]]
      if(!(md_name %in% names(annotation_colors))){
        annotation_colors[[md_name]] <- c("gray80", "black")
        names(annotation_colors[[md_name]]) <- c("FALSE", "TRUE")
      }
    }
  }

  # add the colors for the column
  for(ann in list(annotation_col, annotation_row)){
    for(md_name in setdiff(colnames(ann), names(annotation_colors))){
      md <- ann[[md_name]]
      if(is.numeric(md)){
        annotation_colors[[md_name]] <- c("white", "black")
      }else{
        my_levels <- levels(droplevels(md))
        annotation_colors[[md_name]] <- rainbow(length(my_levels))
        names(annotation_colors[[md_name]]) <- my_levels
      }
    }
  }

  # add the colors for shared meta data between genes and clusters
  # you might consider having shared meta data when certain genes are markers for certain clusters
  shared_meta_data_names <- intersect(colnames(annotation_col),
                                      colnames(annotation_row))
  for(md_name in setdiff(shared_meta_data_names, names(annotation_colors))){
    if(!is.factor(annotation_col[[md_name]]) |
       !is.factor(annotation_row[[md_name]])){
      stop("shared meta data must both be factors")
    }
    md <- forcats::fct_c()

    if(md_name %in% colnames(annotation_col)){
      md <- forcats::fct_c(md, annotation_col[[md_name]])
    }
    if(md_name %in% colnames(annotation_row)){
      md <- forcats::fct_c(md, annotation_row[[md_name]])
    }
    my_levels <- levels(droplevels(md))
    annotation_colors[[md_name]] <- rainbow(length(my_levels))
    names(annotation_colors[[md_name]]) <- my_levels
  }


  args <- list(
    mat[, order(cell_type)],
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    ...
  )

  if(is.data.frame(annotation_row) & (ncol(annotation_row)>0)){
    args[["annotation_row"]] <- annotation_row
  }

  if(is.factor(cell_type)) {
    args[["gap_cols"]] <-
      which(c(diff(as.numeric(cell_type[order(cell_type)])), 0) == 1)
  }
  do.call(pheatmap::pheatmap, args)
}

