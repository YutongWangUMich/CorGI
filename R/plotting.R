#' @export
plot_dimensionality_reduction <- function(emb,batch,cell_type,gs_name,emb_name){
  qplot(emb[, 1], emb[, 2], color = cell_type, shape = batch) +
    ggtitle(paste0(gs_name, ", ", round(corgi::batch_mixing(emb, batch), 2))) +
    theme_bw() +
    xlab(paste0(emb_name, 1)) +
    ylab(paste0(emb_name, 2)) +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid = element_blank()
    )
}
