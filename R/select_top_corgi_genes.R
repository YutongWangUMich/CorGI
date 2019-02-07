#' Get the gene set after running corgi
#'
#' @param corgi_output the output of run_corgi
#' @param n number of genes in the returned gene set
#' @param phase_use which phase of the corgi_output to use, defaults to the last phase
#' @param SVG_use which singular value gap (SVG) to use, defaults to the last SVG
#' @return gene set
#' @export
#' @examples
#' data("res_list_liver_bud")
#' corgi_gene_set <- select_top_corgi_genes(corgi_output_toy_example,10)
#' corgi_gene_set
select_top_corgi_genes <- function(corgi_output,n,phase_use = NULL,SVG_use = NULL){
  if(is.null(phase_use)){
    phase_use <- length(corgi_output)
  }

  res <- corgi_output[[phase_use]]

  if(is.null(SVG_use)){
    SVG_use <- ncol(res)-1
  }
  alignment_score <- res[,SVG_use]/res[,"num.times.sampled"]

  names(
    head(
      sort(alignment_score,
           decreasing = T),
      n)
  )
}
