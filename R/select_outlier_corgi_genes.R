#' Get the gene set after running corgi
#'
#' @param corgi_output the output of run_corgi
#' @param cut_off how many standard deviations above for the corgi gene score to be kept
#' @param phase_use which phase of the corgi_output to use, defaults to the last phase
#' @param SVG_use which singular value gap (SVG) to use, defaults to the last SVG
#' @return gene set
#' @export
#' @examples
#' data("res_list_liver_bud")
#' corgi_gene_set <- select_outlier_corgi_genes(corgi_output_toy_example)
#' corgi_gene_set
select_outlier_corgi_genes <- function(corgi_output,cut_off=2,phase_use = NULL,SVG_use = NULL){
  if(is.null(phase_use)){
    phase_use <- length(corgi_output)
  }

  res <- corgi_output[[phase_use]]

  if(is.null(SVG_use)){
    SVG_use <- ncol(res)-1
  }
  res[,SVG_use]/res[,ncol(res)] -> x
  names(x)[which((abs(x - median(x)) / mad(x)) > cut_off)]
}
