#' rSirem
#' 
#'  MSI data deconvolution
#' 
#' @docType package
#' @author Esteban del Castillo
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @useDynLib rSirem
#' @name rSirem
NULL  

.onUnload <- function (libpath)
{
  library.dynam.unload("rSirem", libpath)
}