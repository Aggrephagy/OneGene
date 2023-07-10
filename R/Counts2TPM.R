#' Counts2TPM
#'
#' @param path
#' @param index
#'
#' @return
#' @export
#'
#' @examples
Counts2TPM <- function(expr) {
  library(IOBR)
  exp = IOBR::count2tpm(countMat = expr,idType = "SYMBOL")
  exp = as.data.frame(expr)
}
