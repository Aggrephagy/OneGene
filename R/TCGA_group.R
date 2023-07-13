#' Title TCGA_group
#'
#' @param expr
#'
#' @return
#' @export
#'
#' @examples
TCGA_group <- function(expr) {
  group=sapply(strsplit(colnames(expr),"\\-"),"[",4)
  group=sapply(strsplit(group,""),"[",1)
  group_list=ifelse(group=="0",'tumor','normal')
  group_list=factor(group_list,levels = c('normal','tumor'))
}
