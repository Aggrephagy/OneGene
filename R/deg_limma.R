#' @Title deg_limma
#'
#' @param expr
#' @param group_info
#' @param contrast1
#' @param contrast2
#' @param p
#' @param logfc
#'
#' @return
#' @export
#'
#' @examples
deg_limma <- function (expr, group_info, contrast1='tumor',
                       contrast2='normal',p = 0.05, logfc = 1) {
  library(limma)
  library(dplyr)
  compare <<- paste0(contrast1,"-", contrast2) #定义全局变量
  group <- data.frame('id' = colnames(exp.tpm),
                      'group' = group_info)
  row.names(group) <- group$id
  design <- model.matrix(~0 + factor(group$group))
  colnames(design) = levels(factor(group$group))
  contrast.matrix <- makeContrasts(compare,
                                   levels = design)
  fit <- lmFit(expr, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  deg_list <- topTable(fit2, coef = 1, n = Inf) %>% na.omit()
  deg_sig <- deg_list %>% filter(abs(logFC) > logfc & P.Value <
                                   p)
  return(list(deg_list = deg_list, deg_sig = deg_sig))
}
