#' title One_gene_cor
#' @description To calculate correlation with others for one gene.
#' @param expr
#' @param gene
#' @param p.threshold
#'
#' @return a list with all correlation and significant correlation
#' @export
#'
#' @examples cor_res = One_gene_cor(expr = expr,gene = gene,p.threshold = 0.05)
One_gene_cor <- function(expr,gene,p.threshold=0.05,parallel = FALSE) {
  library(parallel)
  library(dplyr)
  gene <<- gene
  expr = t(expr) %>%
    as.data.frame()

  expr[gene] %>% head()

  # 设置并行计算的核心数
  num_cores <- detectCores()

  if (parallel) {
    corResult <- lapply(colnames(expr), function(i) {
      exprc <- cor.test(expr[,gene], expr[,i], method = "pearson")
      pvalue <- exprc$p.value
      data.frame(id = i,
                 R = exprc$estimate,
                 Low95CI = exprc$conf.int[1],
                 High95CI = exprc$conf.int[2],
                 pvalue = pvalue)
    })
  }else{
    # 使用mclapply进行多线程运行
    corResult <- mclapply(colnames(expr), function(i) {
      exprc <- cor.test(expr[,gene], expr[,i], method = "pearson")
      pvalue <- exprc$p.value
      data.frame(id = i,
                 R = exprc$estimate,
                 Low95CI = exprc$conf.int[1],
                 High95CI = exprc$conf.int[2],
                 pvalue = pvalue)
    }, mc.cores = num_cores)
  }

  # 合并多线程运行的结果
  corResult <- do.call(rbind, corResult) %>% as.data.frame()
  cor_sig = subset(corResult,pvalue<p.threshold)
  return(list(corResult,cor_sig))
}
