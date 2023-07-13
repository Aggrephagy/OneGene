#' Title One_gene_ssgsea
#'
#' @param gene
#' @param expr
#' @param geneset
#'
#' @return
#' @export
#'
#' @examples gsea_res <- run_ssgsea(expr = expr.tpm,gene=gene,geneset = gene_set)
One_gene_ssgsea <-  function(gene,expr,geneset) {
  library(GSVA)
  library(ggpubr)
  library(ggplot2)
  gene <<- gene
  # 分组
  group = data.frame('id' = colnames(expr),
                     'group' = c(ifelse(expr[gene,]>median(t(expr[gene,])[,1]),
                                        paste0(gene,'_high'),paste0(gene,'_low')))
  )
  #进行gsva分析
  res <- gsva(as.matrix(expr), geneset, method="ssgsea",
              mx.diff=FALSE, verbose=FALSE) #注意表达谱exp载入后需转化为matrix，前面已转换
  res <- t(res)
  res <- as.data.frame(res)
  names(res) = geneset$term %>% unique()
  gene_score = cbind(group,res)
}
