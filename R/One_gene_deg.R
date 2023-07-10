#' @title One_gene_deg
#' @description  To caculate the differtial expressed gene by limma.
#' @param expr
#' @param gene
#' @param p
#' @param logfc
#'
#' @return A list with deg_list and deg_sig.
#' @export
#'
#' @examples res = diffExprAnalysis(expr = expr, gene = gene,p = 0.01, logfc = 0.5)
One_gene_deg <- function(expr, gene, p = 0.05, logfc = 1) {
  library(limma)
  library(dplyr)
  gene <<- gene
  # 分组
  group = data.frame('id' = colnames(expr),
                     'group' = c(ifelse(expr[gene,]>median(t(expr[gene,])[,1]),
                                        paste0(gene,'_high'),paste0(gene,'_low')))
  )

  design <- model.matrix(~0+factor(group$group))
  colnames(design)=levels(factor(group$group))

  contrast.matrix<-makeContrasts(paste0(paste0(gene,'_high'),paste0('-',gene,'_low')),
                                 levels = design)

  fit <- lmFit(expr,design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)

  deg_list <-topTable(fit2, coef=1, n=Inf) %>% na.omit()  ## coef比较分组 n基因数

  deg_sig <- deg_list %>% filter(abs(logFC)>logfc&P.Value<p)

  return(list(deg_list = deg_list, deg_sig = deg_sig))
}
