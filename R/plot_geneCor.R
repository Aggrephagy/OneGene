#' Title plot_geneCor
#'
#' @param expr
#' @param gene1
#' @param gene2
#' @param cor_method
#'
#' @return
#' @export
#'
#' @examples
plot_geneCor <- function(expr,gene1,gene2,cor_method="pearson") {
  library(ggplot2)
  library(dplyr)
  gene1 <<- gene1;  gene2 <<- gene2; cor_method <<- cor_method
  plot_data = t(expr[c(gene1,gene2),]) %>% as.data.frame()
  corT=cor.test(unlist(c(plot_data[gene1])),unlist(c(plot_data[gene2])),method=cor_method)
  cor=corT$estimate
  pValue=corT$p.value
  ggplot(plot_data, aes_string(gene1, gene2)) +
    xlab(gene1)+ylab(gene2)+
    geom_point()+ geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
    stat_cor(method = cor_method, aes(x =unlist(c(plot_data[gene1])),
                                      y =unlist(c(plot_data[gene2]))))
}
