#' Title plot_deg_volcano
#'
#' @param plot_data
#' @param logFC_name
#' @param pvalue_name
#' @param p.cutoff
#' @param label_text_size
#' @param custum_label
#'
#' @return
#' @export
#'
#' @examples plot_deg_volcano(plot_data = deg_res,logFC_name = 'logFC',pvalue_name = "P.Value",label_text_size = 4,custum_label = c('Cd86'),p.cutoff = 0.05)
plot_deg_volcano <- function(plot_data,logFC_name,pvalue_name,p.cutoff=0.05,
                             label_text_size=3,custum_label=NULL) {
  library(ggplot2)
  library(latex2exp)
  library(ggrepel)
  # 全局变量：
  Pvalue <<- pvalue_name; LogFC <<- logFC_name
  # -log10(padj)小于2的设置为灰色
  plot_data <- plot_data %>%
    dplyr::mutate('color' = c(ifelse(-log10(plot_data[Pvalue]) < -log10(p.cutoff), "#bfc0c1",
                                     ifelse(plot_data[LogFC]>0, "#e88182", "#6489b2") )))
  plot_data =plot_data %>% dplyr::mutate(pvalue = c(.[Pvalue][[1]]))
  plot_data =plot_data %>% dplyr::mutate(logFC = c(.[LogFC][[1]]))
  plot_data$label <- rep("", nrow(plot_data))
  if (is.null(custum_label)) {
    plot_data$label[order(plot_data$pvalue)[1:20]] <- rownames(plot_data)[order(plot_data$pvalue)[1:20]]
  }else{
    plot_data$label[which(rownames(plot_data)==custum_label)] <- custum_label
  }


  # 绘图
  ggplot(plot_data) +
    geom_point(aes(logFC, -log10(pvalue),
                   # shape = group,
                   size = -log10(pvalue)),
               color = plot_data$color,
               alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "longdash") +
    geom_hline(yintercept = -log10(p.cutoff), linetype = "longdash") +
    geom_text_repel(aes(x = logFC, y = -log10(pvalue), label = label),
                    size = label_text_size, color = plot_data$color,
                    max.overlaps = 1000000000, key_glyph = draw_key_point) +
    scale_size_continuous(range = c(0.2, 3)) +
    theme_classic() +
    theme(legend.position = "none") +
    labs(x = TeX("$Log_2 \\textit{FC}$"),
         y = TeX("$-Log_{10} \\textit{FDR} $"))

}


