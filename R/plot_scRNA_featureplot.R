#' Title plot_scRNA_featureplot
#'
#' @param data A seurat object
#' @param plot.gene A vector including gene to plot
#' @param pal A vector with color, also you can use the default
#' @param pal.style The color style to plot, you can provide 'A' ... 'H'
#' @param X.title The title of axsis x
#' @param y.title The title of axsis y
#' @param pt.size The point size, the default is 1
#'
#' @return
#' @export
#'
#' @examples plot_scRNA_featureplot(data = seurat_obj,plot.gene = gene )
plot_scRNA_featureplot <- function(data, plot.gene = NULL, pal = NULL,
                                   pal.style = 'D',
                                   X.title = NULL,
                                   y.title = NULL,
                                   pt.size = 1) {

  # if(! require('meme')){install.packages("meme",ask = F,suppressUpdates = T)}
  # if(! require('yyplot')){devtools::install_github("GuangchuangYu/yyplot",ask = F,suppressUpdates = T)}
  if(! require('aplot')){remotes::install_github("YuLab-SMU/aplot",ask = F,suppressUpdates = T)}
  if(! require('viridis')){install.packages("viridis",ask = F,suppressUpdates = T)}
  if(! require('Seurat')){install.packages("Seurat",ask = F,suppressUpdates = T)}
  if(! require('ggplot2')){install.packages("ggplot2",ask = F,suppressUpdates = T)}
  if(! require('tidyverse')){install.packages("tidyverse",ask = F,suppressUpdates = T)}

  if (is.null(pal)) {
    pal <- viridis(n = 3, option =  pal.style)
  } else {
    pal <- pal
  }

  # 字体
  if (length(plot.gene) > 1) {
    gene_to_plot <- plot.gene
    p <- list()
    for (i in 1:length(gene_to_plot)) {
      p[[i]] <- FeaturePlot(object = data, features = gene_to_plot[i], cols = pal,
                            order = TRUE, min.cutoff = 1, max.cutoff = 3,
                            pt.size = pt.size, blend = F) +
        theme_test() +
        guides(fill = guide_legend(override.aes = list(size = 3, alpha = 1))) +
        theme(text = element_text(family="Arial", face="italic"),
              title = element_text(size = 16),
              plot.title = element_text(hjust = 0.5),
              axis.title = element_text(size = 16),
              axis.text = element_text(size = 12),
              panel.grid = element_blank(),
              legend.title = element_blank(),
              legend.text = element_text(size = 12),
              legend.position = "right",
              legend.key.height = unit(1, "cm"),
              legend.key.width = unit(0.5, "cm"))+
        labs(x = X.title,y = y.title)
    }
    aplot::gglist(p)
  } else {
    # 绘图
    FeaturePlot(object = data, features = plot.gene, cols = pal,
                order = TRUE, min.cutoff = 1, max.cutoff = 3,
                pt.size = pt.size, blend = plot.blend) +
      theme_test() +
      guides(fill = guide_legend(override.aes = list(size = 3, alpha = 1))) +
      theme(text = element_text(family="Arial", face="italic"),
            title = element_text(size = 16),
            plot.title = element_text(hjust = 0.5),
            axis.title = element_text(size = 16),
            axis.text = element_text(size = 12),
            panel.grid = element_blank(),
            legend.title = element_blank(),
            legend.text = element_text(size = 12),
            legend.position = "right",
            legend.key.height = unit(1, "cm"),
            legend.key.width = unit(0.5, "cm"))+
      labs(x = X.title,y = y.title)
  }
}
