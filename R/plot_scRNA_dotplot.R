#' Title plot_scRNA_dotplot
#'
#' @param scRNA A seurat object
#' @param features A gene vector to plot
#' @param legend.position right left
#' @param plot.title
#' @param title.x
#' @param title.y
#'
#' @return
#' @export
#'
#' @examples plot_scRNA_dotplot(scRNA = scRNA_fibro,features = cell_marker)
plot_scRNA_dotplot <- function(scRNA,features,legend.position = "right",
                               plot.title="",title.x="",title.y="") {
  library(Seurat)
  library(ggplot2)
  DotPlot(object = scRNA_fibro, features = cell_marker, assay = "RNA",
          cols = c("#ffffff", "#448444")) +
    guides(color = guide_colorbar(title = 'Scaled Average Expression')) +
    coord_flip() +
    theme(legend.position = legend.position,
          legend.title  = element_text (size = 10),
          panel.border = element_rect(color = "black", fill = NA, size = 1.5),
          # 显示主网格线
          panel.grid.major = element_line(size = 0.2),
          axis.title = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.1)) +
    scale_x_discrete(position = "top") +
    scale_y_discrete(position = "right")+
    labs(title=plot.title,x=title.x,y=title.y)
}
