#' Title plot_boxblot
#'
#' @param dataset
#' @param group
#' @param values
#' @param title
#' @param x.title
#' @param y.title
#' @param color
#'
#' @return
#' @export
#'
#' @examples plot_boxplot(dataset, group, values, title = "Boxplot of Sepal Length", x.title = "Species", y.title = "Sepal Length")
plot_boxplot <- function(dataset, group, values, title = "", x.title = "", y.title = "", color = c("#5CB85C", "#337AB7", "#F0AD4E", "#D9534F")) {

  color <- color  # 设置颜色

  # 使用ggplot2库绘制箱线图
  ggboxplot(dataset, x = group, y = values, color = group, palette = color) +

    # 设置背景为空白
    theme(panel.background = element_blank()) +

    # 设置坐标轴线为黑色
    theme(axis.line = element_line(colour = "black")) +

    # 设置x轴标题为空白
    theme(axis.title.x = element_blank()) +

    # 隐藏图例
    theme(legend.position = "NA") +

    # 设置x轴标题
    labs(x = x.title, y = y.title, title = title) +

    # 设置y轴标题字体大小
    theme(axis.title.y = element_text(size = 15)) +

    # 设置x轴刻度标签字体大小、角度、垂直对齐方式和水平对齐方式
    theme(axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust = 1)) +

    # 添加均值比较标记
    stat_compare_means()
}







