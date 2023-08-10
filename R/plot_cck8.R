#' Title plot_cck8
#'
#' @param data A data frame, to see the example please run load(system.file("data", "cck8.Rdata", package = "OneGene"))
#' @param time The x.axis to plot, it alse can be the drug concentration.
#' @param groups A vector contains all groups in your experiment.
#' @param clean Whether the data is clean.
#' @param plot_title The title to plot.
#' @param x_title The title in x.
#' @param y_title The title in y.
#' @param legend.position A vector for the significant mark position.
#' @param size The text size.
#' @param hide_ns Whether to hide the no significant, default is true.
#' @param stat_method The statistical method to us, default is wilcox.test.

#'
#' @return
#' @export
#'
#' @examples plot_cck8(data, time = "时间", stat_method ='anova',groups = c('对照组', '实验组1','实验组2'))
plot_cck8 <- function(data = NULL, time, groups,clean = F,
                      plot_title = 'CCK8 assay',
                      x_title = "time", y_title = "OD450",
                      legend.position = c(0.1, 0.5),
                      size = 16,
                      hide_ns = TRUE, stat_method = "wilcox.test") {
  if (is.null(data)) {
    path <- system.file("data", "cck8.Rdata", package = "OneGene")
    load(path)
    message('===> No data provide, it will use the example data')
    groups = c('对照组', '实验组1','实验组2')
    time = "时间"
  }


  library(ggpubr)
  library(ggsignif)
  library(tidyverse)

  if (is.null(time) || is.null(groups)) {
    stop("Please provide time and groups.")
  }

  if (clean == F) {
    plot_data <- data %>%
      select(time, groups) %>%
      pivot_longer(cols = -time, names_to = 'group', values_to = 'OD450')
  } else {
    plot_data <- data
  }



  ggline(plot_data, x = time, y = "OD450", group = 'group',size = 1.25,
         color = "group", shape = 'group', title = plot_title,
         add = "mean_sd", palette = 'lancet',
         xlab = x_title, ylab = y_title, legend = legend.position,
         ggtheme = theme_minimal(base_size = size)) +
    stat_compare_means(label = "p.signif", aes(group = group),
                       hide.ns = hide_ns, method = stat_method)
}
