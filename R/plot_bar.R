#' Title plot_bar
#'
#' @param data A dataframe which include at least one group and the corresponding value
#' @param Group1 The first group you focuse on
#' @param Group2 The second group you focuse on
#' @param values The values of your data
#' @param compare_list The group to compare, the default is to compare all group
#' @param signif.method The statistical method to use, default is 't.test', also 'wilcox.test'
#' @param signif.position.y The significant mark position
#' @param pal The color matching plate, you can provide a color vector to change the color
#'
#' @return
#' @export
#'
#' @examples plot_bar(data,Group1  = "group1" ,Group2 = 'group2', values  =  "value")
#'
plot_bar <- function(data,Group1=NULL,Group2=NULL,values=NULL,
                     compare_list=NULL,signif.method='t.test',
                     legend_position = 'right',plot_line = F,
                     plot.title = NULL,x.title = NULL,y.title = NULL,
                     signif.position.y,pal = NULL,map_signif = T ) {
  library(ggplot2);library(tidyverse);library(ggrepel);library(ggprism)

  group_1 <- Group1
  group_2 <- Group2
  values <- values

  # 配色
  if (is.null(pal)) {
    pal <- c("#F39B7FFF","#91D1C2FF", "#8491B4FF",
             "#E64B35FF", "#4DBBD5FF", "#00A087FF")
  }


  # -------------1.数据准备-----------------
  if (is.null(group_2)) {
    df <- data[,c(group_1,values)]
    names(df) <- c('group1','value')
  } else {
    df <- data[,c(group_1,group_2,values)]
    names(df) <- c('group1','group2','value')
  }

  # -------------2.差异分析-----------------
  if (is.null(compare_list)) {
    # 识别group1中的唯一值
    compare <- as.character(unique(df$group1))

    # 创建空的compare_list列表
    compare_list <- list()

    # 循环遍历group1并创建compare_list
    for (i in 1:(length(compare)-1)) {
      for (j in (i+1):length(compare)) {
        compare_list[[length(compare_list)+1]] <- c(compare[i], compare[j])
      }
    }
    compare_list
  } else {
    compare_list <- compare_list
  }
  # -------------3.显著性标记位置-----------------
  signif.position.y <- c(max(df$value)+max(df$value)*0.2,
                         max(df$value)+max(df$value)*0.3,
                         max(df$value)+max(df$value)*0.4)

  # --------------4.绘图---------------------------
  # 重新映射数据
  plot_subdata <- df %>% group_by(group1) %>%
    summarise(across(.cols = value, mean))

  if (plot_line == F) {
    p <- ggplot(df,aes(group1,value,fill=group1,color = group1))+
      geom_bar(stat="summary",fun=mean,position="dodge")+ #绘制柱状图
      stat_summary(geom = "errorbar",fun.data = 'mean_sd', width = 0.3)+#误差棒
      labs(title = plot.title, x=x.title,y=y.title)+#标题
      theme_prism(palette = "candy_bright",
                  base_fontface = "bold", # 字体样式，可选 bold, plain, italic
                  base_family = "Arial", # 字体格式，可选 serif, sans, mono, Arial等
                  base_size = 16,  # 图形的字体大小
                  base_line_size = 0.8, # 坐标轴的粗细
                  axis_text_angle = 45)+ # 可选值有 0，45，90，270
      scale_fill_manual(values = pal)+
      scale_color_manual(values = rep('black',100))+
      theme(text = element_text(size = 14,family = 'Arial',face = 'bold'),
            legend.position = legend_position)+
      geom_signif(comparisons = compare_list,
                  map_signif_level = map_signif , #是否使用星号显示
                  test = signif.method, ##计算方法
                  y_position = signif.position.y,#图中横线位置设置
                  tip_length = c(c(0.01,0.01),
                                 c(0.01,0.01),
                                 c(0.01,0.01)),#横线下方的竖线设置
                  size=0.8,color="black")+
      geom_jitter(data=df,aes(group1,value),size=2,pch=20,color="black")
  } else {
    p <- ggplot(df, aes(x = group1, y = value, fill = group1,color = group1)) +
      # 绘制均值点
      geom_point(stat = 'summary', fun = 'mean') +
      # 绘制均值条形图
      geom_bar(stat = 'summary', fun = 'mean') +
      # 绘制抖动点图
      geom_jitter(size=2,pch=20,color="black",alpha=0.8)+
      # 绘制均值的误差线
      stat_summary(fun.data = 'mean_sd', geom = 'errorbar', width = 0.1) +
      # 绘制分组的虚线
      geom_line(data=plot_subdata, aes(x=group1,  y=value,group =1),lty=2,
                inherit.aes = F)+
      labs(title = plot.title, x=x.title,y=y.title)+#标题
      theme_prism(palette = "candy_bright",
                  base_fontface = "bold", # 字体样式，可选 bold, plain, italic
                  base_family = "Arial", # 字体格式，可选 serif, sans, mono, Arial等
                  base_size = 16,  # 图形的字体大小
                  base_line_size = 0.8, # 坐标轴的粗细
                  axis_text_angle = 45)+ # 可选值有 0，45，90，270
      scale_fill_manual(values = pal)+
      scale_color_manual(values = rep('black',100))+
      theme(text = element_text(size = 14,family = 'Arial',face = 'bold'),
            legend.position = legend_position)+
      geom_signif(comparisons = compare_list,
                  map_signif_level = map_signif , #是否使用星号显示
                  test = signif.method, ##计算方法
                  y_position = signif.position.y,#图中横线位置设置
                  tip_length = c(c(0.01,0.01),
                                 c(0.01,0.01),
                                 c(0.01,0.01)),#横线下方的竖线设置
                  size=0.8,color="black")
  }


  if (is.null(group_2)) {
    plot(p)
  } else if( plot_line == T){
    message('===> ====Error!==== -----Only One Group supported plot_line !-----')
  } else {
    p + facet_grid(~group2,scales = 'free')
  }

}
