#' Title plot_scRNA_volcano
#'
#' @param data A seurat object
#' @param plot.gene A vector including gene to plot in the volcano plot
#' @param plot.gene.number The default will display the top 5 up and down regulated genes
#' @param label.text.size The size of the gene, default is 80
#' @param logfc.cut The threshold of the locfc, the default is 1
#' @param p.cut The threshold of the pvalue, the default is 0.05
#' @param label.diff_pct Whether to use the pct.1 - pct.2 (diff_pct) as the label
#' @param diff_pct.cutoff The diff_pct threshold
#' @param label_postion_x_left The axis x position of the label in left
#' @param label_postion_x_right The axis x position of the label in right
#' @param label_postion_y The axis y position
#' @param base_size The text size, the default is 14
#' @param legend_size The legend size, the default is 14
#' @param legend_position The legend position, the default is 'top'
#'
#' @return
#' @export
#'
#' @examples
plot_scRNA_volcano <- function(data,plot.gene=NULL,plot.gene.number = 5,
                               label.text.size=80,logfc.cut=1,p.cut=0.05,
                               label.diff_pct = T,
                               diff_pct.cutoff = 0.4,
                               label_postion_x_left = -2.5,
                               label_postion_x_right = 2.5,
                               label_postion_y = 0.1,
                               base_size = 14,
                               legend_size = 14,
                               legend_position = 'top'
                               ) {
  library(tidyverse)
  library(ggrepel)
  options(ggrepel.max.overlaps = Inf)

  data <- data %>% mutate('gene'=row.names(.))

  if (label.diff_pct == T) {
    data <- data %>% select(gene,avg_log2FC, p_val, pct.1, pct.2) %>%
      mutate(diff_pct = pct.1 - pct.2) %>%
      mutate(group = if_else(p_val < p.cut & avg_log2FC > logfc.cut, "sig_up",
                             if_else(p_val < p.cut & avg_log2FC < -(logfc.cut), "sig_down", "no_sig"))) %>%
      mutate('label' = NA)
  } else {
    data <- data %>% select(gene,avg_log2FC, p_val, pct.1, pct.2) %>%
      mutate(diff_pct = pct.1 - pct.2) %>%
      mutate(group = if_else(p_val < p.cut & avg_log2FC > logfc.cut, "sig_up",
                             if_else(p_val < p.cut & avg_log2FC < -(logfc.cut), "sig_down", "no_sig"))) %>%
      mutate('label' = NA)
  }



  if (!is.null(plot.gene)) {
    # 标记基因
    gene_label <- as.vector(plot.gene)

    # 索引基因位置
    index <- match(gene_label,data$gene)

    # 增加基因名称为标记准备
    data$label[index] <- gene_label
  } else if(label.diff_pct == T) {
    # 前五个上调
    up_gene <- data %>% arrange(desc(diff_pct)) %>% .[,'gene'] %>% c() %>% .[1:plot.gene.number]

    # 前五个下调
    down_gene <- data %>% arrange(diff_pct) %>% .[,'gene'] %>% c() %>% .[1:plot.gene.number]

    # 标记基因
    gene_label <- c(up_gene,down_gene)

    # 索引基因位置
    index <- match(gene_label,data$gene)

    # 增加基因名称为标记准备
    data$label[index] <- gene_label

  } else {
    # 前五个上调
    up_gene <- data %>% arrange(desc(avg_log2FC)) %>% .[,'gene'] %>% c() %>% .[1:plot.gene.number]

    # 前五个下调
    down_gene <- data %>% arrange(avg_log2FC) %>% .[,'gene'] %>% c() %>% .[1:plot.gene.number]

    # 标记基因
    gene_label <- c(up_gene,down_gene)

    # 索引基因位置
    index <- match(gene_label,data$gene)

    # 增加基因名称为标记准备
    data$label[index] <- gene_label
  }


  # ------------------绘图准备--------------------
  colnames(data) <- c("gene",'log2FoldChange',
                      'pvalue',
                      "pct.1","pct.2",
                      "diff_pct", 'Group',"label"
  )

  data$Group <- factor(data$Group,levels = c("sig_down","no_sig","sig_up"))
  data$Group
  # ------------------绘图-------------------------
  ggplot(data,aes(log2FoldChange,-log10(pvalue),fill=Group))+
    geom_point(shape=21,alpha=0.5,aes(size=-log10(pvalue),color=Group))+
    scale_fill_manual(values=c("#6697ea", "grey60","#b02428"))+ #差异分组填充颜色
    scale_color_manual(values=c("#6697ea", "grey60","#b02428"))+ #边界颜色 # "midnightblue", "grey60","darkred"
    guides(size=F)+  #去掉size大小图例
    geom_label_repel(aes(x = log2FoldChange,y = -log10(pvalue),size=label.text.size,label =label),fill = "white")+
    ###坐标轴名称
    xlab(expression('log'[2]*'(FoldChange)'))+
    ylab(expression('-log'[10]*'(pvalue)'))+
    ###画出分界线
    geom_vline(xintercept=c(-logfc.cut,logfc.cut),lty=2,col="black",lwd=0.5) +
    geom_hline(yintercept = -log10(p.cut),lty=2,col="black",lwd=0.5) +
    ###指定分类颜色
    theme_test()+
    ####设置标题位置及字体大小
    theme(text = element_text(size = base_size),
          legend.text = element_text(size = legend_size),
          legend.position = legend_position,
          legend.title = element_blank())+
    ###添加文本注释
    annotate('text',x= label_postion_x_left,y=label_postion_y*max(-log10(data$pvalue)),fontface="bold.italic",size=4,color='red',
             label='Genes \n Down Regulated')+
    annotate('text',x= label_postion_x_right,y=label_postion_y*max(-log10(data$pvalue)),fontface="bold.italic",size=4,color='red',
             label='Genes \n  Up Regulated')
}

