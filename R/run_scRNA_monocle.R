#' Title run_scRNA_monocle
#'
#' @param scRNA A Seurat Object
#' @param save_path A path to save your results
#'
#' @return
#' @export
#'
#' @examples
#'
run_scRNA_monocle <- function(scRNA, save_path = "./output_data/monocle2.Rdata") {
  if (!dir.exists('./output_data')) {
    dir.create('./output_data')
  }
  library(monocle)
  library(tidyverse)

  # 将scRNA对象转换为矩阵
  data <- as.matrix(scRNA@assays$RNA@counts)

  # 创建phenoData的AnnotatedDataFrame
  pd <- new('AnnotatedDataFrame', data = scRNA@meta.data)

  # 创建featureData的AnnotatedDataFrame
  fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  fd <- new('AnnotatedDataFrame', data = fData)

  # 创建CellDataSet对象
  mycds <- newCellDataSet(data,
                          phenoData = pd,
                          featureData = fd,
                          expressionFamily = negbinomial.size())

  # 估计大小因子和离散度
  mycds <- estimateSizeFactors(mycds)
  mycds <- estimateDispersions(mycds, cores = 14, relative_expr = TRUE)

  # 获取高变量基因
  disp_table <- dispersionTable(mycds)
  disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id

  # 设置排序过滤器
  mycds <- setOrderingFilter(mycds, disp.genes)

  # 绘制排序基因
  plot_ordering_genes(mycds)

  # 降维和排序细胞
  mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
  mycds <- orderCells(mycds)

  # 保存mycds对象
  save(mycds, file = save_path)

  # 返回结果
  result <- mycds
  return(result)
}

