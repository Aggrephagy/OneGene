#' Title download_xena
#'
#' @param project
#' @param data_type
#'
#' @return a Rdata file with clinical, survival, and mRNA expression.
#' @export
#'
#' @examples download_xena(project = 'LUAD',data_type = 'counts')
download_xena <- function(project,data_type = 'counts') {
  # 加载R包
  library(tidyverse)
  library(UCSCXenaTools)
  # 创建输出文件夹
  if (!dir.exists("output_xena")) {
    dir.create("output_xena")
  }

  # Project种类
  message(knitr::kable(data.frame(cancer = c('ACC','BLCA','BRCA','CESC','CHOL','COAD',
                                             'COADREAD', 'DLBC','ESCA','FPPP','GBM', 'GBMLGG',
                                             'HNSC','KICH','KIPAN','KIRC',
                                             'KIRP','LAML','LGG','LIHC',
                                             'LUAD', 'LUSC','MESO', 'OV',
                                             'PAAD','PCPG','PRAD','READ', 'SARC',
                                             'SKCM','STAD', 'STES',
                                             'TGCT','THCA','THYM',
                                             'UCEC', 'UCS', 'UVM'))))
  if (!data_type %in% c('counts', 'fpkm')) {
    message(knitr::kable(c('counts','fpkm')))
  }
  else {
    # 检索数据
    message("=> Starting query! \n")
    project_data <- XenaScan(pattern = project)
    # 只要GDC
    cohort <- str_extract(project_data$XenaCohorts,
                          pattern = "^GDC.*") %>%
      na.omit() %>% unique()

    # 下载表达矩阵
    dataset <- project_data %>%
      XenaGenerate(subset = XenaCohorts == cohort &
                     DataSubtype == "gene expression RNAseq")
    message("=> Starting download! \n")
    exp <- dataset %>%
      XenaFilter(filterDatasets = dataset@datasets[str_detect(dataset@datasets,data_type)]) %>% # 筛选数据集
      XenaQuery() %>% # 查询
      XenaDownload() %>% # 下载
      XenaPrepare() # 整理

    # 下载临床数据
    # 检索
    sur_pheno <- project_data %>% XenaGenerate(subset = XenaCohorts == cohort & DataSubtype == "phenotype")

    dataset2 <- sur_pheno@datasets # 索引

    if (project == "BRCA") {
      # BRCA数据集
      # 下载生存数据
      surv <- sur_pheno %>%
        XenaFilter(filterDatasets = dataset2[1]) %>% # 筛选数据集
        XenaQuery() %>% # 查询
        XenaDownload() %>% # 下载
        XenaPrepare() # 整理

      cli <- sur_pheno %>%
        XenaFilter(filterDatasets = dataset2[2]) %>% # 筛选数据集
        XenaQuery() %>% # 查询
        XenaDownload() %>% # 下载
        XenaPrepare() # 整理
    } else {
      # 其他数据集
      cli <- sur_pheno %>%
        XenaFilter(filterDatasets = dataset2[1]) %>% # 筛选数据集
        XenaQuery() %>% # 查询
        XenaDownload() %>% # 下载
        XenaPrepare() # 整理

      # 下载生存数据
      surv <- sur_pheno %>%
        XenaFilter(filterDatasets = dataset2[2]) %>% # 筛选数据集
        XenaQuery() %>% # 查询
        XenaDownload() %>% # 下载
        XenaPrepare() # 整理
    }

    # 保存数据
    save(exp, cli, surv, file = paste0("./output_xena/",data_type,"_xena_", project, ".Rdata"))
    message("=>Success! 数据已保存")
  }
}



