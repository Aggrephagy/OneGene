#' Title download_xena
#'
#' @param project
#' @param data_type
#'
#' @return a Rdata file with clinical, survival, and mRNA expression.
#' @export
#'
#' @examples download_xena(project = 'LUAD',data_type = 'counts')
download_xena <- function(project, data_type=NULL){
  library(tidyverse)
  library(UCSCXenaTools)
  # -----------1.下载前检查及环境配置-----------------------
  options(timeout=10000)
  Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 10)
  check_data_type <- c('counts','fpkm')

  ## ---------data_type 检查
  if (is.null(data_type) || !data_type %in% check_data_type) {
    data_type = 'counts'
    message('===> Warning: data_type not match, it will download the counts data!')
  }

  ## ---------输出文件夹检查
  if (!dir.exists("output_xena")) {
    dir.create("output_xena")
  }

  ## ---------下载项目检查
  project_type <- c('ACC','BLCA','BRCA','CESC','CHOL','COAD',
                    'COADREAD', 'DLBC','ESCA','FPPP','GBM', 'GBMLGG',
                    'HNSC','KICH','KIPAN','KIRC',
                    'KIRP','LAML','LGG','LIHC',
                    'LUAD', 'LUSC','MESO', 'OV',
                    'PAAD','PCPG','PRAD','READ', 'SARC',
                    'SKCM','STAD', 'STES',
                    'TGCT','THCA','THYM',
                    'UCEC', 'UCS', 'UVM','PANCAN','GTEX')
  if (!project %in% project_type) {
    message(knitr::kable(data.frame(project_type)))
    message('===> Error: project not match, please check your input!')
  }

  # -----------2.数据开始下载----------------

  # -----------开始下载表达矩阵
  if (project=='PANCAN') {
    # 检索数据
    message("=> Starting query! \n")
    project_data <- XenaScan() %>% XenaGenerate()
    target_data <- "tcga_target_no_normal_RSEM_hugo_norm_count"
    message("=> Starting download! \n")
    exp <- project_data %>%
      XenaFilter(filterDatasets = target_data) %>% # 筛选数据集
      XenaQuery() %>% # 查询
      XenaDownload(method = "libcurl") %>% # 下载
      XenaPrepare() # 整理

    # -----------下载临床数据
    cli <- project_data %>%
      XenaFilter(filterDatasets = 'TCGA_TARGET_phenotype') %>% # 筛选数据集
      XenaQuery() %>% # 查询
      XenaDownload() %>% # 下载
      XenaPrepare() # 整理

    # -----------下载生存数据
    surv <- project_data %>%
      XenaFilter(filterDatasets = 'TCGA_survival_data_2.txt') %>% # 筛选数据集
      XenaQuery() %>% # 查询
      XenaDownload() %>% # 下载
      XenaPrepare() # 整理
  } else if (project=='GTEX') {
    # 检索数据
    message("=> Starting query! \n")
    project_data <- XenaScan() %>% XenaGenerate()
    target_data <- "gtex_RSEM_Hugo_norm_count"
    message("=> Starting download! \n")
    exp <- project_data %>%
      XenaFilter(filterDatasets = target_data) %>% # 筛选数据集
      XenaQuery() %>% # 查询
      XenaDownload(method = "libcurl") %>% # 下载
      XenaPrepare() # 整理

    # -----------下载临床数据
    cli <- project_data %>%
      XenaFilter(filterDatasets = 'GTEX_phenotype') %>% # 筛选数据集
      XenaQuery() %>% # 查询
      XenaDownload() %>% # 下载
      XenaPrepare() # 整理

    # -----------无生存数据
    surv <- NULL
  } else {
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
      XenaDownload(method = "libcurl") %>% # 下载
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
        XenaFilter(filterDatasets = dataset[2]) %>% # 筛选数据集
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
  }

  # 保存数据
  save(exp, cli, surv, file = paste0("./output_xena/",data_type,"_xena_", project, ".Rdata"))
  message("=>Success! 数据已保存")
}



