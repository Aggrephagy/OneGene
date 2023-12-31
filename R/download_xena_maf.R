#' Title download_Xena_maf
#'
#' @param project
#'
#' @return A dataframe for MuTect2 data frome xena, which can be deal with maftools.
#' @export
#'
#' @examples download_maf(project = 'LUAD')
download_Xena_maf <- function(project="eg:STAD") {
  library(tidyverse)
  library(UCSCXenaTools)
  # 创建输出文件夹
  if (!dir.exists("output_xena")) {
    dir.create("output_xena")
  }

  project_data <- UCSCXenaTools::XenaScan(pattern = project)
  index <- str_detect(project_data$Label, pattern = "MuTect2*")
  table(index)
  cohort <- project_data[index,]$XenaCohorts %>% unique()
  cohort
  DataSubtype <- project_data[index,]$DataSubtype
  DataSubtype

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
  # 下载表达矩阵
  dataset <- project_data %>%
    XenaGenerate(subset = XenaCohorts == cohort &
                   DataSubtype == DataSubtype)

  message("=> Starting download! \n")

  maf <- dataset %>%
    XenaFilter(filterDatasets = dataset@datasets[str_detect(dataset@datasets,'mutect2')]) %>%
    XenaQuery() %>%
    XenaDownload() %>%
    XenaPrepare()

  # 预处理开始
  message("=> Start preprocessing! \n")

  colnames(maf) <- c("Tumor_Sample_Barcode", "Hugo_Symbol",
                     "Chromosome", "Start_Position",
                     "End_Position", "Reference_Allele", "Tumor_Seq_Allele2",
                     "HGVSp_Short" , 'effect' ,"Consequence",
                     "vaf" )
  maf$Entrez_Gene_Id <- 1
  maf$Center <- 'ucsc'
  maf$NCBI_Build <- 'GRCh38'
  maf$Strand <- '+'
  maf$Variant_Classification <- maf$effect
  maf$Tumor_Seq_Allele1 <- maf$Reference_Allele
  maf$Variant_Type <- ifelse(
    maf$Reference_Allele %in% c('A','C','T','G') & maf$Tumor_Seq_Allele2 %in% c('A','C','T','G'),
    'SNP','INDEL'
  )
  data.table::fwrite(maf,file = paste0("./output_xena/",project,'_maf.txt'))
  # 成功
  message("=> Success! \n")
  return(maf)
}
