#' @title xena_clean
#'
#' @param Rdata.path
#' @param project
#' @param data_type
#'
#' @return
#' @export
#'
#' @examples xena_clean(Rdata.path = "./output_xena/counts_xena_LUAD.Rdata",project = 'LUAD',data_type = 'counts')
xena_clean <- function(Rdata.path,project,data_type = 'counts') {
  library(tidyverse)
  load(Rdata.path)
  # 加载注释文件
  # coding_gene = data.table::fread('../HCC/protein_coding.csv',
  #                                 data.table = F) %>%
  #   filter(gene_type == "protein_coding") %>%
  #   distinct(gene_name,.keep_all = T)
  # save(coding_gene,file = './data/protein_coding.Rdata')

  coding_gene <- system.file("data", "protein_coding.Rdata", package = "OneGene")
  load(coding_gene)

  coding_gene <- coding_gene%>%
    filter(gene_type == "protein_coding") %>%
    distinct(gene_name,.keep_all = T)
  # 共有基因
  com_gene = intersect(exp$Ensembl_ID,coding_gene$gene_id)


  # 过滤（只保留编码基因）
  coding_gene = coding_gene[coding_gene$gene_id%in%com_gene,]
  exp = exp[exp$Ensembl_ID%in%com_gene,]

  exp = exp %>% column_to_rownames(var = "Ensembl_ID")
  exp = exp[coding_gene$gene_id,]
  identical(row.names(exp),coding_gene$gene_id)# TRUE
  row.names(exp) = NULL
  row.names(exp) = coding_gene$gene_name

  #对表达矩阵提取编码基因#####################

  com_gene = intersect(colnames(exp),cli$submitter_id.samples)#407
  com_gene = intersect(com_gene,surv$sample)

  exp = exp[,com_gene]
  exp = 2^exp -1
  cli = cli[match(com_gene,cli$submitter_id.samples),]
  surv = surv[match(com_gene,surv$sample),]

  # 保存数据###########################
  identical(colnames(exp),surv$sample)#TRUE
  identical(colnames(exp),cli$submitter_id.samples)#TRUE
  # 所有的顺序，数量一致
  save(exp,cli,surv,file =  paste0("./output_xena/",data_type,"_xena_clean_", project,".Rdata"))
}


