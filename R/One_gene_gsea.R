#' @Title One_gene_gsea
#'
#' @param deg_list
#' @param logFC
#'
#' @return
#' @export
#'
#' @examples gsea_res = run_GSEA(deg_list = res[[1]])
One_gene_gsea <- function(deg_list,logFC = 'logFC',
                          orgDB=as.name('org.Hs.eg.db'),
                          p.Cutoff = 0.05,gene_type="SYMBOL") {
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(enrichplot)
  # library(GseaVis)
  geneList = deg_list$logFC;names(geneList)=row.names(deg_list)
  geneList = sort(geneList,decreasing = T);head(geneList)
  enrich <- gseGO(geneList,
                  ont = "BP",    # 可选"BP"、"MF"、"CC"三大类或"ALL"
                  OrgDb = orgDB,    # 使用人的OrgDb
                  keyType = gene_type,    # 基因id类型
                  pvalueCutoff = p.Cutoff,
                  pAdjustMethod = "BH"  ) # p值校正方法
}
