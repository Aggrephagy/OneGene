#' Title run_scRNA_CellChat
#'
#' @param data A seurat object
#' @param idents It is what group you focus on to calculate the interactions, the default is 'celltype'
#' @param save_path The path you want to save, the the default is './'
#' @param save_name The file name you want to save, the the default is 'CellChat_res.Rdata'
#' @param CellChat_DB 'CellChatDB.human' or 'CellChatDB.mouse', the default is 'CellChatDB.human'
#' @param
#' @param
#'
#' @return
#' @export
#'
#' @examples run_scRNA_CellChat(scRNA = scRNA_chat,idents = 'celltype',save_path = './output_data/',save_name = 'CellChat.Rdata')

run_scRNA_CellChat <- function(data, idents = 'celltype',
                               save_path = NULL, save_name = NULL,
                               CellChat_DB = NULL) {

  # 检查和安装所需的包
  if (!requireNamespace("CellChat", quietly = TRUE)) {
    devtools::install_github("sqjin/CellChat", ask = FALSE, suppressUpdates = TRUE)
  }
  if (!requireNamespace("tidyverse", quietly = TRUE)) {
    install.packages("tidyverse", ask = FALSE, suppressUpdates = TRUE)
  }
  library(Seurat)
  library(CellChat)
  library(tidyverse)
  # library(future)
  # future::plan("multiprocess", workers = 12)
  # message(paste0('==> ----------',nbrOfWorkers(),' cores will be used !'))



  # 数据准备
  scRNA_chat <- data
  Idents(scRNA_chat) <- idents
  meta <- scRNA_chat@meta.data
  data_input <- as.matrix(scRNA_chat@assays$RNA@data)

  message(identical(colnames(data_input), rownames(meta)))
  message("==> Start Create CellChat Object!")

  cellchat <- createCellChat(object = data_input, meta = meta, group.by = idents)

  database <- c('CellChatDB.human', 'CellChatDB.mouse')
  if (is.null(CellChat_DB)) {
    CellChatDB <- CellChatDB.human
  } else if (!CellChat_DB %in% database) {
    message("--------\nError\n--------\t==> -------Only 'CellChatDB.human' and 'CellChatDB.mouse' supported!")
  } else if(CellChat_DB == database[1]){
    CellChatDB <- CellChatDB.human
  } else if (CellChat_DB == database[2]) {
    CellChatDB <- CellChatDB.mouse
  }

  groupSize <- as.numeric(table(cellchat@idents))
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
  cellchat@DB <- CellChatDB.use

  # 运行CellChat
  message("==> Start Run CellChat-----Please Wait--------")

  dplyr::glimpse(CellChatDB$interaction)  # 配体-受体分析
  cellchat <- subsetData(cellchat)  # 提取数据库支持的数据子集
  cellchat <- identifyOverExpressedGenes(cellchat)  # 识别过表达基因
  cellchat <- identifyOverExpressedInteractions(cellchat)  # 识别配体-受体对
  cellchat <- projectData(cellchat, PPI.human)  # 将配体、受体投射到PPI网络
  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = 10)  # 如果某些细胞群中的细胞数量很少，则过滤掉细胞间通信
  cellchat <- computeCommunProbPathway(cellchat)
  df.net <- subsetCommunication(cellchat)
  cellchat <- aggregateNet(cellchat)
  groupSize <- as.numeric(table(cellchat@idents))
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

  message("==> CellChat Done-----Start Save the data--------")

  if (is.null(save_path)) {
    path <- './'
  } else {
    path <- save_path
  }

  if (is.null(save_name)) {
    name <- 'CellChat_res.Rdata'
  } else {
    name <- save_name
  }

  # 保存文件
  # plan(sequential)#恢复单线程
  path_to_save <- file.path(path, name)
  save(cellchat, file = path_to_save)

  message("==> Success!--------")
}








