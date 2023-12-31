% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_scRNA_CellChat.R
\name{run_scRNA_CellChat}
\alias{run_scRNA_CellChat}
\title{Title run_scRNA_CellChat}
\usage{
run_scRNA_CellChat(
  data,
  idents = "celltype",
  save_path = NULL,
  save_name = NULL,
  CellChat_DB = NULL
)
}
\arguments{
\item{data}{A seurat object}

\item{idents}{It is what group you focus on to calculate the interactions, the default is 'celltype'}

\item{save_path}{The path you want to save, the the default is './'}

\item{save_name}{The file name you want to save, the the default is 'CellChat_res.Rdata'}

\item{CellChat_DB}{'CellChatDB.human' or 'CellChatDB.mouse', the default is 'CellChatDB.human'}
}
\description{
Title run_scRNA_CellChat
}
\examples{
run_scRNA_CellChat(scRNA = scRNA_chat,idents = 'celltype',save_path = './output_data/',save_name = 'CellChat.Rdata')
}
