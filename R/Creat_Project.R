#' @Title Creat Project
#'
#' @param ProjectName Your project name.
#'
#' @return
#' @export
#'
#' @examples creat_project(ProjectName = "ATAC")
creat_project <- function(ProjectName) {
  dir.create(ProjectName)
  dir.create(paste0(ProjectName,"/script"))
  dir.create(paste0(ProjectName,"/input_data"))
  dir.create(paste0(ProjectName,"/output_figure"))
  dir.create(paste0(ProjectName,"/output_data"))
  dir.create(paste0(ProjectName,"/other_resourse"))
  a=c("Version: 1.0",
      "RestoreWorkspace: Default",
      "SaveWorkspace: Default",
      "AlwaysSaveHistory: Default",
      "EnableCodeIndexing: Yes",
      "UseSpacesForTab: Yes",
      "NumSpacesForTab: 2",
      "Encoding: UTF-8",
      "RnwWeave: Sweave",
      "LaTeX: pdfLaTeX")
  writeLines(a,con=paste0(ProjectName,"/R.Rproj"), sep="\n")
  # file.show(paste0(ProjectName,"/code/R.Rproj")) #start project
}



