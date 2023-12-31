% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_cck8.R
\name{plot_cck8}
\alias{plot_cck8}
\title{Title plot_cck8}
\usage{
plot_cck8(
  data = NULL,
  time,
  groups,
  clean = F,
  plot_title = "CCK8 assay",
  x_title = "time",
  y_title = "OD450",
  legend.position = c(0.1, 0.5),
  size = 16,
  hide_ns = TRUE,
  stat_method = "wilcox.test"
)
}
\arguments{
\item{data}{A data frame, to see the example please run load(system.file("data", "cck8.Rdata", package = "OneGene"))}

\item{time}{The x.axis to plot, it alse can be the drug concentration.}

\item{groups}{A vector contains all groups in your experiment.}

\item{clean}{Whether the data is clean.}

\item{plot_title}{The title to plot.}

\item{x_title}{The title in x.}

\item{y_title}{The title in y.}

\item{legend.position}{A vector for the significant mark position.}

\item{size}{The text size.}

\item{hide_ns}{Whether to hide the no significant, default is true.}

\item{stat_method}{The statistical method to us, default is wilcox.test.}
}
\description{
Title plot_cck8
}
\examples{
plot_cck8(data, time = "时间", stat_method ='anova',groups = c('对照组', '实验组1','实验组2'))
}
