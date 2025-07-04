% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot_sankey}
\alias{plot_sankey}
\title{This function plots a sankey plot for selected genes or cell types}
\usage{
plot_sankey(
  lrobj_tbl,
  target = NULL,
  ligand_cluster = NULL,
  receptor_cluster = NULL,
  plt_name = NULL,
  threshold = 50,
  tfflag = TRUE,
  score_col = "LRScore",
  fil_col = "LRScore"
)
}
\arguments{
\item{lrobj_tbl}{LRobject table with all data}

\item{target}{gene}

\item{ligand_cluster}{Ligand Clusters}

\item{receptor_cluster}{Receptor Clusters}

\item{plt_name}{plot title}

\item{threshold}{top_n n value}

\item{tfflag}{if true, input includes transcription factors}

\item{score_col}{column name for the LR Score}

\item{fil_col}{score column name to filter}
}
\value{
R default plot
}
\description{
This function plots a sankey plot for selected genes or cell types
}
\examples{
paths <- c(
  "CTR" = system.file("extdata",
    "CTR_LR.csv",
    package = "CrossTalkeR"
  ),
  "EXP" = system.file("extdata",
    "EXP_LR.csv",
    package = "CrossTalkeR"
  )
)
output <- system.file("extdata", package = "CrossTalkeR")
genes <- c("TGFB1")

data <- generate_report(paths,
  genes,
  out_path = paste0(output, "/"),
  threshold = 0,
  out_file = "vignettes_example.html",
  output_fmt = "html_document",
  report = FALSE
)
}
