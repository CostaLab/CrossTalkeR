% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot_cci}
\alias{plot_cci}
\title{Plot Cell Cell Interaction}
\usage{
plot_cci(
  graph,
  plt_name,
  emax = NULL,
  leg = FALSE,
  low = 25,
  high = 75,
  ignore_alpha = FALSE,
  log = FALSE,
  efactor = 8,
  vfactor = 12,
  vnames = TRUE,
  pg = NULL,
  vnamescol = NULL,
  colors,
  coords,
  col_pallet = NULL,
  standard_node_size = 20,
  pg_node_size_low = 10,
  pg_node_size_high = 60,
  arrow_size = 0.4,
  arrow_width = 0.8,
  node_label_position = 1.2,
  node_label_size = 0.6,
  score_filter = 0,
  cell_name_filter = NULL
)
}
\arguments{
\item{graph}{Paths of single condition LR data}

\item{plt_name}{Plot Name (Title)}

\item{emax}{Max MeanLR across the all inputs, if its not defined,
the method going to consider the max find within a sample}

\item{leg}{Set color legend}

\item{low}{Lower threshold: This parameter low and high defines the edges}

\item{high}{Higher threshould
which will be filtered. Edges within the interval [low\,high] are filtered.}

\item{ignore_alpha}{not include transparency on the plot}

\item{log}{logscale the interactions}

\item{efactor}{edge scale factor}

\item{vfactor}{edge scale factor}

\item{vnames}{remove vertex labels}

\item{pg}{pagerank values}

\item{colors}{Cell type (Cluster) Colors}

\item{coords}{object coordinates}

\item{col_pallet}{Custom color pallet for the Edges}

\item{standard_node_size}{Node size if no Pagerank values are given}

\item{pg_node_size_low}{Smallest node size if Pagerank values are given}

\item{pg_node_size_high}{Largest node size if Pagerank values are given}

\item{arrow_size}{Scale value for the arrow size}

\item{arrow_width}{Scale value for the arrow width}

\item{node_label_position}{Scale Factor to move the node labels}

\item{node_label_size}{Scale Factor to change the node label size}

\item{score_filter}{Filter Graph by LR Score symmetrically around 0}

\item{cell_name_filter}{Filter incoming and outgoing interactions by defined cell types}
}
\value{
R default plot
}
\description{
This function does a CCI plot
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

genes <- c("TGFB1")

output <- system.file("extdata", package = "CrossTalkeR")
data <- generate_report(paths,
  genes,
  out_path = paste0(output, "/"),
  threshold = 0,
  out_file = "vignettes_example.html",
  output_fmt = "html_document",
  report = FALSE
)
plot_cci(
  graph = data@graphs$CTR,
  colors = data@colors,
  plt_name = "Example 1",
  coords = data@coords[igraph::V(data@graphs$CTR)$name, ],
  emax = NULL,
  leg = FALSE,
  low = 0,
  high = 0,
  ignore_alpha = FALSE,
  log = FALSE,
  efactor = 8,
  vfactor = 12,
  vnames = TRUE
)
}
