% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot_cci}
\alias{plot_cci}
\title{Plot Cell Cell Interaction}
\usage{
plot_cci(
  graph,
  colors,
  plt_name,
  coords,
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
  vnamescol = NULL
)
}
\arguments{
\item{graph}{Paths of single condition LR data}

\item{colors}{Cell type (Cluster) Colors}

\item{plt_name}{Plot Name (Title)}

\item{coords}{object coordinates}

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
}
\value{
R default plot
}
\description{
This function do a CCI plot
}
\examples{
paths <- c('CTR' = system.file("extdata",
                              "CTR_LR.csv",
                              package = "CrossTalkeR"),
          'EXP' = system.file("extdata",
                              "EXP_LR.csv",
                              package = "CrossTalkeR"))

genes <- c('TGFB1')

output_path <- system.file("extdata", package = "CrossTalkeR")
data <- generate_report(paths,
                       genes,
                       out_path = output_path,
                       threshold = 0,
                       out_file = 'vignettes_example.html',
                       output_fmt = "html_document",
                       report = FALSE)
plot_cci(graph = data@graphs$CTR,
       colors = data@colors,
       plt_name = 'Example 1',
       coords = data@coords[igraph::V(data@graphs$CTR)$name,],
       emax = NULL,
       leg = FALSE,
       low = 0,
       high = 0,
       ignore_alpha = FALSE,
       log = FALSE,
       efactor = 8,
       vfactor = 12,
       vnames = TRUE)
}
