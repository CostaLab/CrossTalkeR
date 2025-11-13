# Plot Cell Cell Interaction

This function does a CCI plot

## Usage

``` r
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
```

## Arguments

- graph:

  Paths of single condition LR data

- plt_name:

  Plot Name (Title)

- emax:

  Max MeanLR across the all inputs, if its not defined, the method going
  to consider the max find within a sample

- leg:

  Set color legend

- low:

  Lower threshold: This parameter low and high defines the edges

- high:

  Higher threshould which will be filtered. Edges within the interval
  \[low\\high\] are filtered.

- ignore_alpha:

  not include transparency on the plot

- log:

  logscale the interactions

- efactor:

  edge scale factor

- vfactor:

  edge scale factor

- vnames:

  remove vertex labels

- pg:

  pagerank values

- colors:

  Cell type (Cluster) Colors

- coords:

  object coordinates

- col_pallet:

  Custom color pallet for the Edges

- standard_node_size:

  Node size if no Pagerank values are given

- pg_node_size_low:

  Smallest node size if Pagerank values are given

- pg_node_size_high:

  Largest node size if Pagerank values are given

- arrow_size:

  Scale value for the arrow size

- arrow_width:

  Scale value for the arrow width

- node_label_position:

  Scale Factor to move the node labels

- node_label_size:

  Scale Factor to change the node label size

- score_filter:

  Filter Graph by LR Score symmetrically around 0

- cell_name_filter:

  Filter incoming and outgoing interactions by defined cell types

## Value

R default plot

## Examples

``` r
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
#> Warning: file("") only supports open = "w+" and open = "w+b": using the former
#> Error in read.table(file = file, header = header, sep = sep, quote = quote,     dec = dec, fill = fill, comment.char = comment.char, ...): no lines available in input
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
#> Error in data@graphs: no applicable method for `@` applied to an object of class "function"
```
