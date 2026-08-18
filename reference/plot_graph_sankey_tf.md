# This function creates a sankey plot for a selected gene including transcription factor interactions.

This function creates a sankey plot for a selected gene including
transcription factor interactions.

## Usage

``` r
plot_graph_sankey_tf(
  lrobj_tbl,
  pagerank_table,
  target = NULL,
  cluster = NULL,
  target_type = NULL,
  plt_name = NULL,
  threshold = 50
)
```

## Arguments

- lrobj_tbl:

  LRobject table with all data

- pagerank_table:

  table with pagerank scores

- target:

  gene

- cluster:

  cluster

- target_type:

  type of target

- plt_name:

  plot title

- threshold:

  top_n n value

## Value

A ggraph/ggplot object; print it to display, or pass it to
ggplot2::ggsave() to save. In the not-found / no-target cases a
diagnostic message is printed and the message string is returned
invisibly.

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
#> Warning: file("") only supports open = "w+" and open = "w+b": using the former
#> Error in read.table(file = file, header = header, sep = sep, quote = quote,     dec = dec, fill = fill, comment.char = comment.char, ...): no lines available in input
```
