# This function plots a sankey plot for selected genes or cell types

This function plots a sankey plot for selected genes or cell types

## Usage

``` r
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
```

## Arguments

- lrobj_tbl:

  LRobject table with all data

- target:

  gene

- ligand_cluster:

  Ligand Clusters

- receptor_cluster:

  Receptor Clusters

- plt_name:

  plot title

- threshold:

  top_n n value

- tfflag:

  if true, input includes transcription factors

- score_col:

  column name for the LR Score

- fil_col:

  score column name to filter

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
