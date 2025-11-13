# This function is a proxy to the PCA plot

This function is a proxy to the PCA plot

## Usage

``` r
plot_pca(lrobj_tblPCA, curr, dims = c(1, 2), ret = F, ggi = TRUE)
```

## Arguments

- curr:

  table entry

- dims:

  PCA dims

- ret:

  return plot

- ggi:

  GGI mode

- lrobj_tbl:

  LRobject table with all data

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
