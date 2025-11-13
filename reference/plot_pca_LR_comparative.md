# This function is a proxy to the PCA plot in comparative conditions

This function is a proxy to the PCA plot in comparative conditions

## Usage

``` r
plot_pca_LR_comparative(
  lrobj_tblPCA,
  pca_table,
  dims = c(1, 2),
  ret = F,
  ggi = TRUE,
  include_tf = TRUE,
  gene_types = "all"
)
```

## Arguments

- lrobj_tblPCA:

  LRobject table with all data

- pca_table:

  table entry

- dims:

  PCA dims

- ret:

  return plot

- ggi:

  GGI mode

- include_tf:

  intracellular option

- gene_types:

  filter option of genes

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
