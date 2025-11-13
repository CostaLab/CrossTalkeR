# Generate a report for given LRObj

Generate a report for given LRObj

## Usage

``` r
make_report(
  genes = NULL,
  tf_genes = NULL,
  out_path,
  threshold = 0,
  colors = NULL,
  out_file = NULL,
  output_fmt = "html_document",
  LRObj = NULL,
  sel_columns = c("source", "target", "gene_A", "gene_B", "type_gene_A", "type_gene_B",
    "MeanLR")
)
```

## Arguments

- genes:

  list of genes to be considered in the sankey plots

- out_path:

  output directory path

- threshold:

  percentage of edges to be pruned

- colors:

  celltypes colorscheme

- out_file:

  output file names

- output_fmt:

  rmarkdown render output format parameter

- LRObj:

  rmarkdown render output format parameter

- sel_columns:

  columns from data

## Value

Rmarkdown report all objects from each step

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
data <- generate_report(
  lrpaths = paths,
  genes = genes,
  out_path = paste0(output, "/"),
  threshold = 0,
  out_file = "report.html"
)
#> Warning: file("") only supports open = "w+" and open = "w+b": using the former
#> Error in read.table(file = file, header = header, sep = sep, quote = quote,     dec = dec, fill = fill, comment.char = comment.char, ...): no lines available in input
```
