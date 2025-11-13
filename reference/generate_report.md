# Run all LR Downstream analysis

This function loads the single conditions LR outputs and return the LR
network based analysis. It assumes that the table present the following
columns
('source','target','gene_A','gene_B','type_gene_A','type_gene_B','MeanLR')
measure

## Usage

``` r
generate_report(
  lrpaths,
  genes = NULL,
  tf_genes = NULL,
  out_path,
  sep = ",",
  threshold = 0,
  colors = NULL,
  out_file = NULL,
  report = TRUE,
  output_fmt = "html_document",
  sel_columns = c("source", "target", "gene_A", "gene_B", "type_gene_A", "type_gene_B",
    "MeanLR"),
  org = "hsa",
  comparison = NULL,
  filtered_net = F,
  score_col = "LRScore",
  p_val = 0.05
)
```

## Arguments

- lrpaths:

  Paths of single condition LR data

- genes:

  list of genes to be considered in the sankey plots

- out_path:

  output directory path

- sep:

  character used on csv

- threshold:

  percentage of edges to be pruned

- colors:

  celltypes colorscheme

- out_file:

  output file names

- report:

  decide if a report is generated or not

- output_fmt:

  rmarkdown render output format parameter

- sel_columns:

  columns from data

- org:

  organism to be used for annotation, default is "hsa" (human)

- comparison:

  condition pairs to be used for differential analysis

- filtered_net:

  if TRUE, filter the CCI network based on p-value

- score_col:

  column name for the score used in the analysis, default is "LRScore"

- p_val:

  p-value threshold for filtering the network, default is 0.05

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
