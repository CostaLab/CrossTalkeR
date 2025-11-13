# CrossTalkeR Cookbook

## Version Info

``` r
suppressPackageStartupMessages({
  require(CrossTalkeR)
})
suppressPackageStartupMessages({
  require(igraph)
})
suppressPackageStartupMessages({
  require(ggraph)
})
suppressPackageStartupMessages({
  require(ggplot2)
})
suppressPackageStartupMessages({
  require(EnhancedVolcano)
})
```

**R version**: R version 4.5.2 (2025-10-31)

**Package version**: 2.0.0

## Generate Report Example

In our vignette we provide examples on how to analyse cell interactions
from a human myelofibrosis single cell RNA-seq dataset.

``` r
data(CTR)
data(EXP)
paths <- list('CTR'=CTR,
              'EXP'=EXP)

genes <- c('TGFB1|L')
output <- system.file("extdata", package = "CrossTalkeR")
data <- generate_report(paths,
                        genes,
                        out_path=file.path(output),
                        threshold=0,
                        out_file = 'vignettes_example.html',
                        output_fmt = "html_document",
                        report = TRUE)
```

## Individual Visualization

### CCI

``` r
plot_cci(
  graph = data@graphs$CTR,
  colors = data@colors,
  plt_name = "Example 1",
  coords = data@coords[V(data@graphs$CTR)$name, ],
  emax = NULL,
  leg = FALSE,
  low = 0,
  high = 0,
  ignore_alpha = FALSE,
  log = FALSE,
  efactor = 8,
  vfactor = 12
)
```

### Sankey plot

``` r
plot_sankey(
  lrobj_tbl = data@tables$EXP_x_CTR,
  target = c("TGFB1|L"),
  ligand_cluster = NULL,
  receptor_cluster = NULL,
  plt_name = "TGFB1"
)
```

### Report Samples

[Single](Single_vignettes_example.md)

[Comparative](Comparative_vignettes_example.md)

## Session information
