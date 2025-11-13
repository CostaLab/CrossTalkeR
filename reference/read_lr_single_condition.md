# Read single condition tables

This function loads the single conditions LR outputs and use it to
generate the report data and it\`s object It assumes that the table
presents the following columns Ligand, Ligand.Cluster,
Receptor,Receptor.Cluster and MeanLR/another measure

## Usage

``` r
read_lr_single_condition(
  lrpaths,
  sel_columns,
  out_path = "/tmp/",
  sep = ",",
  colors = NULL
)
```

## Arguments

- lrpaths:

  Named vector with the lrpaths of each output

- sel_columns:

  selected columns

- out_path:

  Path to deposit the results

- sep:

  character used to divide the columns on input file

- colors:

  colorlist

## Value

LRObject
