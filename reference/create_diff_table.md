# Read the lrobject and generate the comparative tables

Read the lrobject and generate the comparative tables

## Usage

``` r
create_diff_table(data, out_path, comparison = NULL, score_col = "LRScore")
```

## Arguments

- data:

  LRObj with single condition

- out_path:

  output path

- comparison:

  condition pairs to be used for differential analysis

- score_col:

  column name for the score to be used in the comparison, default is
  "LRScore"

## Value

LRObject
