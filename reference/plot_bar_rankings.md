# This function generates the barplot for a given network ranking on the CGI level. Further, the genes can be filtered by selected gene types to filter the plot.

This function generates the barplot for a given network ranking on the
CGI level. Further, the genes can be filtered by selected gene types to
filter the plot.

## Usage

``` r
plot_bar_rankings(
  data_object,
  table_name,
  ranking,
  type = NULL,
  filter_sign = NULL,
  mode = "cci",
  top_num = 10
)
```

## Arguments

- data_object:

  LRobject with all data

- table_name:

  name of the ranking table

- ranking:

  name of the network ranking to use

- type:

  gene type (L,R,TF, LR/RL, RTF/TFR, LTF/TFL)

- filter_sign:

  show all (NULL), only positive (pos), or only negativ (neg) results

- mode:

  mode of the ranking (cci, ggi)

- top_num:

  number of top genes to show in positive or negative direction

## Value

R default plot
