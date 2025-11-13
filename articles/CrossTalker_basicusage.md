# CrossTalkeR - Basic Usage

## CrossTalkeR - Basic Execution

Here we like to give a short introduction on the basic usage of
CrossTalkeR.

### Basic Execution

Here we show the most basic way to execute CrossTalkeR. All you need to
provide are the paths to the ligand-receptor interaction analysis
results and a path to save the results to:

``` r
library(CrossTalkeR)

paths <- c(
  "Condition1" = "/path/to/condition1/LR-interactions.csv",
  "Condition2" = "/path/to/condition2/LR-interactions.csv"
)

output_path <- "/path/to/output/folder/"
data <- generate_report(paths,
                out_path=output,
                out_file = 'vignettes_example.html',
                output_fmt = "html_document",
                report = TRUE,
                org = "hsa", 
                filtered_net=TRUE)
```

In case the **tables are already obtained inside your R session**, one
can use the following approach(This also apply for the functions below).

``` r
obj1 <- "a condition1 data frame object"
obj2 <- "a condition2 data frame object"

paths <- list('Condition1' = obj1,
           'Condition2' = obj2)

output <- "/path/to/output/folder/"
data <- generate_report(paths,
                out_path=output,
                out_file = 'vignettes_example.html',
                output_fmt = "html_document",
                report = TRUE,
                org = "hsa", 
                filtered_net=TRUE)
```

### Performing only the analysis

To just **perform the analysis** part of CrossTalker, **without
generating the report**:

``` r
library(CrossTalkeR)

paths <- c(
  "Condition1" = "/path/to/condition1/LR-interactions.csv",
  "Condition2" = "/path/to/condition2/LR-interactions.csv"
)

output_path <- "/path/to/output/folder/"
data <- analise_LR(paths,
                out_path=output,
                org = "hsa", 
                filtered_net=TRUE)
```

### Getting Report from an existing project

An addictional possibility is to **generate the reports from existing
CrossTalkeR objects**. In this case the out_path should point to the
folder with the ‘LR_data_final.Rds’ file.

``` r
library(CrossTalkeR)

output <- "/path/to/output/folder/"
data <- make_report(out_path=output,
                out_file = 'vignettes_example.html',
                output_fmt = "html_document",
                report = TRUE,
                org = "hsa", 
                filtered_net=TRUE)
```
