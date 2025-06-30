#'
#' Run all LR Downstream analysis
#'
#' This function loads the single conditions LR outputs and return the LR network based analysis.
#' It assumes that the table present the following columns ('source','target','gene_A','gene_B','type_gene_A','type_gene_B','MeanLR')
#' measure
#' @param lrpaths Paths of single condition LR data
#' @param genes list of genes to be considered in the sankey plots
#' @param out_path output directory path
#' @param sep character used on csv
#' @param threshold percentage of edges to be pruned
#' @param colors celltypes colorscheme
#' @param out_file output file names
#' @param report decide if a report is generated or not
#' @param output_fmt rmarkdown render output format parameter
#' @param sel_columns columns from data
#' @param org organism to be used for annotation, default is "hsa" (human)
#' @param comparison condition pairs to be used for differential analysis
#' @param filtered_net if TRUE, filter the CCI network based on p-value
#' @param score_col column name for the score used in the analysis, default is "LRScore"
#' @param p_val p-value threshold for filtering the network, default is 0.05
#' @import org.Mm.eg.db
#' @import org.Hs.eg.db
#' @importFrom tidyr %>%
#' @return Rmarkdown report all objects from each step
#' @export
#' @examples
#' paths <- c(
#'   "CTR" = system.file("extdata",
#'     "CTR_LR.csv",
#'     package = "CrossTalkeR"
#'   ),
#'   "EXP" = system.file("extdata",
#'     "EXP_LR.csv",
#'     package = "CrossTalkeR"
#'   )
#' )
#' output <- system.file("extdata", package = "CrossTalkeR")
#' genes <- c("TGFB1")
#' data <- generate_report(
#'   lrpaths = paths,
#'   genes = genes,
#'   out_path = paste0(output, "/"),
#'   threshold = 0,
#'   out_file = "report.html"
#' )
#'
generate_report <- function(lrpaths,
                            genes = NULL,
                            tf_genes = NULL,
                            out_path,
                            sep = ",",
                            threshold = 0,
                            colors = NULL,
                            out_file = NULL,
                            report = TRUE,
                            output_fmt = "html_document",
                            sel_columns = c("source", "target", "gene_A", "gene_B", "type_gene_A", "type_gene_B", "MeanLR"),
                            org = "hsa",
                            comparison = NULL,
                            filtered_net = F,
                            score_col = "LRScore",
                            p_val = 0.05) {
  data <- analise_LR(
    lrpaths,
    genes,
    tf_genes,
    out_path,
    sep,
    threshold,
    colors,
    out_file,
    output_fmt,
    sel_columns,
    org,
    comparison,
    filtered_net,
    score_col,
    p_val
  )
  if (report) {
    make_report(
      genes = genes,
      tf_genes = tf_genes,
      out_path = out_path,
      threshold = 0,
      colors = colors,
      out_file = out_file,
      output_fmt = "html_document",
      sel_columns = sel_columns
    )
  }
  message("Analysis Complete")
  return(data)
}







#' Run all LR Downstream analysis
#'
#' Core engine to generate report. Here we perform all the computation related to CrossTalkeR
#'
#'
#' @param lrpaths Paths of single condition LR data
#' @param genes list of genes to be considered in the sankey plots
#' @param out_path output directory path
#' @param sep character used on csv
#' @param threshold percentage of edges to be pruned
#' @param colors celltypes colorscheme
#' @param out_file output file names
#' @param report decide if a report is generated or not
#' @param output_fmt rmarkdown render output format parameter
#' @param sel_columns columns from data
#' @param org organism to be used for annotation, default is "hsa" (human)
#' @param comparison condition pairs to be used for differential analysis
#' @param filtered_net if TRUE, filter the CCI network based on p-value
#' @param score_col column name for the score used in the analysis, default is "LRScore"
#' @param p_val p-value threshold for filtering the network, default is 0.05
#' @import org.Mm.eg.db
#' @import org.Hs.eg.db
#' @importFrom tidyr %>%
#' @return Rmarkdown report all objects from each step
#' @export
#' @examples
#' paths <- c(
#'   "CTR" = system.file("extdata",
#'     "CTR_LR.csv",
#'     package = "CrossTalkeR"
#'   ),
#'   "EXP" = system.file("extdata",
#'     "EXP_LR.csv",
#'     package = "CrossTalkeR"
#'   )
#' )
#' output <- system.file("extdata", package = "CrossTalkeR")
#' genes <- c("TGFB1")
#' data <- generate_report(
#'   lrpaths = paths,
#'   genes = genes,
#'   out_path = paste0(output, "/"),
#'   threshold = 0,
#'   out_file = "report.html"
#' )
#'
analise_LR <- function(lrpaths,
                       genes = NULL,
                       tf_genes = NULL,
                       out_path,
                       sep = ",",
                       threshold = 0,
                       colors = NULL,
                       out_file = NULL,
                       output_fmt = "html_document",
                       sel_columns = c("source", "target", "gene_A", "gene_B", "type_gene_A", "type_gene_B", "MeanLR"),
                       org = "hsa",
                       comparison = NULL,
                       filtered_net = F,
                       score_col = "LRScore",
                       p_val = 0.05) {
  data <- read_lr_single_condition(lrpaths,
    sel_columns,
    out_path,
    sep = ",",
    colors
  )
  # Obtaining the differential table
  message("Create a Differential Table")
  if (length(lrpaths) > 1) {
    data <- create_diff_table(data, out_path, comparison, score_col)
  }
  # Generating the single condition report
  lrobj_path1 <- file.path(out_path, "LR_data_final.Rds")

  if (length(lrpaths) > 1) {
    data <- fisher_test_cci(data, "LRScore", out_path = out_path, comparison)
  }

  if (filtered_net) {
    data <- filtered_graphs(data, out_path, p_val)
  }

  message("Calculating CCI Ranking")
  data <- suppressWarnings({
    ranking(data, out_path, sel_columns = sel_columns, slot = "graphs")
  })
  message("Calculating GCI Ranking")
  data <- suppressWarnings({
    ranking(data, out_path, sel_columns = sel_columns, slot = "graphs_ggi")
  })
  message("Annotating the top Cell Genes")
  if (org == "hsa") {
    data <- suppressWarnings({
      kegg_annotation(
        data = data,
        slot = "rankings", out_path = out_path, database = org.Hs.eg.db::org.Hs.eg.db, org = org
      )
    })} else if (org == "mmu") {
    data <- suppressWarnings({
      kegg_annotation(
        data = data,
        slot = "rankings", out_path = out_path, database = org.Mm.eg.db::org.Mm.eg.db, org = org
      )
    })
  } else{
    data <- suppressWarnings({
      kegg_annotation(
        data = data,
        slot = "rankings", out_path = out_path, database = org[2], org = org[1]
      )
    })
  }

  message("Network Analysis Done")
  return(data)
}


#'
#' Generate a report for given LRObj
#'
#'
#'
#'
#' @param genes list of genes to be considered in the sankey plots
#' @param out_path output directory path
#' @param threshold percentage of edges to be pruned
#' @param colors celltypes colorscheme
#' @param out_file output file names
#' @param output_fmt rmarkdown render output format parameter
#' @param LRObj rmarkdown render output format parameter
#' @param sel_columns columns from data
#' @import org.Mm.eg.db
#' @import org.Hs.eg.db
#' @importFrom tidyr %>%
#' @return Rmarkdown report all objects from each step
#' @export
#' @examples
#' paths <- c(
#'   "CTR" = system.file("extdata",
#'     "CTR_LR.csv",
#'     package = "CrossTalkeR"
#'   ),
#'   "EXP" = system.file("extdata",
#'     "EXP_LR.csv",
#'     package = "CrossTalkeR"
#'   )
#' )
#' output <- system.file("extdata", package = "CrossTalkeR")
#' genes <- c("TGFB1")
#' data <- generate_report(
#'   lrpaths = paths,
#'   genes = genes,
#'   out_path = paste0(output, "/"),
#'   threshold = 0,
#'   out_file = "report.html"
#' )
#'
make_report <- function(genes = NULL,
                        tf_genes = NULL,
                        out_path,
                        threshold = 0,
                        colors = NULL,
                        out_file = NULL,
                        output_fmt = "html_document",
                        LRObj = NULL,
                        sel_columns = c("source", "target", "gene_A", "gene_B", "type_gene_A", "type_gene_B", "MeanLR")) {
  # Creating the single condition Object
  index_single <- system.file("templates",
    "FinalReport_Single.Rmd",
    package = "CrossTalkeR"
  )
  index <- system.file("templates",
    "FinalReport.Rmd",
    package = "CrossTalkeR"
  )
  lrobj_path1 <- paste0(out_path, "LR_data_final.Rds")
  if (!is.null(LRObj)) {
    if (!is.null(out_path)) {
      lrobj_path1 <- paste0(out_path, "LR_data_final.Rds")
      saveRDS(LRObj, paste0(out_path, "LR_data_final.Rds"))
    } else {
      message("missing output path (out_path)")
    }
  }
  if (file.exists(lrobj_path1)) {
    tmpdata <- readRDS(lrobj_path1)
    message("Defining templates")
    param_single <- list(
      obj1 = lrobj_path1,
      obj2 = genes,
      thr = threshold,
      sel = sel_columns
    )
    param_comp <- list(
      obj1 = lrobj_path1,
      obj2 = genes,
      obj3 = tf_genes,
      thr = threshold,
      sel = sel_columns
    )
    conds <- names(tmpdata@graphs)[!grepl("_x_", names(tmpdata@graphs))]
    message("Generating Report")
    if (length(conds) > 1) {
      message("Preparing Single Phenotype Report")
      rmarkdown::render(index_single,
        output_format = output_fmt,
        output_dir = out_path,
        output_file = paste0("Single_", out_file),
        intermediates_dir = out_path,
        params = param_single, quiet = TRUE
      )
      message("Preparing Comparative Phenotype Report")
      rmarkdown::render(index,
        output_format = output_fmt,
        output_dir = out_path,
        output_file = paste0("Comparative_", out_file),
        intermediates_dir = out_path,
        params = param_comp, quiet = TRUE
      )
    } else {
      message("Preparing Single Phenotype Report")
      rmarkdown::render(index_single,
        output_format = output_fmt,
        output_dir = out_path,
        output_file = paste0("Single_", out_file),
        intermediates_dir = out_path,
        params = param_single, quiet = TRUE
      )
    }
    message("Report Done!")
  } else {
    message(glue::glue("LRObj not found:{out_path}"))
  }
}
