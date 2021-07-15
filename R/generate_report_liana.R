#'
#'Run all LR Downstream analysis
#'
#'This function loads the single conditions LR outputs and return the LR network based analysis.
#'It assumes that the table present the following columns Ligand, Ligand.Cluster, Receptor,Receptor.Cluster and MeanLR/another
#'measure
#'@param lianalr liana object
#'@param genes list of genes to be considered in the sankey plots
#'@param out_path output directory path
#'@param sep character used on csv
#'@param threshold percentage of edges to be pruned
#'@param colors celltypes colorscheme
#'@param out_file output file names
#'@param report decide if a report is generated or not
#'@param output_fmt rmarkdown render output format parameter
#'@importFrom tidyr %>%
#'@return Rmarkdown report all objects from each step
#'@export
#'@examples
#'paths <- c('CTR' = system.file("extdata",
#'                               "ctr_nils_bm_human.csv",
#'                               package = "CrossTalkeR"),
#'           'EXP' = system.file("extdata",
#'                               "exp_nils_bm_human.csv",
#'                               package = "CrossTalkeR"))
#'output =  system.file("extdata", package = "CrossTalkeR")
#'genes <- c('TGFB1')
#'data <- generate_report(lrpaths = paths,
#'                        genes = genes,
#'                        out_path = paste0(output,'/'),
#'                        threshold = 0,
#'                          out_file = "report.html")
#'
generate_report_liana <- function(lianalr,
                            genes = NULL,
                            out_path,
                            threshold = 0,
                            colors = NULL,
                            out_file = NULL,
                            report = TRUE,
                            output_fmt = "html_document",
                            sel_columns = c("source",
                                             "target",
                                             "Ligand",
                                             "Receptor",
                                             'meanlr')) {
  # Creating the single condition Object
  index_single <- system.file("templates",
                              "FinalReport_Single.Rmd",
                              package = "CrossTalkeR")
  index <- system.file("templates",
                       "FinalReport.Rmd",
                       package = "CrossTalkeR")
  message("Reading Data")
  data <- read_lr_single_condiction_liana(lianalr,
                                    out_path,
                                    colors,
                                    sel_columns)
  # Generating the single condition report
  lrobj_path1 <- paste0(out_path, "LR_data_final.Rds")
  message("Calculating CCI Ranking")
  data <- suppressWarnings({ ranking(data, out_path, slot = "graphs")})
  message("Calculating GCI Ranking")
  data <- suppressWarnings({ ranking(data, out_path, slot = "graphs_ggi")})
  message("Annotating the top Cell Genes")
  data <- suppressWarnings({kegg_annotation(data=data,
                          slot='rankings',out_path=out_path)})
  message("Defining templates")
  param_single <- list(obj1 = lrobj_path1,
                       obj2 = genes,
                       thr = threshold)
 if (report) {
    message("Generating Report")
      message('Preparing Single Phenotype Report')
      rmarkdown::render(index_single,
                      output_format = output_fmt,
                      output_dir = out_path,
                      output_file = paste0("Single_",out_file),
                      intermediates_dir	= out_path,
                      params = param_single,quiet = TRUE)
  }
  message("DONE!")
  return(data)
}
