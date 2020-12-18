#'
#'Run all LR Downstream analysis
#'
#'This function loads the single conditions LR outputs and return the LR network based analysis.
#'It assumes that the table present the following columns Ligand, Ligand.Cluster, Receptor,Receptor.Cluster and MeanLR/another
#'measure
#'@param lrpaths Paths of single condition LR data
#'@param sep character used on csv
#'@return Rmarkdown report all objects from each step
#'@importFrom tidyr %>%
#'@export
generate_report <- function(lrpaths,
                            genes = NULL,
                            out_path,
                            sep = ",",
                            threshold = 0,
                            colors = NULL,
                            out_file = NULL,
                            report = TRUE,
                            output_fmt = "html_document") {
  # Creating the single condition Object
  index_single <- system.file("templates",
                              "FinalReport_Single.Rmd",
                              package = "CrossTalkeR")
  index <- system.file("templates",
                       "FinalReport.Rmd",
                       package = "CrossTalkeR")
  comp <- system.file("templates",
                      "Comparative_Condition.Rmd",
                      package = "CrossTalkeR")
  message("Reading Files")
  data <- read_lr_single_condiction(lrpaths,
                                    out_path,
                                    sep = ",",
                                    colors)
  # Obtaining the differential table
  message("Create a Differential Table")
  if (length(lrpaths) > 1) {
    data <- create_diff_table(data, out_path)
  }
  message("Defining templates")
  # Generating the single condition report
  lrobj_path1 <- paste0(out_path, "LR_data_final.Rds")
  data <- ranking(data, out_path, slot = "graphs")
  data <- ranking(data, out_path, slot = "graphs_ggi")

  param_single <- list(single = single,
                       obj1 = lrobj_path1,
                       obj2 = genes,
                       thr = threshold)
  param_comp <- list(single = single,
                     comp = comp,
                     obj1 = lrobj_path1,
                     obj2 = genes,
                     thr = threshold)
  if (report) {
    message("Generating Report")
    rmarkdown::render(index_single,
                    output_format = output_fmt,
                    output_dir = out_path,
                    output_file = paste0("Single_",
                    out_file),
                    params = param_single)
    rmarkdown::render(index,
                      output_format = output_fmt,
                      output_dir = out_path,
                      output_file =  paste0("Comparative_", out_file),
                      params = param_comp)

  }
  return(data)
}
