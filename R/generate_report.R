#'
#'Run all LR Downstream analysis
#'
#'This function loads the single conditions LR outputs and return the LR network based analysis.
#'It assumes that the table present the following columns Ligand, Ligand.Cluster, Receptor,Receptor.Cluster and MeanLR/another
#'measure
#'@param lrpaths Paths of single condition LR data
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
#'paths <- c('CTR' = file.path(curr_path,"ctr_filtered_corrected.csv"),
#'           'EXP' = file.path(curr_path,"exp_filtered_corrected.csv"))
#'genes <- c('TGFB1','PF4','PPBP')
#'data <- generate_report(lrpaths = paths,
#'                        genes = genes,
#'                        out_path = curr_path,
#'                        threshold = 0,
#'                          out_file = "report.html")
#'
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
  message("Reading Files")
  data <- read_lr_single_condiction(lrpaths,
                                    out_path,
                                    sep = ",",
                                    colors)
  # Obtaining the differential table
  message("Create a Differential Table")
  if (length(lrpaths) > 1) {
    data <- create_diff_table_wip(data, out_path)
    data <- kegg_annotation(data=data,
                            slot='tables')

  }
  message("Defining templates")
  # Generating the single condition report
  lrobj_path1 <- paste0(out_path, "LR_data_final.Rds")
  data <- ranking(data, out_path, slot = "graphs")
  data <- ranking(data, out_path, slot = "graphs_ggi")

  param_single <- list(obj1 = lrobj_path1,
                       obj2 = genes,
                       thr = threshold)
  param_comp <- list(obj1 = lrobj_path1,
                     obj2 = genes,
                     thr = threshold)
  if (report) {
    message("Generating Report")
    if (length(lrpaths) > 1) {
        rmarkdown::render(index_single,
                      output_format = output_fmt,
                      output_dir = out_path,
                      output_file = paste0("Single_",out_file),
                      intermediates_dir	= out_path,
                      params = param_single)
        rmarkdown::render(index,
                          output_format = output_fmt,
                          output_dir = out_path,
                          output_file =  paste0("Comparative_", out_file),
                          intermediates_dir	= out_path,
                          params = param_comp)
    }else{
      rmarkdown::render(index_single,
                      output_format = output_fmt,
                      output_dir = out_path,
                      output_file = paste0("Single_",out_file),
                      intermediates_dir	= out_path,
                      params = param_single)
    }
  }
  return(data)
}
