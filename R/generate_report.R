#'
#'Run all LR Downstream analysis
#'
#'This function loads the single conditions LR outputs and return the LR network based analysis.
#'It assumes that the table present the following columns Ligand, Ligand.Cluster, Receptor,Receptor.Cluster and MeanLR/another
#'measure
#'@param LRpaths Paths of single condition LR data
#'@param sep character used on csv
#'@return Rmarkdown report all objects from each step
#'@importFrom tidyr %>%
#'@export
#'
generate_report <- function(LRpaths,genes=NULL,out_path,sep=',',threshold=50,colors=NULL,out_file=NULL,report=TRUE){
  # Creating the single condition Object
  single <- system.file('templates','Single_Condition.Rmd', package = 'LRAnalytics')
  comp <- system.file('templates','Comparative_Condition.Rmd', package = 'LRAnalytics')
  sankey <- system.file('templates','SankeyPlots.Rmd', package = 'LRAnalytics')
  tbl <- system.file('templates','Tables.Rmd', package = 'LRAnalytics')
  index <- system.file('templates','FinalReport.Rmd', package = 'LRAnalytics')

  message('Read Files')
  data <- read_lr_single_condiction(LRpaths,out_path,sep=',',colors)
  # Obtaining the differential table
  message('Create a Differential Table')
  if(length(LRpaths)>1){
    data <- create_diff_table(data,out_path)
  }
  message('Defining templates')
  # Generating the single condition report
  lrObj_path1 <- paste0(out_path,'LR_data_final.Rds')
  data <- ranking_cci(data,out_path)

  param <- list(single=single,
                comp = comp,
                obj1 = lrObj_path1,
                obj2 = genes,
                thr = threshold)
  if(report){
  message('Generating Report')
  rmarkdown::render(index,
                    output_format = 'html_document',
                    output_dir = out_path,output_file = out_file,params = param)
  }
  return(data)
}






