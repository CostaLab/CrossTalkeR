
#'
#'Run and Generate all LR Downstream analysis
#'
#'This function loads the single conditions LR outputs and return the LR network
#`based analysis.
#'It assumes that the table present the following columns Ligand,
#` Ligand.Cluster,Receptor,Receptor.Cluster and MeanLR/another
#'measure
#'@slot graphs All Cell Cell Interaction Networks
#'@slot tables All tables from single condition
#'@slot max_iter  Max meanLR from all
#'@slot max_nodes All Celltype in the experiment
#'@slot coords  Cell Cell Interaction Plots
#'@slot colors  Cell type colors
#'@slot rankings Ranking of cells and Genes
#'@slot loadings  CCI values to remove multiple times genes
lrbject <- setClass("LRObj",
                    slots = list(graphs = "list",
                               graphs_ggi = "list",
                               tables = "list",
                               max_iter = "numeric",
                               max_nodes = "numeric",
                               coords = "array",
                               colors = "character",
                               rankings = "list",
                               loadings = "list"
                               )
                    )
