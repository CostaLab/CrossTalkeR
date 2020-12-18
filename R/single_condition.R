#'Read single condition tables
#'
#'This function loads the single conditions LR outputs and use it to generate
#' the report data and it`s object
#'It assumes that the table presents the following columns Ligand, Ligand.Cluster, Receptor,Receptor.Cluster and MeanLR/another
#'measure
#'
#'
#'@param LRpaths Named vector with the paths of each output
#'@param sep character used to divide the columns on input file
#'@param out_path Path to deposit the results
#'@return LRObject
#'@importFrom tidyr %>%
#'@export
read_lr_single_condiction <- function(paths,
                                      out_path = "/tmp/",
                                      sep = ",",
                                      colors = NULL,
                                      measure = "MeanLR") {
  data <- list()
  graphs <- list()
  conds <- names(paths)
  max <- 0
  max_nodes <- 0
  unif_celltypes <- c()
  sel_columns <- c("Ligand.Cluster",
                   "Receptor.Cluster",
                   "Ligand",
                   "Receptor",
                   measure)
  for (i in seq_len(length(paths))) {
    print(paths[i])
    data1 <- read.csv(paths[i], sep = sep) # Reading csv
    data1 <- data1[, sel_columns]
    data1$cellpair <- paste(data1$Ligand.Cluster,
                            data1$Receptor.Cluster,
                            sep = "_")
    data1$ligpair <- paste(data1$Ligand,
                           data1$Ligand.Cluster,
                           sep = "/")
    data1$recpair <- paste(data1$Receptor,
                           data1$Receptor.Cluster,
                          sep = "/")
    data1$allpair <- paste(data1$ligpair,
                           data1$recpair,
                           sep = "_")
    unif_celltypes <- unique(c(data1$Ligand.Cluster,
                               data1$Receptor.Cluster,
                               unif_celltypes)
                             )
    data1$MeanLR <- data1[,measure]
    data1 <- tibble::as_tibble(data1)
    final <- data1 %>%
      dplyr::group_by(.data$cellpair) %>%
      dplyr::summarise(MeanLR = sum(.data[[measure]]))
    aux <- final$cellpair
    clusters_num <- unique(c(unique(data1$Ligand.Cluster),
                             unique(data1$Receptor.Cluster)
                            )
                          )
    print(length(clusters_num))
    final <- final %>%
      tidyr::separate(.data$cellpair, c("u", "v"), "_")
    final$pair <- aux
    freq <- table(data1$cellpair) / max(table(data1$cellpair))
    final$freq <- as.array(freq)[final$pair]
    final <- dplyr::arrange(final, abs(final$MeanLR))
    graph1 <- igraph::graph_from_data_frame(final[, c("u", "v", measure)])
    igraph::E(graph1)$inter <- final$freq#setting thickness and weight
    igraph::E(graph1)$weight <- igraph::E(graph1)$MeanLR
    igraph::E(graph1)$mean <- igraph::E(graph1)$MeanLR
    data[[conds[i]]] <- data1
    graphs[[conds[i]]] <- graph1
    if (max(igraph::E(graph1)$mean) > max) {
      max <- max(igraph::E(graph1)$mean)
    }
    if (length(igraph::V(graph1)) > max_nodes) {
      max_nodes <- length(igraph::V(graph1))
    }
  }
  template <- igraph::make_full_graph(n = length(unif_celltypes),
                                      directed = T,
                                      loops = T)
  c <- igraph::layout.circle(template)
  if (is.null(colors)) {
    colors <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))
    colors <- colors(length(unif_celltypes))
    names(colors) <- sort(unif_celltypes)
  }
  for (g in names(graphs)) {
    sel <- match(unif_celltypes,
                 unique(igraph::V(graphs[[g]])$name),
                 nomatch = F)
    sel <- sel == 0
    if (sum(sel) != 0) {
      nodes <- seq_len(length(unif_celltypes[sel]))
      names(nodes) <- unif_celltypes[sel]
      graphs[[g]] <- igraph::add.vertices(graphs[[g]],
                                          length(nodes),
                                          attr = list(name = names(nodes)))
    }
  }
  rownames(c) <- sort(unif_celltypes)
  lr <- new("LRObj",
            graphs = graphs,
            tables = data,
            max_iter = max,
            max_nodes = max_nodes,
            coords = c,
            colors = colors,
            rankings = list())
  saveRDS(lr, paste0(out_path, "/LR_data.Rds"))
  return(lr)
}
