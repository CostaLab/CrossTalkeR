#'Read single condition tables
#'
#'This function loads the single conditions LR outputs and use it to generate
#' the report data and it`s object
#'It assumes that the table presents the following columns Ligand, Ligand.Cluster, Receptor,Receptor.Cluster and MeanLR/another
#'measure
#'
#'
#'@param lrpaths Named vector with the lrpaths of each output
#'@param sel selected columns
#'@param out_path Path to deposit the results
#'@param sep character used to divide the columns on input file
#'@param colors colorlist
#'@importFrom tidyr %>%
#'@import tibble dplyr assertthat stringr
#'@return LRObject
read_lr_single_condition <- function(lrpaths,
                                     sel_columns,
                                     out_path = "/tmp/",
                                     sep = ",",
                                     colors = NULL) {
  assertthat::assert_that(assertthat::not_empty(lrpaths))
  assertthat::assert_that(!is.null(names(lrpaths)))

  data <- graphs <-graphs_ggi <- list()
  conds <- names(lrpaths)
  max <- max_nodes <-0
  unif_celltypes <- c()
  for (i in seq_len(length(lrpaths))){
      ## Reading from tables
      if(is.character(lrpaths[[conds[i]]])){
        data1 <- tibble::tibble(utils::read.csv(lrpaths[[conds[i]]], 
                                                sep = sep)) # Reading csv
      }else if(is.data.frame(lrpaths[[conds[i]]]) | 
               tibble::is_tibble(lrpaths[[conds[i]]])){ ## Reading From DF
        data1 <- tibble(lrpaths[[conds[i]]])
      }
      if(!(sum(str_detect(data1$gene_A, '\\|')) > 0)){
         data1 <- add_node_type(data1)
      }
      message('NO Factor')
      data1 <- data1 %>%
          dplyr::mutate_if(is.character,as.factor) %>%
          tidyr::unite("cellpair", 
                       c(sel_columns[1],sel_columns[2]), 
                       remove = FALSE,sep='_')
      if('Ligand' %in% data1[[sel_columns[5]]]){
        data1 <- data1 %>%
                 tidyr::unite("ligpair",
                              c(sel_columns[1],sel_columns[3]),
                              remove = FALSE,
                              sep='/')%>%
                 tidyr::unite("recpair",
                              c(sel_columns[2],sel_columns[4]),
                              remove = FALSE,
                              sep='/')  %>%
                 dplyr::mutate(allpair=paste(.data$ligpair,.data$recpair,sep='_'))%>%
                 tidyr::separate(.data$allpair,
                                 c('ligpair','recpair'),
                                 sep = '_',remove = F)
      }else{
        tmp <- data1[[sel_columns[5]]]
        data1[[sel_columns[5]]] <- data1[[sel_columns[6]]]
        data1[[sel_columns[6]]] <- tmp
        tmp <- data1[[sel_columns[4]]]
        data1[[sel_columns[4]]] <- data1[[sel_columns[3]]]
        data1[[sel_columns[3]]] <- tmp
        data1 <- data1 %>%
                 tidyr::unite("ligpair",
                              c(sel_columns[1],sel_columns[3]),
                              remove = FALSE,
                              sep='/')%>%
                 tidyr::unite("recpair",
                              c(sel_columns[2],sel_columns[4]),
                              remove = FALSE,
                              sep='/')  %>%
                 dplyr::mutate(allpair=paste(.data$ligpair,.data$recpair,sep='_')) %>%
                 tidyr::separate(.data$allpair ,
                                 c('ligpair','recpair'),
                                 sep = '_',remove = F)


      }
      unif_celltypes <- unique(c(data1[[sel_columns[1]]],
                                 data1[[sel_columns[2]]],
                                 unif_celltypes)
                               )
      data1 <-data1 %>% dplyr::mutate(LRScore=data1[[sel_columns[length(sel_columns)]]])

      final <- data1 %>%  
              dplyr::mutate(ccitype = paste(data1[[sel_columns[5]]],data1[[sel_columns[6]]])) %>%
              dplyr::filter(!(str_detect(.data$ccitype,"Transcription Factor"))) %>%
              dplyr::group_by(.data$cellpair) %>%
              dplyr::summarise(LRScore=sum(.data$LRScore))
      aux <- final$cellpair
      clusters_num <- unique(c(unique(data1[[sel_columns[1]]]),
                               unique(data1[[sel_columns[2]]])))
      final <- final %>%
        tidyr::separate(.data$cellpair, c("u", "v"), "_")
      final$pair <- aux
      filtervar <- grepl('Transcription',data1[[sel_columns[5]]])| grepl('Transcription',data1[[sel_columns[6]]])
      raw_inter <- table(data1$cellpair[!filtervar])
      freq <- (as.array(raw_inter)[final$pair]- min(as.array(raw_inter)[final$pair]))
      freq <- freq/(max(as.array(raw_inter)[final$pair]) -min(as.array(raw_inter)[final$pair]))+0.1
      final$freq <- as.array(freq)[final$pair]
      final <- dplyr::arrange(final, abs(final$LRScore))
      graph1 <- igraph::graph_from_data_frame(final[, c("u", "v", "LRScore")])
      igraph::E(graph1)$inter <- final$freq#setting thickness and weight
      igraph::E(graph1)$weight <- igraph::E(graph1)$LRScore
      data[[conds[i]]] <- data1
      graphs[[conds[i]]] <- graph1
      graph2 <- igraph::graph_from_data_frame(data1[, c("ligpair",
                                                             "recpair",
                                                             "LRScore")],directed=TRUE)
      igraph::E(graph2)$weight <- igraph::E(graph2)$LRScore
      igraph::E(graph2)$mean <- igraph::E(graph2)$LRScore
      graphs_ggi[[conds[i]]] <- graph2
      if (max(igraph::E(graph1)$weight) > max) {
        max <- max(igraph::E(graph1)$weight)
      }
      if (length(igraph::V(graph1)) > max_nodes) {
        max_nodes <- length(igraph::V(graph1))
      }
    }
    template <- igraph::make_full_graph(n = length(unif_celltypes),
                                        directed = TRUE,
                                        loops = TRUE)
    c <- igraph::layout.circle(template)
    if (is.null(colors)) {
      colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))
      colors <- colors(length(unif_celltypes))
      names(colors) <- sort(unif_celltypes)
    }
    for (g in names(graphs)) {
      sel <- match(unif_celltypes,
                   unique(igraph::V(graphs[[g]])$name),
                   nomatch = FALSE)
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
              graphs_ggi = graphs_ggi,
              tables = data,
              max_iter = max,
              max_nodes = max_nodes,
              coords = c,
              colors = colors,
              rankings = list(),
              pca = list(),
              stats=list())
    saveRDS(lr,file.path(out_path, "LR_data.Rds"))
    return(lr)
}

