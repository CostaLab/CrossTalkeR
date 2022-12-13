#'Read the lrobject and generate the comparative tables
#'
#'@param data LRObj with single condition
#'@param out_path output path
#'@return LRObject
#'@import dplyr
#'@importFrom tidyr %>%
create_diff_table1 <- function(data, out_path) {
  ctr_name <- names(data@tables)[1]
  ctr_table <- data@tables[[ctr_name]]
  for (i in 2:length(data@tables)) {
    exp_name <- names(data@tables)[i]
    cmp_name <- paste0(exp_name, "_x_", ctr_name)
    exp_table <- data@tables[[exp_name]]
    tmp_data <- merge(exp_table,ctr_table,by='allpair',all = TRUE)
    tmp_data <- tmp_data %>%
                tidyr::separate(.data$allpair ,
                                c('ligpair','recpair'),
                                sep = '_',remove = F) %>%
                tidyr::separate(.data$ligpair ,
                                c('Ligand.Cluster','Ligand'),
                                sep = '/',remove = F) %>%
                tidyr::separate(.data$recpair, 
                                c('Receptor.Cluster','Receptor'),
                                sep = '/',remove = F)
    tmp_data$LRScore.x[is.na(tmp_data$LRScore.x)] <- 0
    tmp_data$LRScore.y[is.na(tmp_data$LRScore.y)] <- 0
    final_data <- tmp_data %>%
                  dplyr::mutate(LRScore = .data$LRScore.x-.data$LRScore.y) %>%
                  dplyr::mutate(cellpair = paste0(.data$Ligand.Cluster,"_",
                                                  .data$Receptor.Cluster)) %>%
                  dplyr::select(c(.data$Ligand,
                                  .data$Ligand.Cluster,
                                  .data$Receptor,
                                  .data$Receptor.Cluster,
                                  .data$LRScore,
                                  .data$cellpair,
                                  .data$ligpair,
                                  .data$recpair,
                                  .data$allpair,
                                  .data$type_gene_A.x,
                                  .data$type_gene_B.x,
                                  .data$type_gene_A.y,
                                  .data$type_gene_B.y))
    final_data <- final_data[final_data$LRScore!=0,]
    final_data <- final_data %>%
      dplyr::mutate(type_gene_A = coalesce(.data$type_gene_A.x, 
                                           .data$type_gene_A.y)) %>%
      dplyr::mutate(type_gene_B = coalesce(.data$type_gene_B.x, 
                                           .data$type_gene_B.y)) %>%
      dplyr::select(-.data$type_gene_A.x, 
                    -.data$type_gene_A.y, 
                    -.data$type_gene_B.x, 
                    -.data$type_gene_B.y)
    data@tables[[cmp_name]] <- final_data
    final <- final_data %>%
             dplyr::mutate(ccitype = paste(.data$type_gene_A,.data$type_gene_B)) %>%
             dplyr::filter(!(str_detect(.data$ccitype,"Transcription Factor"))) %>%
             dplyr::group_by(.data$cellpair) %>%
             dplyr::summarise(LRScore=sum(.data$LRScore))
    final <- final %>%
      tidyr::separate(.data$cellpair, c("u", "v"), sep="_",remove = F)
    filtervar <- grepl('Transcription',final_data[['type_gene_A']]) | grepl('Transcription',final_data[['type_gene_B']])
    raw_inter <- table(final_data$cellpair[!filtervar])
    freq <- as.array(raw_inter)[final$cellpair]- min(as.array(raw_inter)[final$cellpair])
    freq <- freq/(max(as.array(raw_inter)[final$cellpair]) -min(as.array(raw_inter)[final$cellpair]))+0.1
    final$freq <- as.array(freq)[final$cellpair]
    final$pair <- final$cellpair
    final <- dplyr::arrange(final, abs(final$LRScore))
    graph1 <- igraph::graph_from_data_frame(final[, c("u", "v", "LRScore")])
    igraph::E(graph1)$inter <- final$freq #setting thickness and weight
    igraph::E(graph1)$inter.raw <- as.array(raw_inter)[final$cellpair] #setting thickness and weight
    igraph::E(graph1)$weight <- igraph::E(graph1)$LRScore
    igraph::E(graph1)$mean <- igraph::E(graph1)$LRScore
    data@graphs[[cmp_name]] <- graph1
    graph1 <- igraph::graph_from_data_frame(final_data[, c("ligpair",
                                                           "recpair",
                                                           "LRScore")])
    igraph::E(graph1)$weight <- igraph::E(graph1)$LRScore
    igraph::E(graph1)$mean <- igraph::E(graph1)$LRScore
    data@graphs_ggi[[cmp_name]] <- graph1
  }
  saveRDS(data,file.path(out_path, "LR_data_step2.Rds"))
  return(data)
}
