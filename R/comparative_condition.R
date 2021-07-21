#'Read the lrobject and generate the comparative tables
#'
#'@param data LRObj with single condition
#'@param out_path output path
#'@return LRObject
#'@importFrom tidyr %>%
create_diff_table <- function(data, out_path) {
  ctr_name <- names(data@tables)[1]
  ctr_table <- data@tables[[ctr_name]]
  for (i in 2:length(data@tables)) {
    exp_name <- names(data@tables)[i]
    cmp_name <- paste0(exp_name, "_x_", ctr_name)
    exp_table <- data@tables[[exp_name]]
    tmp_data <- merge(exp_table,ctr_table,by='allpair',all = TRUE)
    tmp_data <- tmp_data %>%
                tidyr::separate(allpair ,c('ligpair','recpair'),sep = '_',remove = F) %>%
                tidyr::separate(ligpair ,c('Ligand.Cluster','Ligand'),sep = '/',remove = F) %>%
                tidyr::separate(recpair ,c('Receptor.Cluster','Receptor'),sep = '/',remove = F)
    tmp_data$LRScore.x[is.na(tmp_data$LRScore.x)] <- 0
    tmp_data$LRScore.y[is.na(tmp_data$LRScore.y)] <- 0
    final_data <- tmp_data %>%
                  dplyr::mutate(LRScore = LRScore.x-LRScore.y) %>%
                  dplyr::mutate(cellpair = paste0(Ligand.Cluster,"_",Receptor.Cluster)) %>%
                  dplyr::select(c(Ligand,Ligand.Cluster,
                                  Receptor,
                                  Receptor.Cluster,
                                  LRScore,
                                  cellpair,
                                  ligpair,
                                  recpair,
                                  allpair))
    final_data <- final_data[final_data$LRScore!=0,]
    data@tables[[cmp_name]] <- final_data
    final <- final_data %>%
      dplyr::group_by(.data$cellpair) %>%
      dplyr::summarise(LRScore = sum(.data$LRScore))
    final <- final %>%
      tidyr::separate(.data$cellpair, c("u", "v"), sep="_",remove = F)
    raw_inter <- table(final_data$cellpair)
    freq <- table(final_data$cellpair) / max(table(final_data$cellpair))
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
