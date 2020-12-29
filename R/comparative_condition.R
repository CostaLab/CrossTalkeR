



#'Read the lrobject and generate the comparative tables
#'
#'@param lrpaths Paths of single condition LR data
#'@param sep character used on csv
#'@return LRObject
#'@importFrom tidyr %>%
#'@export
create_diff_table <- function(data, out_path) {
  ctr_name <- names(data@tables)[1]
  ctr_table <- data@tables[[ctr_name]]
  for (i in 2:length(data@tables)) {
    exp_name <- names(data@tables)[i]
    cmp_name <- paste0(exp_name, "_x_", ctr_name)
    exp_table <- data@tables[[exp_name]]
    exp_tbl_len = length(exp_table$allpair)
    lclust <- vector("character", length = exp_tbl_len)
    rclust <- vector("character", length = exp_tbl_len)
    lgene <- vector("character", length = exp_tbl_len)
    rgene <- vector("character", length = exp_tbl_len)
    lpair <- vector("character", length = exp_tbl_len)
    cpair <- vector("character", length = exp_tbl_len)
    rpair <- vector("character", length = exp_tbl_len)
    apair <- vector("character", length = exp_tbl_len)
    mlr <- vector("numeric", length = exp_tbl_len)


    for (i in 1:length(exp_table$allpair)) {
      m = exp_table$allpair[i]
      if (m %in% ctr_table$allpair) {
        idx_e <- match(m, exp_table$allpair)
        idx_c <- match(m, ctr_table$allpair)
        lclust[i] <- exp_table[idx_e, ]$Ligand.Cluster
        rclust[i] <- exp_table[idx_e, ]$Receptor.Cluster
        lgene[i] <- exp_table[idx_e, ]$Ligand
        rgene[i] <- exp_table[idx_e, ]$Receptor
        lpair[i] <- exp_table[idx_e, ]$ligpair
        cpair[i] <- exp_table[idx_e, ]$cellpair
        rpair[i] <- exp_table[idx_e, ]$recpair
        apair[i] <- exp_table[idx_e, ]$allpair
        mlr[i] <- exp_table[idx_e, ]$MeanLR - ctr_table[idx_c, ]$MeanLR
      }
      else{
        idx_e <- match(m, exp_table$allpair)
        lclust[i] <- exp_table[idx_e, ]$Ligand.Cluster
        rclust[i] <- exp_table[idx_e, ]$Receptor.Cluster
        lgene[i] <- exp_table[idx_e, ]$Ligand
        rgene[i] <- exp_table[idx_e, ]$Receptor
        lpair[i] <- exp_table[idx_e, ]$ligpair
        cpair[i] <- exp_table[idx_e, ]$cellpair
        rpair[i] <- exp_table[idx_e, ]$recpair
        apair[i] <- exp_table[idx_e, ]$allpair
        mlr[i] <- exp_table[idx_e, ]$MeanLR
      }
    }
    for (m in ctr_table$allpair) {
      if (!m %in% exp_table$allpair) {
        idx_c <- match(m, ctr_table$allpair)
        lclust[i] <- ctr_table[idx_c, ]$Ligand.Cluster
        rclust[i] <- ctr_table[idx_c, ]$Receptor.Cluster
        lgene[i] <- ctr_table[idx_c, ]$Ligand
        rgene[i] <- ctr_table[idx_c, ]$Receptor
        lpair[i] <- ctr_table[idx_c, ]$ligpair
        cpair[i] <- ctr_table[idx_c, ]$cellpair
        rpair[i] <- ctr_table[idx_c, ]$recpair
        apair[i] <- ctr_table[idx_c, ]$allpair
        mlr[i] <- 0 - ctr_table[idx_c, ]$MeanLR
      }
    }
    final_data <- tibble::tibble(
      "Ligand.Cluster" = lclust,
      "Receptor.Cluster" = rclust,
      "Ligand" = lgene,
      "Receptor" = rgene,
      "cellpair" = cpair,
      "ligpair" = lpair,
      "recpair" = rpair,
      "allpair" = apair,
      "MeanLR" = mlr
    )
    data@tables[[cmp_name]] <- final_data
    final <- final_data %>%
      dplyr::group_by(.data$cellpair) %>%
      dplyr::summarise(MeanLR = sum(.data$MeanLR))
    aux <- final$cellpair
    final <- final %>%
      tidyr::separate(.data$cellpair, c("u", "v"), "_")
    final$pair <- aux
    final <- final[final$MeanLR != 0, ]
    freq <- table(final_data$cellpair) / max(table(final_data$cellpair))
    final$freq <- as.array(freq)[final$pair]
    final <- dplyr::arrange(final, abs(final$MeanLR))
    graph1 <- igraph::graph_from_data_frame(final[, c("u", "v", "MeanLR")])
    igraph::E(graph1)$inter <- final$freq #setting thickness and weight
    igraph::E(graph1)$weight <- igraph::E(graph1)$MeanLR
    igraph::E(graph1)$mean <- igraph::E(graph1)$MeanLR
    data@graphs[[cmp_name]] <- graph1
    data@tables[[cmp_name]] <- final_data
    graph1 <- igraph::graph_from_data_frame(final_data[, c("ligpair",
                                                           "recpair",
                                                           "MeanLR")])
    igraph::E(graph1)$weight <- igraph::E(graph1)$MeanLR
    igraph::E(graph1)$mean <- igraph::E(graph1)$MeanLR
    data@graphs_ggi[[cmp_name]] <- graph1

  }
  saveRDS(data,file.path(out_path, 'LR_data_step2.Rds'))
  return(data)
}
