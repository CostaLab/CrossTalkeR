
#'Read the lrobject and generate the comparative tables
#'
#'@param data LRObj with single condition
#'@param out_path output path
#'@return LRObject
#'@importFrom tidyr %>%
create_diff_table1 <- function(data, out_path) {
  ctr_name <- names(data@tables)[1]
  ctr_table <- data@tables[[ctr_name]]

  for (i in 2:length(data@tables)) {
    exp_name <- names(data@tables)[i]
    cmp_name <- paste0(exp_name, "_x_", ctr_name)
    exp_table <- data@tables[[exp_name]]
    possible_nodes <- union(ctr_table$allpair,exp_table$allpair)
    tbl_len <- length(possible_nodes)
    lclust <- rclust <- lgene <- rgene <- lpair <- cpair <- rpair <- apair <-
      vector("character", length = tbl_len)
    mlr <- vector("numeric", length = tbl_len)
    for(k in 1:tbl_len){
      m=possible_nodes[k]
      if (m %in% ctr_table$allpair) {
        idx_e <- match(m, exp_table$allpair)
        idx_c <- match(m, ctr_table$allpair)
        lclust[k] <- exp_table[idx_e, ]$Ligand.Cluster
        rclust[k] <- exp_table[idx_e, ]$Receptor.Cluster
        lgene[k] <- exp_table[idx_e, ]$Ligand
        rgene[k] <- exp_table[idx_e, ]$Receptor
        lpair[k] <- exp_table[idx_e, ]$ligpair
        cpair[k] <- exp_table[idx_e, ]$cellpair
        rpair[k] <- exp_table[idx_e, ]$recpair
        apair[k] <- exp_table[idx_e, ]$allpair
        mlr[k] <- exp_table[idx_e, ]$MeanLR - ctr_table[idx_c, ]$MeanLR
      }
      else{
        idx_e <- match(m, exp_table$allpair)
        lclust[k] <- exp_table[idx_e, ]$Ligand.Cluster
        rclust[k] <- exp_table[idx_e, ]$Receptor.Cluster
        lgene[k] <- exp_table[idx_e, ]$Ligand
        rgene[k] <- exp_table[idx_e, ]$Receptor
        lpair[k] <- exp_table[idx_e, ]$ligpair
        cpair[k] <- exp_table[idx_e, ]$cellpair
        rpair[k] <- exp_table[idx_e, ]$recpair
        apair[k] <- exp_table[idx_e, ]$allpair
        mlr[k] <- exp_table[idx_e, ]$MeanLR
      }
      if (!m %in% exp_table$allpair) {
        idx_c <- match(m, ctr_table$allpair)
        lclust[k] <- ctr_table[idx_c, ]$Ligand.Cluster
        rclust[k] <- ctr_table[idx_c, ]$Receptor.Cluster
        lgene[k] <- ctr_table[idx_c, ]$Ligand
        rgene[k] <- ctr_table[idx_c, ]$Receptor
        lpair[k] <- ctr_table[idx_c, ]$ligpair
        cpair[k] <- ctr_table[idx_c, ]$cellpair
        rpair[k] <- ctr_table[idx_c, ]$recpair
        apair[k] <- ctr_table[idx_c, ]$allpair
        mlr[k] <- 0 - ctr_table[idx_c, ]$MeanLR
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
    final_data <- final_data[final_data$MeanLR!=0,]
    data@tables[[cmp_name]] <- final_data
    final <- final_data %>%
      dplyr::group_by(.data$cellpair) %>%
      dplyr::summarise(MeanLR = sum(.data$MeanLR))
    aux <- final$cellpair
    final <- final %>%
      tidyr::separate(.data$cellpair, c("u", "v"), "_")
    final$pair <- aux
    raw_inter <- table(final_data$cellpair)
    freq <- table(final_data$cellpair) / max(table(final_data$cellpair))
    final$freq <- as.array(freq)[final$pair]
    final <- dplyr::arrange(final, abs(final$MeanLR))
    graph1 <- igraph::graph_from_data_frame(final[, c("u", "v", "MeanLR")])
    igraph::E(graph1)$inter <- final$freq #setting thickness and weight
    igraph::E(graph1)$inter.raw <- as.array(raw_inter)[final$pair] #setting thickness and weight
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
  saveRDS(data,file.path(out_path, "LR_data_step2.Rds"))
  return(data)
}



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
                tidyr::separate(ligpair ,c('Ligand','Ligand.Cluster'),sep = '/',remove = F) %>%
                tidyr::separate(recpair ,c('Receptor','Receptor.Cluster'),sep = '/',remove = F)
    tmp_data$MeanLR.x[is.na(tmp_data$MeanLR.x)] <- 0
    tmp_data$MeanLR.y[is.na(tmp_data$MeanLR.y)] <- 0
    final_data <- tmp_data %>%
                  dplyr::mutate(MeanLR = MeanLR.x-MeanLR.y) %>%
                  dplyr::mutate(cellpair = paste0(Ligand.Cluster,"_",Receptor.Cluster)) %>%
                  dplyr::select(c(Ligand,Ligand.Cluster,
                           Receptor,
                           Receptor.Cluster,
                           MeanLR,
                           cellpair,
                           ligpair,
                           recpair,
                           allpair))
    final_data <- final_data[final_data$MeanLR!=0,]
    data@tables[[cmp_name]] <- final_data
    final <- final_data %>%
      dplyr::group_by(.data$cellpair) %>%
      dplyr::summarise(MeanLR = sum(.data$MeanLR))
    final <- final %>%
      tidyr::separate(.data$cellpair, c("u", "v"), sep="_",remove = F)
    raw_inter <- table(final_data$cellpair)
    freq <- table(final_data$cellpair) / max(table(final_data$cellpair))
    final$freq <- as.array(freq)[final$cellpair]
    final$pair <- final$cellpair
    final <- dplyr::arrange(final, abs(final$MeanLR))
    print(final)
    graph1 <- igraph::graph_from_data_frame(final[, c("u", "v", "MeanLR")])
    igraph::E(graph1)$inter <- final$freq #setting thickness and weight
    igraph::E(graph1)$inter.raw <- as.array(raw_inter)[final$cellpair] #setting thickness and weight
    igraph::E(graph1)$weight <- igraph::E(graph1)$MeanLR
    igraph::E(graph1)$mean <- igraph::E(graph1)$MeanLR
    data@graphs[[cmp_name]] <- graph1
    graph1 <- igraph::graph_from_data_frame(final_data[, c("ligpair",
                                                           "recpair",
                                                           "MeanLR")])
    igraph::E(graph1)$weight <- igraph::E(graph1)$MeanLR
    igraph::E(graph1)$mean <- igraph::E(graph1)$MeanLR
    data@graphs_ggi[[cmp_name]] <- graph1
  }
  saveRDS(data,file.path(out_path, "LR_data_step2.Rds"))
  return(data)
}
