



#'Read the LRObject and generate the comparative tables
#'
#'@param LRpaths Paths of single condition LR data
#'@param sep character used on csv
#'@return LRObject
#'@importFrom tidyr %>%
#'@export
create_diff_table <- function(data,out_path){
  ctr_name <- names(data@tables)[1]
  ctr_table <- data@tables[[ctr_name]]
  for(i in 2:length(data@tables)){
    exp_name <-names(data@tables)[i]
    cmp_name <- paste0(exp_name,'_x_',ctr_name)
    exp_table <- data@tables[[exp_name]]
    lclust <- list()
    rclust <- list()
    lgene <- list()
    rgene <- list()
    lpair <- list()
    cpair <- list()
    rpair <- list()
    apair <- list()
    mlr <- list()
    for(m in exp_table$allpair){
      if(m %in% ctr_table$allpair){
        idx_e <- match(m,exp_table$allpair)
        idx_c <- match(m,ctr_table$allpair)
        lclust[[m]] <- exp_table[idx_e,]$Ligand.Cluster
        rclust[[m]] <- exp_table[idx_e,]$Receptor.Cluster
        lgene[[m]] <- exp_table[idx_e,]$Ligand
        rgene[[m]] <- exp_table[idx_e,]$Receptor
        lpair[[m]] <- exp_table[idx_e,]$ligpair
        cpair[[m]] <- exp_table[idx_e,]$cellpair
        rpair[[m]] <- exp_table[idx_e,]$recpair
        apair[[m]] <- exp_table[idx_e,]$allpair
        mlr[[m]] <- exp_table[idx_e,]$MeanLR-ctr_table[idx_c,]$MeanLR
      }
      else{
        idx_e <- match(m,exp_table$allpair)
        lclust[[m]] <- exp_table[idx_e,]$Ligand.Cluster
        rclust[[m]] <- exp_table[idx_e,]$Receptor.Cluster
        lgene[[m]] <- exp_table[idx_e,]$Ligand
        rgene[[m]] <- exp_table[idx_e,]$Receptor
        lpair[[m]] <- exp_table[idx_e,]$ligpair
        cpair[[m]] <- exp_table[idx_e,]$cellpair
        rpair[[m]] <- exp_table[idx_e,]$recpair
        apair[[m]] <- exp_table[idx_e,]$allpair
        mlr[[m]] <- exp_table[idx_e,]$MeanLR
      }
    }
    for(m in ctr_table$allpair){
      if(!m %in% exp_table$allpair){
        idx_c <- match(m,ctr_table$allpair)
        lclust[[m]] <- ctr_table[idx_c,]$Ligand.Cluster
        rclust[[m]] <- ctr_table[idx_c,]$Receptor.Cluster
        lgene[[m]] <- ctr_table[idx_c,]$Ligand
        rgene[[m]] <- ctr_table[idx_c,]$Receptor
        lpair[[m]] <- ctr_table[idx_c,]$ligpair
        cpair[[m]] <- ctr_table[idx_c,]$cellpair
        rpair[[m]] <- ctr_table[idx_c,]$recpair
        apair[[m]] <- ctr_table[idx_c,]$allpair
        mlr[[m]] <- 0-ctr_table[idx_c,]$MeanLR
      }
    }
    final_data <- tibble::tibble('Ligand.Cluster'=unlist(lclust),
                                 'Receptor.Cluster'=unlist(rclust),
                                 'Ligand'=unlist(lgene),
                                 'Receptor'=unlist(rgene),
                                 'cellpair'=unlist(cpair),
                                 'ligpair'=unlist(lpair),
                                 'recpair'=unlist(rpair),
                                 'allpair'=unlist(apair),
                                 'MeanLR'=unlist(mlr))
    data@tables[[cmp_name]] <- final_data
    final <- final_data %>%
      dplyr::group_by(.data$cellpair) %>%
      dplyr::summarise(MeanLR=sum(.data$MeanLR))
    aux <- final$cellpair
    final <- final %>%
      tidyr::separate(.data$cellpair, c("u", "v"), "_")
    #clusters_num <- unique(c(unique(final_data$Ligand.Cluster),unique(final_data$Receptor.Cluster)))
    #for(cnt in 1:dim(final)[1]){
    # exp_l <- ifelse(!is.null(data@loadings[[exp_name]][[final$u[cnt]]][[final$v[cnt]]]),data@loadings[[exp_name]][[final$u[cnt]]][[final$v[cnt]]], 0.0)
    # ctr_l <- ifelse(!is.null(data@loadings[[ctr_name]][[final$u[cnt]]][[final$v[cnt]]]),data@loadings[[ctr_name]][[final$u[cnt]]][[final$v[cnt]]], 0.0)
    # final$MeanLR[cnt]=(final$MeanLR[cnt]+exp_l+ctr_l)/(length(clusters_num)*length(clusters_num))
    # final$MeanLR[cnt]=(final$MeanLR[cnt])/(length(clusters_num)*length(clusters_num))
    #}
    final$pair=aux
    final_data <- final_data[final_data$MeanLR!=0,]
    freq = table(final_data$cellpair)/max(table(final_data$cellpair))
    final$freq <- as.array(freq)[final$pair]
    graph1 <- igraph::graph_from_data_frame(final[,c('u','v',"MeanLR")])
    igraph::E(graph1)$inter <- final$freq #setting thickness and weight
    igraph::E(graph1)$weight <- igraph::E(graph1)$MeanLR
    igraph::E(graph1)$mean <- igraph::E(graph1)$MeanLR
    data@graphs[[cmp_name]] <- graph1
    data@tables[[cmp_name]] <- final_data
    graph1 <- igraph::graph_from_data_frame(final_data[,c('ligpair','recpair',"MeanLR")])
    igraph::E(graph1)$weight <- igraph::E(graph1)$MeanLR
    igraph::E(graph1)$mean <- igraph::E(graph1)$MeanLR
    data@graphs_ggi[[cmp_name]] <- graph1

  }
  saveRDS(data,paste0(out_path,'/LR_data_step2.Rds'))

  return(data)
}
