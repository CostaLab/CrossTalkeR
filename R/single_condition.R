#'Read single condition tables
#'
#'This function loads the single conditions LR outputs
#'
#'It assumes that the table present the following columns Ligand, Ligand.Cluster, Receptor,Receptor.Cluster and MeanLR/another
#'measure
#'
#'@param LRpaths Paths of single condition LR data
#'@param sep character used on csv
#'@return LRObject
#'@importFrom tidyr %>%
#'@export
read_lr_single_condiction <- function(LRpaths,out_path,sep=',',colors=NULL){
  data <- list()
  graphs <- list()
  load <- list()
  conds <- names(LRpaths)
  max <- 0
  max_nodes <- 0
  unif_celltypes <- c()
  for(i in 1:length(LRpaths)){
    print(LRpaths[i])
    data1 <- read.csv(LRpaths[i],sep=sep) # Reading csv
    data1 <- data1[, c('Ligand.Cluster','Receptor.Cluster','Ligand','Receptor', 'MeanLR')]
    data1$cellpair <- paste(data1$Ligand.Cluster,data1$Receptor.Cluster,sep='_')
    data1$ligpair <- paste(data1$Ligand,data1$Ligand.Cluster,sep='_')
    data1$recpair <- paste(data1$Receptor,data1$Receptor.Cluster,sep='_')
    data1$allpair <- paste(data1$ligpair,data1$recpair,sep='/')
    unif_celltypes <- unique(c(data1$Ligand.Cluster, data1$Receptor.Cluster,unif_celltypes))

    data1 <- tibble::as_tibble(data1)
    final <- data1 %>%
      dplyr::group_by(cellpair) %>%
      dplyr::summarise(MeanLR = sum(MeanLR))
    final2dosage <- data1 %>%
      dplyr::group_by(cellpair) %>% remove_dosage()
    load[[conds[i]]] <- final2dosage
    aux <- final$cellpair
    clusters_num <- unique(c(unique(data1$Ligand.Cluster),unique(data1$Receptor.Cluster)))
    print(length(clusters_num))
    final <- final %>%
      tidyr::separate(cellpair, c("u", "v"), "_")
    #for(cnt in 1:dim(final)[1]){
    #     final$MeanLR[cnt] =  (final$MeanLR[cnt] + final2dosage[[final$u[cnt]]][[final$v[cnt]]])/(length(clusters_num)*length(clusters_num))
    #     final$MeanLR[cnt] =  (final$MeanLR[cnt])/(length(clusters_num)*length(clusters_num))
    #}
    final$pair=aux
    freq = table(data1$cellpair)/sum(table(data1$cellpair))
    final$freq <- as.array(freq)[final$pair]
    graph1 <- igraph::graph_from_data_frame(final[,c('u','v',"MeanLR")])
    igraph::E(graph1)$inter <- final$freq*100 #setting thickness and weight
    igraph::E(graph1)$weight <- igraph::E(graph1)$MeanLR
    igraph::E(graph1)$mean <- igraph::E(graph1)$MeanLR
    data[[conds[i]]] <- data1
    graphs[[conds[i]]] <- graph1
    if(max(igraph::E(graph1)$mean) > max){
      max <- max(igraph::E(graph1)$mean)
    }
    if(length(igraph::V(graph1))> max_nodes){
      max_nodes <- length(igraph::V(graph1))
    }
  }
  template <- igraph::make_full_graph(n=max_nodes, directed = T, loops=T)
  c <- igraph::layout.circle(template)
  if(is.null(colors)){
    colors <- colorRampPalette(RColorBrewer::brewer.pal(12,"Paired"))(max_nodes)
    names(colors) <- sort(igraph::V(graphs[[names(graphs)[1]]])$name)
  }
  for(g in names(graphs)){
    sel = match(unif_celltypes,unique(igraph::V(graphs[[g]])$name), nomatch = F)==0
    if(sum(sel)!=0){
      nodes <- 1:length(unif_celltypes[sel])
      names(nodes) <-unif_celltypes[sel]
      graphs[[g]] <- igraph::add.vertices(graphs[[g]],length(nodes), attr=list(name=names(nodes)))
    }
  }
  rownames(c) <- sort(igraph::V(graphs[[names(graphs)[1]]])$name)
  LR <- new("LRObj",graphs=graphs,
            tables=data,
            max_iter=max,
            max_nodes=max_nodes,
            coords=c,
            colors = colors,
            loadings=load,
            rankings=list())
  saveRDS(LR,paste0(out_path,'/LR_data.Rds'))
  return(LR)
}
