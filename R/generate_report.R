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
generate_report <- function(LRpaths,genes,out_path,sep=',',threshold=50,colors=NULL,out_file=NULL){
  # Creating the single condition Object
  message('Read Files')
  data <- read_lr_single_condiction(LRpaths,out_path,sep=',',colors)
  # Obtaining the differential table
  message('Create a Differential Table')
  if(length(LRpaths)>1){
    data <- create_diff_table(data,out_path)
  }
  message('Defining templates')
  # Generating the single condition report
  single <- system.file('templates','Single_Condition.Rmd', package = 'LRAnalytics')
  comp <- system.file('templates','Comparative_Condition.Rmd', package = 'LRAnalytics')
  sankey <- system.file('templates','SankeyPlots.Rmd', package = 'LRAnalytics')
  tbl <- system.file('templates','Tables.Rmd', package = 'LRAnalytics')
  lrObj_path1 <- paste0(out_path,'LR_data.Rds')
  lrObj_path2 <- paste0(out_path,'LR_data_step2.Rds')
  param <- list(single=single,
                comp = comp,
                obj1 = lrObj_path1,
                obj2 = lrObj_path2,
                obj3 = genes,
                thr = threshold)
  index <- system.file('templates','FinalReport.Rmd', package = 'LRAnalytics')
  message('Generating Report')
  rmarkdown::render(index,
                    output_format = 'html_document',
                    output_dir = out_path,output_file = out_file,params = param)

  return(data)
}

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
    final <- final %>%
      tidyr::separate(cellpair, c("u", "v"), "_")
    for(cnt in 1:dim(final)[1]){
         final$MeanLR[cnt] =  final$MeanLR[cnt] - final2dosage[[final$u[cnt]]][[final$v[cnt]]]
    }
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
            loadings=load)
  saveRDS(LR,paste0(out_path,'/LR_data.Rds'))
  return(LR)
}

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
    final_data <- tibble::tibble(Ligand.Cluster=unlist(lclust),
                                 Receptor.Cluster=unlist(rclust),
                                 Ligand=unlist(lgene),
                                 Receptor=unlist(rgene),
                                 cellpair=unlist(cpair),
                                 ligpair=unlist(lpair),
                                 recpair=unlist(rpair),
                                 allpair=unlist(apair),
                                 MeanLR=unlist(mlr))
    data@tables[[cmp_name]] <- final_data
    final <- final_data %>%
      dplyr::group_by(cellpair) %>%
      dplyr::summarise(MeanLR = sum(MeanLR))
    aux <- final$cellpair
    final <- final %>%
      tidyr::separate(cellpair, c("u", "v"), "_")
    for(cnt in 1:dim(final)[1]){
      exp_l <- ifelse(!is.null(data@loadings[[exp_name]][[final$u[cnt]]][[final$v[cnt]]]),data@loadings[[exp_name]][[final$u[cnt]]][[final$v[cnt]]], 0.0)
      ctr_l <- ifelse(!is.null(data@loadings[[ctr_name]][[final$u[cnt]]][[final$v[cnt]]]),data@loadings[[ctr_name]][[final$u[cnt]]][[final$v[cnt]]], 0.0)
      final$MeanLR[cnt]=final$MeanLR[cnt]-as.double(exp_l)-as.double(ctr_l)
    }

    final$pair=aux
    freq = table(final_data$cellpair)/sum(table(final_data$cellpair))
    final$freq <- as.array(freq)[final$pair]
    graph1 <- igraph::graph_from_data_frame(final[,c('u','v',"MeanLR")])
    igraph::E(graph1)$inter <- final$freq*100 #setting thickness and weight
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

#'Count Gene Dosage
#'
#'@param list Paths of single condition LR data
#'@return list
#'@importFrom tidyr %>%
#'@importFrom foreach %dopar%
remove_dosage<-function(final){
  #doParallel::registerDoParallel(core=4)
  final <- final[final$MeanLR!=0,]
  loadings <- list()
  for(i in unique(final$Ligand.Cluster)){
    print(i)
    loadings[[i]] <- list()
    sel_lig = match(final$Ligand.Cluster,i,nomatch = F) != F
    for(j in unique(final$Receptor.Cluster[sel_lig])){
      loadings[[i]][[j]] <- 0
      sel = sel_lig &(match(final$Receptor.Cluster,j,nomatch = F) != F)
      # Checking Ligands Dosage
      for(k in unique(final$Ligand[sel])){
        sel_lg = sel & (match(final$Ligand,k,nomatch = F)!=F)
        sz = length(final$MeanLR[sel_lg])
        if(sz >=2){
          lr <- sum(final$MeanLR[sel_lg][1:2])
          a <- final$MeanLR[sel_lg][1]
          b <- final$MeanLR[sel_lg][2]
          diff <- a-b
          lig <- (lr-diff)/2
          loadings[[i]][[j]] <- loadings[[i]][[j]]+lig*(sz-1)
        }
      }
      #Checking Receptor Dosage
      for(k in unique(final$Receptor[sel])){
        sel_lg = sel & (match(final$Receptor,k,nomatch = F)!=F)
        sz = length(final$MeanLR[sel_lg])
        if(sz >=2){
          lr <- sum(final$MeanLR[sel_lg][1:2])
          a <- final$MeanLR[sel_lg][1]
          b <- final$MeanLR[sel_lg][2]
          diff <- a-b
          lig <- (lr-diff)/2
          loadings[[i]][[j]] <- loadings[[i]][[j]]+lig*(sz-1)
        }
      }
      print(paste0('     ',j))
    }
  }
  return(loadings)
}


#'
#'Run and Generate all LR Downstream analysis
#'
#'This function loads the single conditions LR outputs and return the LR network based analysis.
#'It assumes that the table present the following columns Ligand, Ligand.Cluster, Receptor,Receptor.Cluster and MeanLR/another
#'measure
#'@slot graphs All Cell Cell Interaction Networks
#'@slot tables All tables from single condition
#'@slot max_iter  Max meanLR from all
#'@slot max_nodes All Celltype in the experiment
#'@slot coords  Cell Cell Interaction Plots
#'@slot colors  Cell type colors
#'@slot loadings  CCI values to remove multiple times genes
LRObject <- setClass("LRObj", slots=list(graphs="list",
                                         graphs_ggi = "list",
                                         tables="list",
                                         max_iter="numeric",
                                         max_nodes="numeric",
                                         coords="array",
                                         colors="character",
                                         loadings="list"))
