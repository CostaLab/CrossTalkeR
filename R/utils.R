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
    # print(i)
    loadings[[i]] <- list()
    sel_lig = match(final$Ligand.Cluster,i,nomatch = F) != F
    for(j in unique(final$Receptor.Cluster[sel_lig])){
      loadings[[i]][[j]] <- 0
      sel_cci = sel_lig &(match(final$Receptor.Cluster,j,nomatch = F) != F)
      # Checking Ligands Dosage
      for(k in unique(final$Ligand[sel_cci])){
        sel_lg = sel_cci & (match(final$Ligand,k,nomatch = F)!=F)
        sz = length(final$MeanLR[sel_lg])
        if(sz >=2){
          lr <- sum(final$MeanLR[sel_lg][1:2])
          a <- final$MeanLR[sel_lg][1]
          b <- final$MeanLR[sel_lg][2]
          diff <- a-b
          lig <- (lr-diff)
          loadings[[i]][[j]] <- loadings[[i]][[j]]-(lig*(sz-1))

        }
      }
      #Checking Receptor Dosage
      for(k in unique(final$Receptor[sel_cci])){
        sel_lg = sel_cci & (match(final$Receptor,k,nomatch = F)!=F)
        sz = length(final$MeanLR[sel_lg])
        if(sz >=2){
          lr <- sum(final$MeanLR[sel_lg][1:2])
          a <- final$MeanLR[sel_lg][1]
          b <- final$MeanLR[sel_lg][2]
          diff <- a-b
          lig <- (lr-diff)
          loadings[[i]][[j]] <- loadings[[i]][[j]]-(lig*(sz-1))
        }
      }
      #print(paste0('     ',j))
    }
  }
  return(loadings)
}


#'Ranking the most interactive gene (ligand or receptor)
#'
#'@param data lrobject
#'@param out path to save the lrobject with ranking
#'@return list
#'@importFrom tidyr %>%
#'@importFrom foreach %dopar%
ranking<-function(data,out_path,slot='graphs_ggi'){
      for(graph in names(slot(data,slot))){
          if(grepl('_x_',graph)){  # Signed Analysis
              up_graph <- igraph::subgraph.edges(slot(data,slot)[[graph]], E(slot(data,slot)[[graph]])[E(slot(data,slot)[[graph]])$MeanLR > 0])
              down_graph <- igraph::subgraph.edges(slot(data,slot)[[graph]], E(slot(data,slot)[[graph]])[E(slot(data,slot)[[graph]])$MeanLR < 0])
              comp <- igraph::components(up_graph)
              all_up <- NULL
              for(i in unique(comp$membership)){
                 subgraph <- igraph::induced.subgraph(up_graph, igraph::V(up_graph)[comp$membership==i])
                 if(is.null(all_up)){
                     all_up <- ranking_net(subgraph)
                 }else{
                     tmp <- ranking_net(subgraph)
                     tmp <- dplyr::bind_rows(all_up,tmp)
                     all_up <-tmp
                 }
              }
              comp <- igraph::components(down_graph)
              all_down <- NULL
              for(i in unique(comp$membership)){
                 subgraph <- igraph::induced.subgraph(down_graph, igraph::V(down_graph)[comp$membership==i])
                 if(is.null(all_down)){
                    all_down <- ranking_net(subgraph)
                 }else{
                   tmp <- ranking_net(subgraph)
                   tmp <- dplyr::bind_rows(all_down,tmp)
                   all_down <-tmp
                 }
              }
              if(grepl('_ggi',slot)){
                data@rankings[[paste0(graph,'_ggi_up')]] <-  all_up
                data@rankings[[paste0(graph,'_ggi_down')]] <-  all_down
              }else{
                data@rankings[[paste0(graph,'_up')]] <-  all_up
                data@rankings[[paste0(graph,'_down')]] <-  all_down
              }
          }else{ # Unsigned
              comp <- igraph::components(slot(data,slot)[[graph]])
              all <- NULL
              for(i in unique(comp$membership)){
                 subgraph <- igraph::induced.subgraph(slot(data,slot)[[graph]], igraph::V(slot(data,slot)[[graph]])[comp$membership==i])
                 if(is.null(all)){
                    all <- ranking_net(subgraph)
                 }else{
                    tmp <- ranking_net(subgraph)
                    tmp <- dplyr::bind_rows(all,tmp)
                    all <-tmp
                 }
              }
              if(grepl('_ggi',slot)){
                data@rankings[[paste0(graph,'_ggi')]] <-  all
              }else{
                data@rankings[[graph]] <-  all
              }
          }
      }
  saveRDS(data,paste0(out_path,'/LR_data_final.Rds'))
  return(data)
}




#' Network Ranking method
#'
#'@param data lrobject
#'@param out path to save the lrobject with ranking
#'@return list
#'@importFrom tidyr %>%
#'@importFrom foreach %dopar%
ranking_net<-function(graph){
    E(graph)$weight <- abs(E(graph)$weight)
    bet <- igraph::betweenness(graph)
    clo <- igraph::closeness(graph)
    eigen <- igraph::eigen_centrality(graph)
    pagerank <- igraph::page.rank(graph)
    ac <- igraph::V(graph) %in% igraph::articulation.points(graph)
    centrality_table <- tibble::tibble(nodes = names(bet),
                                       betweenness=bet,
                                       closeness=clo,
                                       eigenvector=eigen$vector,
                                       pagerank=pagerank$vector,
                                       articulatio_ptn=ac)
  return(centrality_table)
}



#'Ranking the most interactive cell type
#'
#'@param data lrobject
#'@param out path to save the lrobject with ranking
#'@return list
#'@importFrom tidyr %>%
#'@importFrom foreach %dopar%
ranking_cci<-function(data,out){
  for(graph in names(data@graphs)){
    in_deg <- table(data@tables[[graph]]$Ligand.Cluster)
    in_deg <- tibble::tibble(cells=names(in_deg),inter=in_deg)
    out_deg <-table(data@tables[[graph]]$Receptor.Cluster)
    out_deg <- tibble::tibble(cells=names(out_deg),inter=out_deg)
    bet <- igraph::betweenness(data@graphs[[graph]],weights = abs(igraph::E(data@graphs[[graph]])$MeanLR))
    clo <- igraph::closeness(data@graphs[[graph]],weights = abs(igraph::E(data@graphs[[graph]])$MeanLR))
    eigen <- igraph::eigen_centrality(data@graphs[[graph]])
    pagerank <- igraph::page.rank(data@graphs[[graph]])
    lig_order <- rank(-in_deg$inter, ties.method= "first")[match(names(data@colors),in_deg$cells)]
    rec_order <- rank(-out_deg$inter, ties.method= "first")[match(names(data@colors),out_deg$cells)]
    ac <- igraph::V(data@graphs[[graph]]) %in% igraph::articulation.points(data@graphs[[graph]])
    data@rankings[[graph]] <- tibble::tibble(nodes = names(bet),
                                             betweenness=bet,
                                             closeness=clo,
                                             eigenvector=eigen$vector,
                                             pagerank=pagerank$vector,
                                             ligand_count=lig_order,
                                             receptor_count=rec_order,
                                             articulatio_ptn=ac)
  }
  saveRDS(data,paste0(out,'/LR_data_final.Rds'))

  return(data)
}

#'Ranking the most interactive gene (ligand or receptor)
#'
#'@param data lrobject
#'@param out path to save the lrobject with ranking
#'@return list
#'@importFrom tidyr %>%
#'@importFrom foreach %dopar%
ranking_ggi<-function(data,out_path){
  for(graph in names(data@graph_ggi)){
    comp <- igraph::components(data@graphs_ggi[[graph]])
    all <- tibble::tibble()
    for(i in unique(comp$membership)){
      tmp <- tibble::tibble()
      subgraph <- igraph::induced.subgraph(data@graphs_ggi[[graph]], igraph::V(data@graphs_ggi[[graph]])[comp$membership==i])
      bet <- igraph::betweenness(subgraph,weights = abs(igraph::E(subgraph)$MeanLR),normalized = T)
      clo <- igraph::closeness(subgraph,weights = abs(igraph::E(subgraph)$MeanLR),normalized = T)
      eigen <- igraph::eigen_centrality(subgraph,scale = T,weights = abs(igraph::E(subgraph)$MeanLR))
      pagerank <- igraph::page.rank(subgraph)
      ac <- ifelse(igraph::V(subgraph)$name %in%  igraph::articulation.points(subgraph)$name, T,F)
      names(ac) <- igraph::V(subgraph)$name
      in_deg <- igraph::degree(subgraph, mode = 'in',normalized = T)
      out_deg <- igraph::degree(subgraph, mode = 'out',normalized = T)
      tmp <- tibble::tibble(nodes = names(bet),
                            betweenness=bet,
                            closeness=clo[names(bet)],
                            eigenvector=eigen$vector[names(bet)],
                            pagerank=pagerank$vector[names(bet)],
                            articulatio_ptn=ac[names(bet)],
                            indegree=in_deg[names(bet)],
                            outdegre=out_deg[names(bet)])
      all <- rbind(all,tmp)
    }
  }
  data@rankings[[paste0(graph,'_ggi')]] <-  all
  saveRDS(data,paste0(out_path,'/LR_data_final.Rds'))
  return(data)
}
