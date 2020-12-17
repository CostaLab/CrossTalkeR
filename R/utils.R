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
#'@param graph lrobject
#'@return list
#'@importFrom tidyr %>%
#'@importFrom foreach %dopar%
ranking_net<-function(graph){
    E(graph)$weight <- abs(E(graph)$weight)
    bet <- rkg_ties(igraph::betweenness(graph))
    clo <- rkg_ties(igraph::closeness(graph))
    eigen <- rkg_ties(igraph::eigen_centrality(graph)$vector)
    pagerank <- rkg_ties(igraph::page.rank(graph)$vector)
    ac <- igraph::V(graph) %in% igraph::articulation.points(graph)
    centrality_table <- tibble::tibble(nodes = names(bet),
                                       betweenness=bet,
                                       closeness=clo,
                                       eigenvector=eigen,
                                       pagerank=pagerank,
                                       articulatio_ptn=ac,
                                       combined_ranking=bet+
                                                        clo+
                                                        eigen+
                                                        pagerank+ac)
  return(centrality_table)
}

#'@param  lista
#'@return combined ranking with ties.method
rkg_ties <- function(lista){
  x2 <- lista
  rma <- rank(x2, ties.method= "max")  # as used classically
  rmi <- rank(x2, ties.method= "min")  # as in Sports
  return(sort(rma+rmi/2, decreasing=T))
}
