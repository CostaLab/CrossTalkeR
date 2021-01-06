#'Ranking the most interactive gene (ligand or receptor)
#'
#'@param data lrobject
#'@param out_path to save the lrobject with ranking
#'@param slot slot of the networks graphs_ggi to gene cell interaction and abs
#' graphs to cell cell interactions
#'@return list
#'@importFrom tidyr %>%
ranking <- function(data, out_path, slot="graphs_ggi") {
      for (graph in names(slot(data, slot))) {
          if (grepl("_x_", graph)) {  # Signed Analysis
              tmp_g <- slot(data, slot)[[graph]]
              up_graph <- igraph::subgraph.edges(tmp_g,
                                                 E(tmp_g)[E(tmp_g)$MeanLR > 0])
              d_graph <- igraph::subgraph.edges(tmp_g,
                                                   E(tmp_g)[E(tmp_g)$MeanLR < 0]
                                                  )
              comp <- igraph::components(up_graph)
              all_up <- NULL
              for (i in unique(comp$membership)) {
                 tmp_m <- comp$membership == i
                 subgraph <- igraph::induced.subgraph(up_graph,
                                                      igraph::V(up_graph)[tmp_m]
                                                     )
                 if (is.null(all_up)) {
                     all_up <- ranking_net(subgraph)
                 }else{
                     tmp <- ranking_net(subgraph)
                     tmp <- dplyr::bind_rows(all_up, tmp)
                     all_up <- tmp
                 }
              }
              comp <- igraph::components(d_graph)
              all_down <- NULL
              for (i in unique(comp$membership)) {
                 tmp_m <- comp$membership == i
                 subgraph <- igraph::induced.subgraph(d_graph,
                                                      igraph::V(d_graph)[tmp_m]
                                                     )
                 if (is.null(all_down)) {
                    all_down <- ranking_net(subgraph)
                 }else{
                   tmp <- ranking_net(subgraph)
                   tmp <- dplyr::bind_rows(all_down, tmp)
                   all_down <- tmp
                 }
              }
              if (grepl("_ggi", slot)) {
                data@rankings[[paste0(graph, "_ggi_up")]] <-  all_up
                data@rankings[[paste0(graph, "_ggi_down")]] <-  all_down
              }else{
                data@rankings[[paste0(graph, "_up")]] <-  all_up
                data@rankings[[paste0(graph, "_down")]] <-  all_down
              }
          }else{ # Unsigned
              comp <- igraph::components(slot(data, slot)[[graph]])
              memb <- comp$membership
              all <- NULL
              for (i in unique(comp$membership)) {
                tmp_graph <- slot(data, slot)[[graph]]
                tmp_memb <-  igraph::V(tmp_graph)[memb == i]
                subgraph <- igraph::induced.subgraph(tmp_graph,
                                                     tmp_memb)
                 if (is.null(all)) {
                    all <- ranking_net(subgraph)
                 }else{
                    tmp <- ranking_net(subgraph)
                    tmp <- dplyr::bind_rows(all, tmp)
                    all <- tmp
                 }
              }
              if (grepl("_ggi", slot)) {
                data@rankings[[paste0(graph, "_ggi")]] <-  all
              }else{
                data@rankings[[graph]] <-  all
              }
          }
      }
  saveRDS(data, file.path(out_path, "LR_data_final.Rds"))
  return(data)
}




#' Network Ranking method
#'
#'@param graph lrobject
#'@return list
#'@importFrom tidyr %>%
ranking_net <- function(graph) {
  igraph::E(graph)$weight <- abs(igraph::E(graph)$weight)
  bet <- igraph::betweenness(graph)
  clo <- igraph::closeness(graph)
  eigen <- igraph::eigen_centrality(graph)$vector
  pagerank <- igraph::page.rank(graph)$vector
  deg_in <- igraph::degree(graph, mode = "in")
  deg_out <- igraph::degree(graph, mode = "out")
  centrality_table <- tibble::tibble(nodes = names(bet),
                                     indegree = deg_in,
                                     outdegree = deg_out,
                                     betweenness = bet,
                                     closeness = clo,
                                     eigenvector = eigen,
                                     pagerank = pagerank,
                                     combined_ranking = rkg_ties(indegree) +
                                                        rkg_ties(outdegree) +
                                                        rkg_ties(bet) +
                                                        rkg_ties(clo) +
                                                        rkg_ties(eigen) +
                                                        rkg_ties(pagerank)
                                    )
  return(centrality_table)
}

#'@param  lista
#'@return combined ranking with ties.method
rkg_ties <- function(lista) {
  x2 <- lista
  rma <- rank(x2, ties.method = "max")  # as used classically
  rmi <- rank(x2, ties.method = "min")  # as in Sports
  return(sort((rma + rmi) / 2, decreasing = TRUE))
}
