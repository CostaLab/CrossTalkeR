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
ranking_cci<-function(data,out_path){
  for(graph in names(data@graphs)){
    in_deg <- table(data@tables[[graph]]$Ligand.Cluster)
    in_deg <- tibble::tibble(cells=names(in_deg),inter=in_deg)
    out_deg <-table(data@tables[[graph]]$Receptor.Cluster)
    out_deg <- tibble::tibble(cells=names(out_deg),inter=out_deg)
    bet <- igraph::betweenness(data@graphs[[graph]],weights = abs(E(data@graphs[[graph]])$MeanLR))
    clo <- igraph::closeness(data@graphs[[graph]],weights = abs(E(data@graphs[[graph]])$MeanLR))
    eigen <- igraph::eigen_centrality(data@graphs[[graph]])
    pagerank <- igraph::page.rank(data@graphs[[graph]])
    lig_order <- rank(-in_deg$inter, ties.method= "first")[match(names(data@colors),in_deg$cells)]
    rec_order <- rank(-out_deg$inter, ties.method= "first")[match(names(data@colors),out_deg$cells)]
    ac <- igraph::V(data@graphs[[graph]]) %in% igraph::articulation.points(data@graphs[[graph]])
    data@rankings[[graph]] <- tibble::tibble(nodes = names(bet),betweenness=bet,closeness=clo,eigenvector=eigen$vector,pagerank=pagerank$vector,ligand_count=lig_order,receptor_count=rec_order,articulatio_ptn=ac)
  }
  saveRDS(data,paste0(out_path,'/LR_data_final.Rds'))

  return(data)
}

