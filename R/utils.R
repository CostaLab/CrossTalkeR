#'Ranking the most interactive gene (ligand or receptor)
#'
#'@param data lrobject
#'@param out_path to save the lrobject with ranking
#'@param slot slot of the networks graphs_ggi to gene cell interaction and abs
#'@import tibble
#'@import utils
#'@import dplyr
#'@return list
#'@importFrom tidyr %>%
#'@importFrom stats prcomp
ranking <- function(data, out_path, slot="graphs_ggi") {
      for (graph in names(slot(data, slot))) {
          if (grepl("_x_", graph)) {  # Signed Analysis
              message(graph)
              tmp_g <- slot(data, slot)[[graph]]
              comp <- igraph::components(tmp_g)
              memb <- comp$membership
              all_both <- NULL
              for (i in unique(comp$membership)) {
                tmp_memb <-  igraph::V(tmp_g)[memb == i]
                subgraph <- igraph::induced.subgraph(tmp_g,
                                                     tmp_memb)
                igraph::V(subgraph)$name <- tmp_memb$name
                 if (is.null(all_both)) {
                    all_both <- ranking_net(subgraph,mode=FALSE)
                 }else{
                   tmp <- ranking_net(subgraph,mode=FALSE)
                   tmp <- dplyr::bind_rows(all_both, tmp)
                   all_both <- tmp
                 }
              }
              if (grepl("_ggi", slot)) { ## Add dirichelet
                all_both <- comparative_pagerank(data@rankings,slot,graph,all_both)
                all_both<- comparative_med(data@rankings,slot,graph,all_both)
                data@rankings[[paste0(graph, "_ggi")]] <- all_both
                all_both <- all_both[,-1]
                all_both <- all_both[,which(colSums(all_both)!=0)]
                all_both <- all_both[,which(colSums(all_both)/dim(all_both)[1]!=dim(all_both)[1])]
                data@pca[[paste0(graph, "_ggi")]] <- prcomp(all_both, center = TRUE, scale = TRUE)
                rownames(data@pca[[paste0(graph, "_ggi")]]$x) <- data@rankings[[paste0(graph, "_ggi")]]$nodes
                data@pca[[paste0(graph, "_ggi")]]$x <- -data@pca[[paste0(graph, "_ggi")]]$x
                data@pca[[paste0(graph, "_ggi")]]$rotation <- -data@pca[[paste0(graph, "_ggi")]]$rotation
              }else{
                all_both<- comparative_pagerank(data@rankings,slot,graph,all_both)
                all_both<- comparative_med(data@rankings,slot,graph,all_both)
                data@rankings[[graph]] <- all_both
                all_both <- all_both[,-1]
                all_both <- all_both[,which(colSums(all_both)!=0)]
                all_both <- all_both[,which(colSums(all_both)/dim(all_both)[1]!=dim(all_both)[1])]
                data@pca[[graph]] <- prcomp(all_both, center = TRUE, scale = TRUE)
                rownames(data@pca[[graph]]$x) <- data@rankings[[graph]]$nodes
                data@pca[[graph]]$x <- -data@pca[[graph]]$x
                data@pca[[graph]]$rotation <- -data@pca[[graph]]$rotation

              }
            } else{ # Unsigned
              comp <- igraph::components(slot(data, slot)[[graph]])
              memb <- comp$membership
              all <- NULL
              for (i in unique(comp$membership)) {
                tmp_graph <- slot(data, slot)[[graph]]
                tmp_memb <-  igraph::V(tmp_graph)[memb == i]
                subgraph <- igraph::induced.subgraph(tmp_graph,
                                                     tmp_memb)
                igraph::V(subgraph)$name <- tmp_memb$name
                if (is.null(all)) {
                    all <- ranking_net(subgraph)
                }else{
                    tmp <- ranking_net(subgraph)
                    tmp <- dplyr::bind_rows(all,tmp)
                    all <- tmp
                }
              }
              if (grepl("_ggi", slot)) {
                final <- NULL
                table <- slot(data, 'tables')[[graph]]
                cls <- unique(union(table$Ligand.Cluster,table$Receptor.Cluster))
                for(i in cls){
                    all.eq <- unique(union(table$ligpair[table$Ligand.Cluster==i],table$recpair[table$Receptor.Cluster==i]))
                    edges <- t(utils::combn(all.eq,2))
                    df <- tibble::tibble(u=edges[,1],u=edges[,2],MeanLR=rep(0.1,dim(edges)[1]),.name_repair='minimal')
                    if(is.null(all)){
                      final <- df
                    }
                    else{
                      final <- dplyr::bind_rows(final,df)
                    }
                }
                names(final) <- c('ligpair','recpair','MeanLR')
                tmp_tbl = table[,c('ligpair','recpair','MeanLR')]
                tmp_tbl$MeanLR = tmp_tbl$MeanLR+0.1
                all1 <- dplyr::bind_rows(tmp_tbl,final)
                tmp_net <- igraph::graph_from_data_frame(all1)
                pg <- igraph::page.rank(tmp_net,weights=igraph::E(tmp_net)$MeanLR)
                all$Pagerank <- pg$vector[all$nodes]
                data@rankings[[paste0(graph, "_ggi")]] <- all
                all <- all[,-1]
                all <- all[,which(colSums(all)!=0)]
                all <- all[,which(colSums(all)/dim(all)[1]!=dim(all)[1])]
                data@pca[[paste0(graph, "_ggi")]] <- stats::prcomp(all, center = TRUE, scale = TRUE)
                rownames(data@pca[[paste0(graph, "_ggi")]]$x) <- data@rankings[[paste0(graph, "_ggi")]]$nodes
              }else{
                pg <- igraph::page.rank(slot(data, slot)[[graph]],weights=igraph::E(slot(data, slot)[[graph]])$MeanLR)
                all$Pagerank <- pg$vector[all$nodes]
                data@rankings[[graph]] <-  all
                all <- all[,-1]
                all <- all[,which(colSums(all)!=0)]
                all <- all[,which(colSums(all)/dim(all)[1]!=dim(all)[1])]
                data@pca[[graph]] <- stats::prcomp(all, center = TRUE, scale = TRUE)
               rownames(data@pca[[graph]]$x) <- data@rankings[[graph]]$nodes
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
#'@import igraph
#'@importFrom tidyr %>%
ranking_net <- function(graph,mode=TRUE) {
  if(!mode){
        sub_graph <- igraph::subgraph.edges(graph, igraph::E(graph)[igraph::E(graph)$weight>0])
        deg_in_pos <- igraph::degree(sub_graph, mode = "in")
        deg_out_pos <- igraph::degree(sub_graph, mode = "out")
        sub_graph <- igraph::subgraph.edges(graph, igraph::E(graph)[igraph::E(graph)$weight<0])
        igraph::E(graph)$weight <- abs(E(graph)$weight)
        deg_in_neg <- igraph::degree(sub_graph, mode = "in")
        deg_out_neg <- igraph::degree(sub_graph, mode = "out")
        names <- igraph::V(graph)$name
        deg_in_pos  <- deg_in_pos[names]
        deg_out_pos  <- deg_out_pos[names]
        deg_in_neg  <- deg_in_neg[names]
        deg_out_neg  <- deg_out_neg[names]
        deg_in_pos[is.na(deg_in_pos)]  <- 0
        deg_out_pos[is.na(deg_out_pos)]  <-0
        deg_in_neg[is.na(deg_in_neg)]  <- 0
        deg_out_neg[is.na(deg_out_neg)]  <- 0
        deg_in_pos  <- deg_in_pos+1
        deg_out_pos  <- deg_out_pos+1
        deg_in_neg  <- deg_in_neg+1
        deg_out_neg  <- deg_out_neg+1
        centrality_table <- tibble::tibble(nodes = names,
                                           'Influenced' = round(deg_in_pos,6)-round(deg_in_neg,6),
                                           'Influencer' = round(deg_out_pos,6)-round(deg_out_neg,6))

         centrality_table[is.na(centrality_table)] = 0

  }else{
    deg_in_pos <- igraph::degree(graph, mode = "in")
    deg_out_pos <- igraph::degree(graph, mode = "out")
    bet <- igraph::betweenness(graph, weights = abs(igraph::E(graph)$weight))
    names <- igraph::V(graph)$name
    bet  <- bet[names]
    deg_in_pos  <- deg_in_pos[names]
    deg_out_pos  <- deg_out_pos[names]
    centrality_table <- tibble::tibble(nodes = names,
                                       'Influenced' = round(deg_in_pos,6),
                                       'Influencer' = round(deg_out_pos,6),
                                       'Mediator' = round(bet,6))
    centrality_table[is.na(centrality_table)] = 0

  }
  return(centrality_table)
}



#'Annotate Exclusive LR pairs (ligand or receptor)
#'
#'@param data lrobject
#'@param slot table fields
#'@param subslot table
#'@param database complex heatmap database
#'@param org organism to be considered
#'@import clusterProfiler
#'@import org.Hs.eg.db
#'@importFrom tidyr %>%
kegg_annotation <- function(data, slot,out_path,database=org.Hs.eg.db::org.Hs.eg.db, org='hsa') {
  rkg <- slot(data, slot)
  for(x in names(rkg)){
    all = list()
    for(i in names(rkg[[x]])){
        if(i != 'nodes' & grepl('ggi',x) & !grepl('_x_',x)){
          top50 <- rkg[[x]] %>%
                   dplyr::top_n(100,wt=rkg[[x]][[i]])
          top50enrich <- enrich(top50$nodes,name=i)
          all[[i]] <- top50enrich
      }else if(i != 'nodes' & grepl('ggi',x) & grepl('_x_',x)){
           top50 <- rkg[[x]] %>%
                 dplyr::top_n(100,wt=rkg[[x]][[i]])
            top50n <- rkg[[x]] %>%
                    dplyr::top_n(-100,wt=rkg[[x]][[i]])
            top50enrich <- enrich(top50$nodes,name=paste0(i,' up'))
            topn50enrich <- enrich(top50n$nodes,name=paste0(i,' down'))
            all[[i]] <- dplyr::bind_rows(top50enrich,topn50enrich)
        }
    }
    data@annot[[x]] <- dplyr::bind_rows(all)
  }
  saveRDS(data,file.path(out_path, "LR_data_final.Rds"))
  return(data)
}







#'Ranking the most interactive gene (ligand or receptor)
#'
#'@param ranking tables lrobject
#'@param slot slot of the networks graphs_ggi to gene cell interaction and abs
#'@importFrom tidyr %>%
#'@import stringr
#'@return list
comparative_pagerank<- function(rankings,slotname,graphname,curr.rkg){
  p_f1 <- 0.5 # prob to be at disease
  p_f2 <- 0.5
  allnodes <- curr.rkg$nodes
  curr = stringr::str_split(graphname,'_x_')
  p_ctr = curr[[1]][2]
  q_exp = curr[[1]][1]
  if (grepl("_ggi", slotname)){
    p <- rankings[[paste0(p_ctr,'_ggi')]]$Pagerank
    q <- rankings[[paste0(q_exp,'_ggi')]]$Pagerank
  }else{
    p <- rankings[[p_ctr]]$Pagerank[names(rankings[[q_exp]]$Pagerank)]
    q <- rankings[[q_exp]]$Pagerank
  }
  final <-tibble(p.ctr=p[allnodes], p.dis = q[allnodes], names=allnodes)
  final$p.ctr[is.na(final$p.ctr)] = 0
  final$p.dis[is.na(final$p.dis)] = 0
  alpha <- 0.01
  final$p.ctr = final$p.ctr + alpha
  final$p.dis = final$p.dis + alpha
  final$p.ctr = final$p.ctr/sum(final$p.ctr)
  final$p.dis = final$p.dis/sum(final$p.dis)
  p <- final$p.ctr
  q <- final$p.dis
  pc <- p*p_f1 + q*p_f2
  pcontrol <- (p_f1*p)/pc
  pdisease <- (p_f2*q)/pc
  final <- log(pdisease/pcontrol)
  curr.rkg$Pagerank <- final
  return(curr.rkg)
}

#'Delta betweenness the most interactive gene (ligand or receptor)
#'
#'@param ranking tables lrobject
#'@param slot slot of the networks graphs_ggi to gene cell interaction and abs
#'@importFrom tidyr %>%
#'@import stringr
#'@return list
comparative_med<- function(rankings,slotname,graphname,curr.rkg){
  allnodes <- curr.rkg$nodes
  curr = stringr::str_split(graphname,'_x_')
  p_ctr = curr[[1]][2]
  q_exp = curr[[1]][1]
  if (grepl("_ggi", slotname)){
    p <- rankings[[paste0(p_ctr,'_ggi')]]$Mediator
    q <- rankings[[paste0(q_exp,'_ggi')]]$Mediator
  }else{
    p <- rankings[[p_ctr]]$Mediator[names(rankings[[q_exp]]$Mediator)]
    q <- rankings[[q_exp]]$Mediator
  }
  final <-tibble(p.ctr=p[allnodes], p.dis = q[allnodes], names=allnodes)
  final$p.ctr[is.na(final$p.ctr)]  <- 0
  final$p.dis[is.na(final$p.dis)]  <-0
  curr.rkg$Mediator = final$p.dis-final$p.ctr
  return(curr.rkg)
}

#'Ranking the most interactive gene (ligand or receptor)
#'
#'@param list of genes
#'@param rankings's table column name
#'@importFrom tidyr %>%
#'@import stringr
#'@import clusterProfiler
#'@return list
enrich <- function(list,name,org=org.Hs.eg.db, univ=NULL){
  lrdb <- system.file("extdata",
                        "lrDB.csv",
                        package = "CrossTalkeR")
  lr <- read.csv(lrdb)
  univ <- clusterProfiler::bitr(unique(union(lr$ligand,lr$receptor)),
                                fromType="SYMBOL",
                                toType=c("ENTREZID","ENSEMBL"),
                                OrgDb=org.Hs.eg.db)
  fgenes<-list(x=gsub("/.*","",list),y=gsub(".*/","",list))
  nodesentrez <- clusterProfiler::bitr(fgenes$x,
                                       fromType="SYMBOL",
                                       toType=c("ENTREZID","ENSEMBL"),
                                       OrgDb=org)

  enriched <- clusterProfiler::enrichKEGG(nodesentrez$ENTREZID,
                                            organism = 'hsa',
                                            universe=univ$ENTREZID)
  enriched <- enriched@result
  enriched$type <-name
  return(enriched)
}
