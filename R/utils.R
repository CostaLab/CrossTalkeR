#'Ranking the most interactive gene (ligand or receptor)
#'
#'@param data lrobject
#'@param out_path to save the lrobject with ranking
#'@param sel_columns columns to consider
#'@param slot slot of the networks graphs_ggi to gene cell interaction and abs
#'@import tibble
#'@import utils
#'@import dplyr
#'@return list
#'@importFrom tidyr %>%
#'@importFrom stats prcomp
#'@noRd
ranking <- function(data, out_path, sel_columns, slot = "graphs_ggi") {
  for (graph in names(slot(data, slot))) {
    if (grepl("_x_", graph)) {  # Signed Analysis
      message(graph)
      tmp_g <- slot(data, slot)[[graph]]
      comp <- igraph::components(tmp_g)
      memb <- comp$membership
      all_both <- NULL
      for (i in unique(comp$membership)) {
        tmp_memb <- igraph::V(tmp_g)[memb == i]
        subgraph <- igraph::induced.subgraph(tmp_g,
                                             tmp_memb)
        igraph::V(subgraph)$name <- tmp_memb$name
        if (is.null(all_both)) {
          all_both <- ranking_net(subgraph, mode = FALSE)
        }else {
          tmp <- ranking_net(subgraph, mode = FALSE)
          tmp <- dplyr::bind_rows(all_both, tmp)
          all_both <- tmp
        }
      }
      if (grepl("_ggi", slot)) { ## Add dirichelet
        all_both <- comparative_pagerank(data@rankings, slot, graph, all_both)
        all_both <- comparative_med(data@rankings, slot, graph, all_both)
        data@rankings[[paste0(graph, "_ggi")]] <- all_both
        all_both <- all_both[, -1]
        all_both <- all_both[, which(colSums(all_both) != 0)]
        all_both <- all_both[, which(colSums(all_both) / dim(all_both)[1] != dim(all_both)[1])]
        data@pca[[paste0(graph, "_ggi")]] <- prcomp(all_both, center = TRUE, scale = TRUE)
        rownames(data@pca[[paste0(graph, "_ggi")]]$x) <- data@rankings[[paste0(graph, "_ggi")]]$nodes
        data@pca[[paste0(graph, "_ggi")]]$x <- -data@pca[[paste0(graph, "_ggi")]]$x
        data@pca[[paste0(graph, "_ggi")]]$rotation <- -data@pca[[paste0(graph, "_ggi")]]$rotation
      }else {
        all_both <- comparative_pagerank(data@rankings, slot, graph, all_both)
        all_both <- comparative_med(data@rankings, slot, graph, all_both)
        data@rankings[[graph]] <- all_both
        all_both <- all_both[, -1]
        all_both <- all_both[, which(colSums(all_both) != 0)]
        all_both <- all_both[, sapply(all_both, var) != 0]
        message(colnames(all_both)[sapply(all_both, var) == 0])
        data@pca[[graph]] <- prcomp(all_both, center = TRUE, scale = TRUE)
        rownames(data@pca[[graph]]$x) <- data@rankings[[graph]]$nodes
        data@pca[[graph]]$x <- -data@pca[[graph]]$x
        data@pca[[graph]]$rotation <- -data@pca[[graph]]$rotation

      }
    } else { # Unsigned
      comp <- igraph::components(slot(data, slot)[[graph]])
      memb <- comp$membership
      all <- NULL
      for (i in unique(comp$membership)) {
        tmp_graph <- slot(data, slot)[[graph]]
        tmp_memb <- igraph::V(tmp_graph)[memb == i]
        subgraph <- igraph::induced.subgraph(tmp_graph,
                                             tmp_memb)
        igraph::V(subgraph)$name <- tmp_memb$name
        if (is.null(all)) {
          all <- ranking_net(subgraph)
        }else {
          tmp <- ranking_net(subgraph)
          tmp <- dplyr::bind_rows(all, tmp)
          all <- tmp
        }
      }
      if (grepl("_ggi", slot)) {
        final <- NULL
        table <- slot(data, 'tables')[[graph]]
        cls <- unique(union(table[[sel_columns[1]]], table[[sel_columns[2]]]))
        for (i in cls) {
            all.eq <- unique(union(table$ligpair[table[[sel_columns[1]]] == i], table$recpair[table[[sel_columns[2]]] == i]))
            if(length(grep("\\|R$",all.eq))>0){
                edges <- t(utils::combn(all.eq, 2))
                df <- tibble::tibble(u = edges[, 1], v = edges[, 2], MeanLR = rep(0.0, dim(edges)[1]), .name_repair = ~c('u', 'v', 'LRScore'))
                if (is.null(all)) {
                  final <- df
                }
                else {
                  final <- dplyr::bind_rows(final, df)
                }
            }
        }
        names(final) <- c('ligpair', 'recpair', 'LRScore')
        tmp_tbl = table[, c('ligpair', 'recpair', 'LRScore')]
        tmp_tbl[['LRScore']] = tmp_tbl[['LRScore']]
        all1 <- dplyr::bind_rows(tmp_tbl, final)
        tmp_net <- igraph::graph_from_data_frame(all1)
        pg <- igraph::page.rank(tmp_net, weights = igraph::edge_attr(tmp_net, 'LRScore'))
        all$Pagerank <- pg$vector[all$nodes]
        data@rankings[[paste0(graph, "_ggi")]] <- all
        all <- all[, -1]
        all <- all[, which(colSums(all) != 0)]
        all <- all[, which(colSums(all) / dim(all)[1] != dim(all)[1])]
        data@pca[[paste0(graph, "_ggi")]] <- stats::prcomp(all, center = TRUE, scale = TRUE)
        rownames(data@pca[[paste0(graph, "_ggi")]]$x) <- data@rankings[[paste0(graph, "_ggi")]]$nodes
      }else {
        pg <- igraph::page.rank(slot(data, slot)[[graph]], weights = igraph::edge_attr(slot(data, slot)[[graph]], 'LRScore'))
        all$Pagerank <- pg$vector[all$nodes]
        data@rankings[[graph]] <- all
        ver <- all %>% summarise_if(is.numeric, var)
        all <- all[, c(FALSE, as.vector(ver[1,] != 0))]
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
#'@param mode  is TRUE if is comparive mode
#'@return list
#'@import igraph
#'@importFrom tidyr %>%
#'@noRd
ranking_net <- function(graph, mode = TRUE) {
  if (!mode) {
    sub_graph <- igraph::subgraph.edges(graph, igraph::E(graph)[igraph::E(graph)$weight > 0])
    deg_in_pos <- igraph::degree(sub_graph, mode = "in")
    deg_out_pos <- igraph::degree(sub_graph, mode = "out")
    sub_graph <- igraph::subgraph.edges(graph, igraph::E(graph)[igraph::E(graph)$weight < 0])
    igraph::E(graph)$weight <- abs(E(graph)$weight)
    deg_in_neg <- igraph::degree(sub_graph, mode = "in")
    deg_out_neg <- igraph::degree(sub_graph, mode = "out")
    names <- igraph::V(graph)$name
    deg_in_pos <- deg_in_pos[names]
    deg_out_pos <- deg_out_pos[names]
    deg_in_neg <- deg_in_neg[names]
    deg_out_neg <- deg_out_neg[names]
    deg_in_pos[is.na(deg_in_pos)] <- 0
    deg_out_pos[is.na(deg_out_pos)] <- 0
    deg_in_neg[is.na(deg_in_neg)] <- 0
    deg_out_neg[is.na(deg_out_neg)] <- 0
    deg_in_pos <- deg_in_pos + 1
    deg_out_pos <- deg_out_pos + 1
    deg_in_neg <- deg_in_neg + 1
    deg_out_neg <- deg_out_neg + 1
    centrality_table <- tibble::tibble(nodes = names,
                                       'Listener' = round(deg_in_pos, 2) - round(deg_in_neg, 2),
                                       'Influencer' = round(deg_out_pos, 2) - round(deg_out_neg, 2))
    centrality_table[is.na(centrality_table)] = 0
  }else {
    deg_in_pos <- igraph::degree(graph, mode = "in", normalized = FALSE)
    deg_out_pos <- igraph::degree(graph, mode = "out", normalized = FALSE)
    bet <- igraph::betweenness(graph, weights = abs(igraph::E(graph)$weight), normalized = FALSE)
    names <- igraph::V(graph)$name
    bet <- bet[names]
    deg_in_pos <- deg_in_pos[names]
    deg_out_pos <- deg_out_pos[names]
    centrality_table <- tibble::tibble(nodes = names,
                                       'Listener' = round(deg_in_pos, 2),
                                       'Influencer' = round(deg_out_pos, 2),
                                       'Mediator' = round(bet, 2))
    centrality_table[is.na(centrality_table)] = 0
  }
  return(centrality_table)
}


#'Annotate Exclusive LR pairs (ligand or receptor)
#'
#'@param data lrobject
#'@param slot table fields
#'@param out_path save path
#'@param database annotation database
#'@param org organism to be considered
#'@import clusterProfiler
#'@import org.Hs.eg.db
#'@import org.Mm.eg.db
#'@importFrom tidyr %>%
#'@noRd
kegg_annotation <- function(data, slot, out_path, database = org.Hs.eg.db::org.Hs.eg.db, org = 'hsa', n = 100) {
  rkg <- slot(data, slot)
  for (x in names(rkg)) {
    all = list()
    for (i in names(rkg[[x]])) {
      if (i != 'nodes' & grepl('ggi', x) & !grepl('_x_', x)) {
        sel <- rkg[[x]][!grepl("tf-", rkg[[x]]$nodes),]
        top <- sel %>%
          dplyr::top_n(n, wt = sel[[i]])
        topenrich <- enrich(top$nodes, name = i, db = database, org = org)
        all[[i]] <- topenrich
      }else if (i != 'nodes' & grepl('ggi', x) & grepl('_x_', x)) {
        sel <- rkg[[x]][!grepl("tf-", rkg[[x]]$nodes),]
        top <- sel %>%
          dplyr::top_n(n, wt = sel[[i]])
        topn <- sel %>%
          dplyr::top_n(-n, wt = sel[[i]])
        topenrich <- enrich(top$nodes, name = paste0(i, ' up'), db = database, org = org)
        topnenrich <- enrich(topn$nodes, name = paste0(i, ' down'), db = database, org = org)
        all[[i]] <- dplyr::bind_rows(topenrich, topnenrich)
      }
    }
    data@annot[[x]] <- dplyr::bind_rows(all)
  }
  saveRDS(data, file.path(out_path, "LR_data_final.Rds"))
  return(data)
}


#'Ranking the most interactive gene (ligand or receptor)
#'
#'@param rankings tables lrobject
#'@param slotname slot of the networks graphs_ggi to gene cell interaction and abs
#'@param graphname graph comparison name
#'@param curr.rkg ranking table
#'@importFrom tidyr %>%
#'@import stringr
#'@return list
#'@noRd
comparative_pagerank <- function(rankings, slotname, graphname, curr.rkg) {
  p_f1 <- p_f2 <- 0.5 # prob to be at disease
  allnodes <- curr.rkg$nodes
  if (str_detect(graphname, '_filtered', negate = FALSE)){
    curr = str_remove(graphname, '_filtered')
    curr = stringr::str_split(curr, '_x_')
    p_ctr = curr[[1]][2]
    q_exp = curr[[1]][1]
  } else {
    curr = stringr::str_split(graphname, '_x_')
    p_ctr = curr[[1]][2]
    q_exp = curr[[1]][1]
  }
  if (grepl("_ggi", slotname)) {
    p <- rankings[[paste0(p_ctr, '_ggi')]]$Pagerank
    q <- rankings[[paste0(q_exp, '_ggi')]]$Pagerank
  }else {
    p <- rankings[[p_ctr]]$Pagerank[names(rankings[[q_exp]]$Pagerank)]
    q <- rankings[[q_exp]]$Pagerank
  }
  final <- tibble(p.ctr = p[allnodes], p.dis = q[allnodes], names = allnodes)
  final$p.ctr[is.na(final$p.ctr)] = 0
  final$p.dis[is.na(final$p.dis)] = 0
  alpha <- 0.01
  final$p.ctr = final$p.ctr + alpha
  final$p.dis = final$p.dis + alpha
  final$p.ctr = final$p.ctr / sum(final$p.ctr)
  final$p.dis = final$p.dis / sum(final$p.dis)
  p <- final$p.ctr
  q <- final$p.dis
  pc <- p * p_f1 + q * p_f2
  pcontrol <- (p_f1 * p) / pc
  pdisease <- (p_f2 * q) / pc
  final <- log(pdisease / pcontrol)
  curr.rkg$Pagerank <- final
  return(curr.rkg)
}

#'Delta betweenness the most interactive gene (ligand or receptor)
#'
#'@param rankings tables lrobject
#'@param slotname slot of the networks graphs_ggi to gene cell interaction and abs
#'@param graphname graph comparison name
#'@param curr.rkg ranking table
#'@importFrom tidyr %>%
#'@import stringr
#'@return list
#'@noRd
comparative_med <- function(rankings, slotname, graphname, curr.rkg) {
  allnodes <- curr.rkg$nodes
  if (str_detect(graphname, '_filtered', negate = FALSE)){
    curr = str_remove(graphname, '_filtered')
    curr = stringr::str_split(curr, '_x_')
    p_ctr = curr[[1]][2]
    q_exp = curr[[1]][1]
  } else {
    curr = stringr::str_split(graphname, '_x_')
    p_ctr = curr[[1]][2]
    q_exp = curr[[1]][1]
  }
  if (grepl("_ggi", slotname)) {
    p <- rankings[[paste0(p_ctr, '_ggi')]]$Mediator
    q <- rankings[[paste0(q_exp, '_ggi')]]$Mediator
  }else {
    p <- rankings[[p_ctr]]$Mediator[names(rankings[[q_exp]]$Mediator)]
    q <- rankings[[q_exp]]$Mediator
  }
  final <- tibble(p.ctr = p[allnodes], p.dis = q[allnodes], names = allnodes)
  final$p.ctr[is.na(final$p.ctr)] <- 0
  final$p.dis[is.na(final$p.dis)] <- 0
  curr.rkg$Mediator = final$p.dis - final$p.ctr
  return(curr.rkg)
}


#'Ranking the most interactive gene (ligand or receptor)
#'
#'@param list list of genes
#'@param name measure/subject name
#'@param org annotation database default is org.Hs.eg.db
#'@param univ annotation universe
#'@importFrom tidyr %>%
#'@import stringr
#'@import clusterProfiler
#'@return list
#'@noRd
enrich <- function(list, name, db = org.Hs.eg.db, org = 'hsa', univ = NULL) {
  lrdb <- system.file("extdata",
                      "lrDB.csv",
                      package = "CrossTalkeR")
  lr <- read.csv(lrdb)
  ## Univ to be implemented mmu

  if (org == 'hsa') {
    fgenes <- list(x = gsub("/.*", "", list), y = gsub(".*/", "", list))
    fgenes[["y"]] <- gsub("\\|.*", "", fgenes[["y"]])
    nodesentrez <- clusterProfiler::bitr(fgenes$y,
                                         fromType = "SYMBOL",
                                         toType = c("ENTREZID", "ENSEMBL"),
                                         OrgDb = db)
    univ <- clusterProfiler::bitr(unique(union(lr$ligand, lr$receptor)),
                                  fromType = "SYMBOL",
                                  toType = c("ENTREZID", "ENSEMBL"),
                                  OrgDb = db)
    enriched <- clusterProfiler::enrichKEGG(nodesentrez$ENTREZID,
                                            organism = org,
                                            universe = univ$ENTREZID)
  }else {
    fgenes <- list(x = gsub("/.*", "", list), y = gsub(".*/", "", list))
    fgenes[["y"]] <- gsub("\\|.*", "", fgenes[["y"]])
    nodesentrez <- clusterProfiler::bitr(fgenes$y,
                                         fromType = "SYMBOL",
                                         toType = c("ENTREZID", "ENSEMBL"),
                                         OrgDb = db)
    enriched <- clusterProfiler::enrichKEGG(nodesentrez$ENTREZID,
                                            organism = org)
  }
  enriched <- enriched@result
  enriched$type <- name
  return(enriched)
}


#' Evaluate Differences in the edge proportion
#'@param data datafromlian
#'@param measure intensity
#'@param out_path save path
#'@importFrom tidyr %>%
#'@import tibble dplyr rstatix
#'@return tibble
#'@noRd
fisher_test_cci <- function(data, measure, out_path, comparison = NULL) {
  if (!is.null(comparison)) {
    for (pair in comparison) {
      ctr_name <- pair[2]
      exp_name <- pair[1]
      c <- data@tables[[ctr_name]] %>%
        group_by(cellpair) %>%
        select(c(source, target, measure)) %>%
        summarise(measure = n())
      e <- data@tables[[exp_name]] %>%
        dplyr::group_by(cellpair) %>%
        dplyr::select(c(.data$source, .data$target, measure)) %>%
        dplyr::summarise(measure = n())
      joined <- merge(c, e, by.x = 'cellpair', by.y = 'cellpair', keep = 'all')
      pval <- list()
      for (j in 1:length(joined$cellpair)) {
        ctotal <- sum(joined$measure.x) - joined$measure.x[j]
        etotal <- sum(joined$measure.y) - joined$measure.y[j]
        m <- matrix(c(joined$measure.y[j], etotal, joined$measure.x[j], ctotal), nrow = 2)
        rownames(m) <- c(joined$cellpair[exp_name], 'total')
        colnames(m) <- c("EXP", "CTR")
        t <- rstatix::fisher_test(m, detailed = T, B = 1000)
        pval[[joined$cellpair[j]]] <- t
      }
      pval <- bind_rows(pval, .id = 'columns_name') %>%
        mutate(lodds = log2(estimate))
      data@stats[[paste0(exp_name, '_x_', ctr_name)]] <- pval
    }
    saveRDS(data, file.path(out_path, "LR_data_final.Rds"))
    return(data)
  } else {
    if (length(data@tables) >= 2) {
      c <- data@tables[[1]] %>%
        group_by(cellpair) %>%
        select(c(source, target, measure)) %>%
        summarise(measure = n())
      for (i in 2:length(names(data@tables))) {
        if (!str_detect(names(data@tables)[i], '_x_')) {
          e <- data@tables[[i]] %>%
            dplyr::group_by(cellpair) %>%
            dplyr::select(c(.data$source, .data$target, measure)) %>%
            dplyr::summarise(measure = n())
          joined <- merge(c, e, by.x = 'cellpair', by.y = 'cellpair', keep = 'all')
          pval <- list()
          for (j in 1:length(joined$cellpair)) {
            ctotal <- sum(joined$measure.x) - joined$measure.x[j]
            etotal <- sum(joined$measure.y) - joined$measure.y[j]
            m <- matrix(c(joined$measure.y[j], etotal, joined$measure.x[j], ctotal), nrow = 2)
            rownames(m) <- c(joined$cellpair[i], 'total')
            colnames(m) <- c("EXP", "CTR")
            t <- rstatix::fisher_test(m, detailed = T, B = 1000)
            pval[[joined$cellpair[j]]] <- t
          }
          pval <- bind_rows(pval, .id = 'columns_name') %>%
            mutate(lodds = log2(estimate))
          data@stats[[paste0(names(data@tables)[i], '_x_', names(data@tables)[1])]] <- pval
        }
      }
      saveRDS(data, file.path(out_path, "LR_data_final.Rds"))
      return(data)
    }
  }
}

#' Adding genetype to the gene names to distinguish biological function
#'@param df dataframe with interaction data
#'@import tidyr
#'@import tibble dplyr
#'@return df
#'@noRd
add_node_type <- function(df) {
  df = df %>%
    mutate(gene_A = ifelse(type_gene_A == "Ligand", paste0(gene_A, "|L"), gene_A))
  df = df %>%
    mutate(gene_A = ifelse(type_gene_A == "Receptor", paste0(gene_A, "|R"), gene_A))
  df = df %>%
    mutate(gene_A = ifelse(type_gene_A == "Transcription Factor", paste0(gene_A, "|TF"), gene_A))
  df = df %>%
    mutate(gene_B = ifelse(type_gene_B == "Ligand", paste0(gene_B, "|L"), gene_B))
  df = df %>%
    mutate(gene_B = ifelse(type_gene_B == "Receptor", paste0(gene_B, "|R"), gene_B))
  df = df %>%
    mutate(gene_B = ifelse(type_gene_B == "Transcription Factor", paste0(gene_B, "|TF"), gene_B))
  return(df)
}

#' Filter comparison graphs by statistics
#'@param data datafromlian
#'@param out_path save path
#'@importFrom tidyr %>%
#'@import tibble dplyr rstatix
#'@return tibble
#'@noRd
filtered_graphs <- function(data, out_path) {
  for (name in names(data@graphs)) {
    if (str_detect(name, '_x_', negate = FALSE)) {
      h <- head_of(data@graphs[[name]], E(data@graphs[[name]]))$name
      f <- tail_of(data@graphs[[name]], E(data@graphs[[name]]))$name
      curr_net <- subgraph.edges(data@graphs[[name]],
                                 E(data@graphs[[name]])[
                                   match(data@stats[[name]]$columns_name[data@stats[[name]]$p <= 0.05],
                                         paste(h, f, sep = '@'),
                                         nomatch = F)])
      data@graphs[[paste0(name, "_filtered")]] <- curr_net
    }
  }
  saveRDS(data, file.path(out_path, "LR_data_final.Rds"))
  return(data)
}

