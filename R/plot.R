#'Plot Cell Cell Interaction
#'
#'This function do a CCI plot
#'
#'@param graph Paths of single condition LR data
#'@param colors Cell type (Cluster) Colors
#'@param plt_name Plot Name (Title)
#'@param coords object coordinates
#'@param emax Max MeanLR across the all inputs, if its not defined,
#'the method going to consider the max find within a sample
#'@param leg Set color legend
#'@param low Lower threshold: This parameter low and high defines the edges
#'@param high Higher threshould
#' which will be filtered. Edges within the interval \[low\,high\] are filtered.
#'@param ignore_alpha not include transparency on the plot
#'@param log logscale the interactions
#'@param efactor edge scale factor
#'@param vfactor edge scale factor
#'@param pg pagerank values
#'@param vnames remove vertex labels
#'@importFrom tidyr %>%
#'@import colorBlindness
#'@return R default plot
#'@export
#'@examples
#'paths <- c('CTR' = system.file("extdata",
#'                               "ctr_nils_bm_human_newformat.csv",
#'                               package = "CrossTalkeR"),
#'           'EXP' = system.file("extdata",
#'                               "exp_nils_bm_human_newformat.csv",
#'                               package = "CrossTalkeR"))
#'
#'genes <- c('TGFB1')
#'
#'output =  system.file("extdata", package = "CrossTalkeR")
#'data <- generate_report(paths,
#'                        genes,
#'                        out_path=paste0(output,'/'),
#'                        threshold=0,
#'                        out_file = 'vignettes_example.html',
#'                        output_fmt = "html_document",
#'                        report = FALSE)
#'plot_cci(graph = data@graphs$CTR,
#'        colors = data@colors,
#'        plt_name = 'Example 1',
#'        coords = data@coords[igraph::V(data@graphs$CTR)$name,],
#'        emax = NULL,
#'        leg = FALSE,
#'        low = 0,
#'        high = 0,
#'        ignore_alpha = FALSE,
#'        log = FALSE,
#'        efactor = 8,
#'        vfactor = 12,
#'        vnames = TRUE)
plot_cci <- function(graph,
                     colors,
                     plt_name,
                     coords,
                     emax = NULL,
                     leg = FALSE,
                     low = 25,
                     high = 75,
                     ignore_alpha = FALSE,
                     log = FALSE,
                     efactor = 8,
                     vfactor = 12,
                     vnames = T,
                     pg = NULL) {

  # Check Maximal Weight
  if (is.null(emax)) {
    emax <- max(abs(igraph::E(graph)$weight))
  }
  # Using color pallet to up and down regulation
  col_pallet <- colorBlindness::Blue2DarkOrange18Steps
  # Expanding the pallet range
  col_pallet[10] <- '#B8b9ba'
  col_pallet <- grDevices::colorRampPalette(col_pallet)(201)
  # Checking looops
  edge_start <- igraph::ends(graph,
                             es = igraph::E(graph),
                             names = FALSE)
  # Scale nodes coordinates
  if (nrow(coords) != 1) {
    coords_scale <- scale(coords)
  }else {
    coords_scale <- coords
  }
  # It will make the loops in a correct angle
  loop_angle <- ifelse(coords_scale[igraph::V(graph)$name, 1] > 0,
                       -atan(coords_scale[igraph::V(graph)$name, 2] / coords_scale[igraph::V(graph)$name, 1]),
                       pi - atan(coords_scale[igraph::V(graph)$name, 2] / coords_scale[igraph::V(graph)$name, 1])
  )
  # Setting node colors
  igraph::V(graph)$color <- colors[igraph::V(graph)$name]
  ## Color scheme
  we <- round(oce::rescale(igraph::E(graph)$weight,
                           xlow = (-emax),
                           xhigh = emax,
                           rlow = 1,
                           rhigh = 200,
                           clip = TRUE),
              0)
  igraph::E(graph)$color <- col_pallet[we]
  alpha_cond <- (igraph::E(graph)$inter > low) & (igraph::E(graph)$inter < high)
  alpha <- ifelse(alpha_cond, 0, igraph::E(graph)$inter)
  subgraph <- igraph::delete.edges(graph,
                                   igraph::E(graph)[alpha == 0 | is.na(alpha)]
  )
  if (!ignore_alpha) {
    igraph::E(graph)$color <- scales::alpha(igraph::E(graph)$color, alpha)
  }
  ## Thickness and arrow size
  if (is.null(pg)) {
    igraph::V(graph)$size <- 60
  }
  else {
    igraph::V(graph)$size <- scales::rescale(pg, c(1, 60))
  }

  if (log) {
    igraph::E(graph)$width <- ifelse(igraph::E(graph)$inter != 0,
                                     log2(1 + igraph::E(graph)$inter),
                                     0) * efactor
  }else {
    igraph::E(graph)$width <- ifelse(igraph::E(graph)$inter != 0,
                                     igraph::E(graph)$inter,
                                     0) * efactor
  }
  igraph::E(graph)$arrow.size <- 0.4
  igraph::E(graph)$arrow.width <- igraph::E(graph)$width + 0.8
  igraph::E(graph)$loop.angle <- NA
  if (sum(edge_start[, 2] == edge_start[, 1]) != 0) {
    igraph::E(graph)$loop.angle[which(edge_start[, 2] == edge_start[, 1])] <- loop_angle[edge_start[which(edge_start[, 2] == edge_start[, 1]), 1]]
    igraph::E(graph)$loop.angle[which(edge_start[, 2] != edge_start[, 1])] <- 0
  }
  coords_scale[, 1] <- scales::rescale(coords_scale[, 1], from = c(-1, 1), to = c(-2, 2))
  coords_scale[, 2] <- scales::rescale(coords_scale[, 2], from = c(-1, 1), to = c(-2, 2))
  plot(graph,
       layout = coords_scale,
       xlim = c(-4, 4),
       ylim = c(-4, 4),
       rescale = F,
       edge.curved = 0.5,
       vertex.label = NA,
       vertex.shape = "circle",
       margin = 0.0,
       loop.angle = igraph::E(graph)$loop.angle,
       edge.label = NA,
       main = plt_name
  )
  # Thicknesse legend
  amin <- min(igraph::E(graph)$inter[igraph::E(graph)$inter != 0])
  amax <- max(igraph::E(graph)$inter)
  e_wid_sp <- c(amin,
                amin + amax / 2,
                amax)
  graphics::legend("topleft",
                   legend = round(e_wid_sp, 1),
                   col = "black",
                   title = "Percentage of the interactions",
                   pch = NA,
                   bty = "n",
                   cex = 1,
                   lwd = e_wid_sp,
                   lty = c(1, 1, 1),
                   horiz = FALSE)

  v <- igraph::V(graph)$size
  if (!is.null(pg)) {
    a <- graphics::legend('bottomleft',
                          title = "Node Pagerank",
                          legend = c("", "", ""),
                          pt.cex = c(min(v) + 1, mean(v), max(v)) / 12, col = 'black',
                          pch = 21, pt.bg = 'black', box.lwd = 0, y.intersp = 2)
    graphics::text(a$rect$left + a$rect$w, a$text$y,
                   c(round(min(pg), 2), round(mean(pg), 2), round(max(pg), 2)), pos = 2)
  }
  x <- coords_scale[, 1] * 1.2
  y <- coords_scale[, 2] * 1.2
  coord_ratio <- coords_scale[, 1] / coords_scale[, 2]
  angle <- ifelse(
    atan(-coord_ratio) * (180 / pi) < 0,
    90 + atan(-coord_ratio) * (180 / pi),
    270 + atan(-coord_ratio) * (180 / pi))
  if (vnames) {
    for (i in seq_len(length(x))) {
      graphics::text(x = x[i],
                     y = y[i],
                     labels = igraph::V(graph)$name[i],
                     adj = NULL,
                     pos = NULL,
                     cex = 0.8,
                     col = "black",
                     xpd = TRUE)
    }
  }
  if (leg) {
    # Edge Colormap
    if (min(igraph::E(graph)$weight) < 0) {
      netdiffuseR::drawColorKey(seq(1, 200),
                                tick.marks = c(1, 101, 200),
                                color.palette = col_pallet,
                                labels = c(-round(emax, 3), 0, round(emax, 3)),
                                nlevels = 200,
                                main = "Weights",
                                pos = 2,
                                key.pos = c(0.98, 1.0, 0.0, 0.2),
                                border = "transparent")
    }
    else {
      netdiffuseR::drawColorKey(seq(100, 200),
                                tick.marks = c(100, 200),
                                color.palette = col_pallet[100:201],
                                labels = c(0, round(emax, 3)),
                                nlevels = 100,
                                main = "Weights",
                                pos = 2,
                                key.pos = c(0.98, 1.0, 0.0, 0.2),
                                border = "transparent")
    }
  }
}


#'This function selected genes sankey plot
#'
#'@param lrobj_tbl LRobject table with all data
#'@param target gene
#'@param ligand_cluster Ligand Clusters
#'@param receptor_cluster Receptor Clusters
#'@param plt_name plot title
#'@param threshold top_n n value
#'@import ggplot2
#'@import stringr
#'@import dplyr
#'@import colorBlindness
#'@import ggalluvial
#'@importFrom tidyr %>%
#'@importFrom stats reorder
#'@return R default plot
#'@export
#'@examples
#'paths <- c('CTR' = system.file("extdata",
#'                               "ctr_nils_bm_human_newformat.csv",
#'                               package = "CrossTalkeR"),
#'           'EXP' = system.file("extdata",
#'                               "exp_nils_bm_human_newformat.csv",
#'                               package = "CrossTalkeR"))
#'output =  system.file("extdata", package = "CrossTalkeR")
#'genes <- c('TGFB1')
#'
#'data <- generate_report(paths,
#'                        genes,
#'                        out_path=paste0(output,'/'),
#'                        threshold=0,
#'                        out_file = 'vignettes_example.html',
#'                        output_fmt = "html_document",
#'                        report = FALSE)
plot_sankey <- function(lrobj_tbl,
                        target = NULL,
                        ligand_cluster = NULL,
                        receptor_cluster = NULL,
                        plt_name = NULL,
                        threshold = 50, tfflag = TRUE) {

  if (!is.null(target)) {
    if (length(stringr::str_split(target, "\\|")) > 1) {
      target_type = stringr::str_split(target, "\\|")[[1]][[2]]
      if (target_type == "R") {
        if (length(which(grepl("\\|", lrobj_tbl$Receptor))) > 0) {
          data <- lrobj_tbl %>%
            filter(Receptor == target)
        } else {
          target = stringr::str_split(target, "\\|")[[1]][[1]]
          data <- lrobj_tbl %>%
            filter(Receptor == target)
        }
      } else if (target_type == "L") {
        if (length(which(grepl("\\|", lrobj_tbl$Ligand))) > 0) {
          data <- lrobj_tbl %>%
            filter(Ligand == target)
        } else {
          target = stringr::str_split(target, "\\|")[[1]][[1]]
          data <- lrobj_tbl %>%
            filter(Ligand == target)
        }
      }
    } else {
      data <- lrobj_tbl[grepl(target, lrobj_tbl$allpair),]
    }
  }
  else {
    data <- lrobj_tbl
  }
  if (!is.null(ligand_cluster)) {
    tmp_sel <- grepl(ligand_cluster, data$ligpair)
    data <- data[tmp_sel,]
  }
  if (!is.null(receptor_cluster)) {
    tmp_sel <- grepl(receptor_cluster, data$recpair)
    data <- data[tmp_sel,]
  }
  colp <- c(Blue2DarkOrange18Steps[4], Blue2DarkOrange18Steps[14])
  tmp_cols <- c("source", "Ligand", "Receptor", "target")
  names(colp) <- c("FALSE", "TRUE")
  if (dim(data)[1] >= 1) {
    data$freq <- 1
    tmp <- dplyr::top_n(data, ifelse(dim(data)[1] > threshold, threshold,
                                     dim(data)[1]), abs(.data$LRScore))
    print(ggplot2::ggplot(tmp, aes(y = .data$freq, axis1 = .data$Ligand.Cluster,
                                   axis2 = stats::reorder(.data$Ligand, .data$LRScore),
                                   axis3 = stats::reorder(.data$Receptor, .data$LRScore),
                                   axis4 = .data$Receptor.Cluster)) +
            ggalluvial::geom_alluvium(aes(fill = .data$LRScore, color = 'b'),
                                      width = 1 / 12,
                                      discern = FALSE) +
            ggalluvial::geom_stratum(width = 1 / 12) +
            ggplot2::geom_label(stat = ggalluvial::StatStratum,
                                ggplot2::aes(label = ggplot2::after_stat(.data$stratum)),
                                size = 4) +
            ggplot2::scale_x_discrete(limits = tmp_cols, expand = c(.05, .05)) +
            ggplot2::scale_fill_gradient2(low = colorBlindness::Blue2DarkOrange18Steps[4],
                                          mid = colorBlindness::Blue2DarkOrange18Steps[10],
                                          high = colorBlindness::Blue2DarkOrange18Steps[14], midpoint = 0) +
            ggplot2::scale_color_manual(values = c("black")) +
            ggplot2::ggtitle(plt_name) +
            ggplot2::theme(text = element_text(size = 8)) +
            ggplot2::theme_minimal()
    )
  }
  else {
    print(paste0("Gene->", target, "Not Found"))
  }
}

#'This function assigns coordinates in the sankey plots for the genes on the dataframe
#'
#'@param df dataframe with genes contained in the sankey plot
#'@param type string defining the type of the genes inside the dataframe (L, R or TF)
set_coords <- function(df, type) {
  if (length(df$gene) == 1) {
    coords = c(0)
  } else if ((length(df$gene) %% 2) == 0) {
    x = length(df$gene) / 2 * 5
    coords = seq(-x + 2.5, x - 2.5, by = 5)
  } else {
    x = (floor(length(df$gene) / 2)) * 5
    coords = seq(-x, x, by = 5)
  }
  df$y = rev(coords)
  if (type == "R") {
    df$x = 5
  } else if (type == "TF") {
    df$x = 15
  } else {
    df$x = 25
  }
  return(df)
}


#'This function creates a sankey plot for a selected gene including transcription factor interactions.
#'
#'@param lrobj_tbl LRobject table with all data
#'@param target gene
#'@param cluster cluster
#'@param target_type type of target
#'@param plt_name plot title
#'@param threshold top_n n value
#'@import ggplot2
#'@import dplyr
#'@import colorBlindness
#'@import ggalluvial
#'@importFrom tidyr %>%
#'@importFrom stats reorder
#'@return R default plot
#'@export
#'@examples
#'paths <- c('CTR' = system.file("extdata",
#'                               "ctr_nils_bm_human.csv",
#'                               package = "CrossTalkeR"),
#'           'EXP' = system.file("extdata",
#'                               "exp_nils_bm_human_newformat.csv",
#'                               package = "CrossTalkeR"))
#'output =  system.file("extdata", package = "CrossTalkeR")
#'genes <- c('TGFB1')
#'
#'data <- generate_report(paths,
#'                        genes,
#'                        out_path=paste0(output,'/'),
#'                        threshold=0,
#'                        out_file = 'vignettes_example.html',
#'                        output_fmt = "html_document",
#'                        report = FALSE)
plot_graph_sankey_tf <- function(lrobj_tbl,
                                 pagerank_table,
                                 target = NULL,
                                 cluster = NULL,
                                 target_type = NULL,
                                 plt_name = NULL,
                                 threshold = 50,
                                 save_path = NULL) {


  if (!is.null(target)) {
    if (target_type == "TF") {

      data = lrobj_tbl %>%
        filter(type_gene_B == "Transcription Factor" | type_gene_A == "Transcription Factor") %>%
        filter(Receptor == target | Ligand == target)

    } else if (target_type == "R") {

      receptor_interactions = lrobj_tbl %>%
        filter(type_gene_B == "Transcription Factor") %>%
        filter(Ligand == target)
      ligand_interactions = lrobj_tbl %>%
        filter(type_gene_A == "Transcription Factor") %>%
        filter(Ligand %in% receptor_interactions$Receptor)

      data = rbind(receptor_interactions, ligand_interactions)

    } else {

      ligand_interactions = lrobj_tbl %>%
        filter(type_gene_A == "Transcription Factor") %>%
        filter(Receptor == target)
      receptor_interactions = lrobj_tbl %>%
        filter(type_gene_B == "Transcription Factor") %>%
        filter(Receptor %in% ligand_interactions$Ligand)

      data = rbind(receptor_interactions, ligand_interactions)

    }

    data = data %>%
      filter(Ligand.Cluster == cluster) %>%
      subset(
        select = c(
          "Ligand",
          "Receptor",
          "type_gene_A",
          "type_gene_B",
          "Ligand.Cluster",
          "ligpair",
          "recpair",
          "LRScore"
        )
      ) %>%
      filter(!(type_gene_B == "Transcription Factor" & type_gene_A == "Transcription Factor"))

    if (dim(data)[1] > 0) {
      data_RTF = data %>%
        filter(type_gene_B == "Transcription Factor")
      data_RTF$Pagerank_Score <- pagerank_table$Pagerank[match(data_RTF$ligpair, pagerank_table$nodes)]
      data_RTF$TF_Pagerank_Score <- pagerank_table$Pagerank[match(data_RTF$recpair, pagerank_table$nodes)]

      data_TFL = data %>%
        filter(type_gene_A == "Transcription Factor")
      data_TFL$Pagerank_Score <- pagerank_table$Pagerank[match(data_TFL$recpair, pagerank_table$nodes)]
      data_TFL$TF_Pagerank_Score <- pagerank_table$Pagerank[match(data_TFL$ligpair, pagerank_table$nodes)]

      data = rbind(data_RTF, data_TFL)

      data <- as.data.frame(lapply(data, unlist))

      data_group1 <- subset(data, type_gene_A %in% c("Receptor"))
      data_group2 <- subset(data, type_gene_A %in% c("Transcription Factor"))

      if (target_type == "TF") {
        gene_list1 <- data_group1 %>%
          arrange(desc(abs(Pagerank_Score))) %>%
          head(10)

        gene_list2 <- data_group2 %>%
          arrange(desc(abs(Pagerank_Score))) %>%
          head(10)
      } else if (target_type == "R") {
        gene_list1 <- data_group1 %>%
          arrange(desc(abs(TF_Pagerank_Score))) %>%
          head(10)

        gene_list2 <- data_group2 %>%
          arrange(desc(abs(Pagerank_Score))) %>%
          head(10)
      } else {
        gene_list1 <- data_group1 %>%
          arrange(desc(abs(Pagerank_Score))) %>%
          head(10)

        gene_list2 <- data_group2 %>%
          arrange(desc(abs(TF_Pagerank_Score))) %>%
          head(10)
      }

      graph_df = rbind(gene_list1, gene_list2)

      graph1 <- igraph::graph_from_data_frame(graph_df[, c("Ligand", "Receptor", "Pagerank_Score")])

      receptors_coord = gene_list1 %>%
        select(Ligand, Pagerank_Score) %>%
        rename(gene = Ligand, score = Pagerank_Score) %>%
        arrange(desc(score)) %>%
        unique()

      ligands_coord = gene_list2 %>%
        select(Receptor, Pagerank_Score) %>%
        rename(gene = Receptor, score = Pagerank_Score) %>%
        arrange(desc(score)) %>%
        unique()

      tf_coord_r = gene_list1 %>%
        select(Receptor, TF_Pagerank_Score) %>%
        rename(gene = Receptor, score = TF_Pagerank_Score) %>%
        unique()
      tf_coord_l = gene_list2 %>%
        select(Ligand, TF_Pagerank_Score) %>%
        rename(gene = Ligand, score = TF_Pagerank_Score) %>%
        unique()
      tf_coord = rbind(tf_coord_r, tf_coord_l) %>%
        arrange(desc(score)) %>%
        unique()

      if (dim(ligands_coord)[1] > 0) {
        res_df = set_coords(ligands_coord, "L")
      }

      if (dim(tf_coord)[1] > 0) {
        if (exists('res_df')) {
          res_df = rbind(res_df, set_coords(tf_coord, "TF"))
        } else {
          res_df = set_coords(tf_coord, "TF")
        }
      }

      if (dim(receptors_coord)[1] > 0) {
        if (exists('res_df')) {
          res_df = rbind(res_df, set_coords(receptors_coord, "R"))
        } else {
          res_df = set_coords(receptors_coord, "R")
        }
      }

      rownames(res_df) = res_df$gene

      for (vertice in V(graph1)) {
        name = vertex_attr(graph1, "name", vertice)
        if (grepl("|R", name, fixed = TRUE)) {
          vertex_attr(graph = graph1, name = "Gene_Type", index = vertice) <- "Receptor"
          vertex_attr(graph = graph1, name = "Score", index = vertice) <- pagerank_table[paste0(cluster, name),]
          vertex_attr(graph = graph1, name = "clustername", index = vertice) <- paste0(cluster, "/", name)
        } else if (grepl("|L", name, fixed = TRUE)) {
          vertex_attr(graph = graph1, name = "Gene_Type", index = vertice) <- "Ligand"
          vertex_attr(graph = graph1, name = "Score", index = vertice) <- pagerank_table[paste0(cluster, name),]
          vertex_attr(graph = graph1, name = "clustername", index = vertice) <- paste0(cluster, "/", name)
        } else {
          vertex_attr(graph = graph1, name = "Gene_Type", index = vertice) <- "Transcription Factor"
          vertex_attr(graph = graph1, name = "Score", index = vertice) <- pagerank_table[paste0(cluster, name),]
          vertex_attr(graph = graph1, name = "clustername", index = vertice) <- paste0(cluster, "/", name)
        }
      }

      pagerank_list = setNames(as.list(pagerank_table$Pagerank), pagerank_table$nodes)
      test_result_pg = pagerank_list[V(graph1)$clustername]
      igraph::V(graph1)$size <- scales::rescale(unlist(test_result_pg), c(1, 60))

      mat_coords <- matrix(, nrow = length(V(graph1)), ncol = 2)

      for (vertice in V(graph1)) {
        name = vertex_attr(graph1, "name", vertice)
        mat_coords[vertice,] = c(res_df[name,]$x, res_df[name,]$y)
      }


      print(ggraph(graph1, layout = mat_coords) +
              geom_edge_fan(aes(colour = Pagerank_Score), width = 2, arrow = arrow(angle = 30, length = unit(4, "mm"))) +
              geom_node_point(aes(size = size)) +
              geom_node_point(aes(color = Gene_Type, size = size), show.legend = FALSE) +
              scale_size_area(name = "Node Pageank") +
              geom_node_text(aes(label = name), size = 4, nudge_y = 1) +
              scale_edge_color_gradient2(low = colorBlindness::Blue2DarkOrange18Steps[4],
                                         mid = "white",
                                         high = colorBlindness::Blue2DarkOrange18Steps[14], midpoint = 0, name = "Pagerank Score") +
              coord_cartesian(clip = "off") +
              theme_void() +
              annotate(geom = "text", x = 5, y = max(res_df$y) + 5, label = "Receptor", size = 5, fontface = "bold") +
              annotate(geom = "text", x = 15, y = max(res_df$y) + 5, label = "Transcription Factor", size = 5, fontface = "bold") +
              annotate(geom = "text", x = 25, y = max(res_df$y) + 5, label = "Ligand", size = 5, fontface = "bold") +
              theme(plot.margin = unit(rep(30, 4), "points")) +
              ggtitle(plt_name)
      )

    } else {
      print(paste0("Target gene ", target, " not found in selected cluster ", cluster, "!"))
    }

  }else {
    print("Please provide an target gene to filter the interactions!")
  }
}


#'This function is a proxy to the PCA plot
#'
#'@param lrobj_tbl LRobject table with all data
#'@param curr table entry
#'@param dims PCA dims
#'@param ret return plot
#'@param ggi GGI mode 
#'@import ggplot2
#'@import ggrepel
#'@import factoextra
#'@importFrom tidyr %>%
#'@importFrom stats reorder
#'@return R default plot
#'@export
#'@examples
#'paths <- c('CTR' = system.file("extdata",
#'                               "ctr_nils_bm_human_newformat.csv",
#'                               package = "CrossTalkeR"),
#'           'EXP' = system.file("extdata",
#'                               "exp_nils_bm_human_newformat.csv",
#'                               package = "CrossTalkeR"))
#'output =  system.file("extdata", package = "CrossTalkeR")
#'genes <- c('TGFB1')
#'
#'data <- generate_report(paths,
#'                        genes,
#'                        out_path=paste0(output,'/'),
#'                        threshold=0,
#'                        out_file = 'vignettes_example.html',
#'                        output_fmt = "html_document",
#'                        report = FALSE)
plot_pca <- function(lrobj_tblPCA, curr, dims = c(1, 2), ret = F, ggi = TRUE) { ## Pweeseee Test me!!!
  x <- max(abs(lrobj_tblPCA[[curr]]$x[, dims[1]]))
  y <- max(abs(lrobj_tblPCA[[curr]]$x[, dims[2]]))
  if (ggi) {
    z_x <- lrobj_tblPCA[[curr]]$x[, dims[1]]
    z_y <- lrobj_tblPCA[[curr]]$x[, dims[2]]
    ver_zx <- ifelse(abs(z_x) > 2 * lrobj_tblPCA[[curr]]$sdev[1], 1, 0)
    ver_zy <- ifelse(abs(z_y) > 2 * lrobj_tblPCA[[curr]]$sdev[2], 1, 0)
    plt1 <- factoextra::fviz_pca_biplot(lrobj_tblPCA[[curr]],
                                        axes = dims,
                                        pointshape = 21, pointsize = 0.5, labelsize = 3,
                                        repel = TRUE, max.overlaps = 100, label = 'var') +
      ggrepel::geom_label_repel(aes(label = ifelse(ver_zx & ver_zy, rownames(lrobj_tblPCA[[curr]]$x), NA)),
                                hjust = 0, vjust = 0, size = 3, max.overlaps = 100) +
      xlim(-x, x) +
      ylim(-y, y) +
      ggtitle(curr) +
      theme(text = element_text(size = 3),
            axis.title = element_text(size = 7.5),
            axis.text = element_text(size = 12))

  }else {
    plt1 <- factoextra::fviz_pca_biplot(lrobj_tblPCA[[curr]],
                                        axes = dims,
                                        pointshape = 21, pointsize = 0.5, labelsize = 6,
                                        repel = TRUE, max.overlaps = 100, label = 'var') +
      ggrepel::geom_text_repel(aes(label = rownames(lrobj_tblPCA[[curr]]$x))) +
      xlim(-x, x) +
      ylim(-y, y) +
      ggtitle(curr) +
      theme(text = element_text(size = 7.5),
            axis.title = element_text(size = 7.5),
            axis.text = element_text(size = 12))

  }
  if (ret) {
    return(plt1)
  }
}

#'This function is a proxy to the PCA plot in comparative conditions
#'
#'@param lrobj_tblPCA LRobject table with all data
#'@param curr table entry
#'@param dims PCA dims
#'@param ret return plot
#'@param ggi GGI mode
#'@param include_tf intracellular option
#'@param gene_types filter option of genes
#'@import ggplot2
#'@import ggrepel
#'@import factoextra
#'@importFrom tidyr %>%
#'@importFrom stats reorder
#'@return R default plot
#'@export
#'@examples
#'paths <- c('CTR' = system.file("extdata",
#'                               "ctr_nils_bm_human_newformat.csv",
#'                               package = "CrossTalkeR"),
#'           'EXP' = system.file("extdata",
#'                               "exp_nils_bm_human_newformat.csv",
#'                               package = "CrossTalkeR"))
#'output =  system.file("extdata", package = "CrossTalkeR")
#'genes <- c('TGFB1')
#'
#'data <- generate_report(paths,
#'                        genes,
#'                        out_path=paste0(output,'/'),
#'                        threshold=0,
#'                        out_file = 'vignettes_example.html',
#'                        output_fmt = "html_document",
#'                        report = FALSE)
plot_pca_LR_comparative <- function(lrobj_tblPCA, dims = c(1, 2), ret = F, ggi = TRUE, include_tf = TRUE, gene_types = "all") {

  if (ggi) {
    pca_table = names(lrobj_tblPCA@pca)[which(
      grepl("_x_", names(lrobj_tblPCA@pca))
        & grepl("_ggi", names(lrobj_tblPCA@pca)))]

    #Filter for LR or TF
    if (gene_types == "LR") {
      pca_split = lrobj_tblPCA@pca[[pca_table]]$x[, 1]
      pca_split_names = names(pca_split)
      result_split_names = pca_split_names[grepl("|R", pca_split_names, fixed = TRUE) | grepl("|L", pca_split_names, fixed = TRUE)]
      col.num <- which(rownames(lrobj_tblPCA@pca[[pca_table]]$x) %in% result_split_names)
      lrobj_tblPCA@pca[[pca_table]]$x <- lrobj_tblPCA@pca[[pca_table]]$x[sort(c(col.num)),]

    } else if (gene_types == "TF") {
      pca_split = lrobj_tblPCA@pca[[pca_table]]$x[, 1]
      pca_split_names = names(pca_split)
      result_split_names = pca_split_names[grepl("|TF", pca_split_names, fixed = TRUE)]
      col.num <- which(rownames(lrobj_tblPCA@pca[[pca_table]]$x) %in% result_split_names)
      lrobj_tblPCA@pca[[pca_table]]$x <- lrobj_tblPCA@pca[[pca_table]]$x[sort(c(col.num)),]
    }

    #Mapping_Table
    if (include_tf) {
      map_df <- as.data.frame(rownames(lrobj_tblPCA@pca[[pca_table]]$x))
      colnames(map_df) <- c("gene")
      map_df$mapping <- sapply(map_df$gene, function(gene) {
        if (grepl("|R", gene, fixed = TRUE)) {
          txt <- "Receptor"
        } else if (grepl("|L", gene, fixed = TRUE)) {
          txt <- "Ligand"
        } else {
          txt <- "Transcription Factor"
        }
        return(txt)
      })
      color_groups = c("#f8756b", "#00b835", "#619cff")
    }else {
      l_mapping = lrobj_tblPCA@tables[[names(lrobj_tblPCA@tables)[which(
        grepl("_x_", names(lrobj_tblPCA@tables)))]]] %>%
        select(ligpair, type_gene_A) %>%
        rename(gene = ligpair, mapping = type_gene_A) %>%
        distinct()

      r_mapping = lrobj_tblPCA@tables[[names(lrobj_tblPCA@tables)[which(
        grepl("_x_", names(lrobj_tblPCA@tables)))]]] %>%
        select(recpair, type_gene_B) %>%
        rename(gene = recpair, mapping = type_gene_B) %>%
        distinct()

      map_df = rbind(l_mapping, r_mapping)
      map_df = dplyr::filter(map_df, gene %in% rownames(lrobj_tblPCA@pca[[pca_table]]$x))
      map_df = map_df[!duplicated(map_df$gene),]
      color_groups = c("#f8756b", "#00b835")
    }

    rmd_title <- paste0(pca_table, '_tbl')
    rmd_title1 <- paste0(pca_table, '_pca')
    x <- max(abs(lrobj_tblPCA@pca[[pca_table]]$x[, dims[[1]]]))
    y <- max(abs(lrobj_tblPCA@pca[[pca_table]]$x[, dims[[2]]]))
    z_x <- lrobj_tblPCA@pca[[pca_table]]$x[, dims[[1]]]
    z_y <- lrobj_tblPCA@pca[[pca_table]]$x[, dims[[2]]]
    ver_zx <- ifelse(abs(z_x) >= (2 * lrobj_tblPCA@pca[[pca_table]]$sdev[dims[[1]]]), 1, 0)
    ver_zy <- ifelse(abs(z_y) >= (2 * lrobj_tblPCA@pca[[pca_table]]$sdev[dims[[2]]]), 1, 0)
    pca_plot <- fviz_pca_biplot(lrobj_tblPCA@pca[[pca_table]],
                                axes = c(dims[[1]], dims[[2]]),
                                pointshape = 20, pointsize = 2, labelsize = 10,
                                repel = FALSE, max.overlaps = 100, label = 'var', habillage = map_df$mapping, palette = color_groups) +
      geom_label_repel(aes(label = ifelse((ver_zx | ver_zy), rownames(lrobj_tblPCA@pca[[pca_table]]$x), NA)), size = 5) +
      xlim(-x, x) +
      ylim(-y, y) +
      ggtitle(pca_table) +
      theme(text = element_text(size = 7.5),
            axis.title = element_text(size = 7.5),
            axis.text = element_text(size = 7.5))
  } else {
    pca_table = names(lrobj_tblPCA@pca)[which(
      grepl("_x_", names(lrobj_tblPCA@pca))
        & !grepl("_ggi", names(lrobj_tblPCA@pca)))]

    rmd_title <- paste0(pca_table, '_tbl')
    rmd_title1 <- paste0(pca_table, '_pca')
    x <- max(abs(lrobj_tblPCA@pca[[pca_table]]$x[, dims[[1]]]))
    y <- max(abs(lrobj_tblPCA@pca[[pca_table]]$x[, dims[[2]]]))
    pca_plot <- fviz_pca_biplot(lrobj_tblPCA@pca[[pca_table]],
                                axes = c(1, 2),
                                pointshape = 21, pointsize = 0.5, labelsize = 6,
                                repel = TRUE, max.overlaps = 100, label = 'var') +
      geom_text_repel(aes(label = rownames(lrobj_tblPCA@pca[[pca_table]]$x))) +
      xlim(-x, x) +
      ylim(-y, y) +
      ggtitle(pca_table) +
      theme(text = element_text(size = 7.5),
            axis.title = element_text(size = 7.5),
            axis.text = element_text(size = 7.5))
  }


  if (ret) {
    return(pca_plot)
  }
}

