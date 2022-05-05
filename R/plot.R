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
#'                               "ctr_nils_bm_human.csv",
#'                               package = "CrossTalkeR"),
#'           'EXP' = system.file("extdata",
#'                               "exp_nils_bm_human.csv",
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


#'This function do a ggi plot and higest degree nodes
#'
#'@param graph graph
#'@param color cluster color
#'@param name plot header
#'@import ggplot2
#'@import ggraph
#'@import igraph
#'@import colorBlindness
#'@import graphlayouts
#'@return R default plot
#'@importFrom tidyr %>%
#'@export
#'@examples
#'paths <- c('CTR' = system.file("extdata",
#'                               "ctr_nils_bm_human.csv",
#'                               package = "CrossTalkeR"),
#'           'EXP' = system.file("extdata",
#'                               "exp_nils_bm_human.csv",
#'                               package = "CrossTalkeR"))
#'output =  system.file("extdata", package = "CrossTalkeR")
#'genes <- c('TGFB1')
#'data <- generate_report(paths,
#'                        genes,
#'                        out_path=paste0(output,'/'),
#'                        threshold=0,
#'                        out_file = 'vignettes_example.html',
#'                        output_fmt = "html_document",
#'                        report = FALSE)
#'plot_ggi(graph = data@graphs_ggi$EXP_x_CTR,
#'         color = data@colors,name="EXP_x_CTR")
plot_ggi <- function(graph, color, name) {
  deg <- igraph::degree(graph)
  v_names <- igraph::V(graph)$name
  v_names <- tibble::as_tibble(v_names)
  v_names <- v_names %>%
    tidyr::separate(.data$value, c("genes", "cluster"), "/")
  igraph::V(graph)$genes <- v_names$genes
  igraph::V(graph)$cluster <- v_names$cluster
  lr2col <- round(oce::rescale(igraph::E(graph)$LRScore,
                               xlow = min(igraph::E(graph)$LRScore),
                               rlow = 1,
                               rhigh = 25
  )
  )
  ecolors <- rev(grDevices::colorRampPalette(colorBlindness::Blue2DarkOrange18Steps)(25))[lr2col]
  igraph::E(graph)$colors <- ecolors
  cls <- names(color)
  names(color) <- NULL
  ewidth <- abs(igraph::E(graph)$LRScore) - mean(abs(igraph::E(graph)$LRScore))
  ewidth <- ewidth / stats::sd(abs(igraph::E(graph)$LRScore))
  ewidth <- (ewidth - min(ewidth)) / (max(ewidth) - min(ewidth))
  print(ggraph::ggraph(graph, layout = "stress") +
          ggraph::geom_edge_link0(ggplot2::aes(edge_width = ewidth,
                                               alpha = ewidth)) +
          ggraph::geom_node_point(size = ((deg / max(deg)) * 10),
                                  alpha = 1,
                                  ggplot2::aes(color = igraph::V(graph)$cluster)
          ) +
          ggplot2::scale_colour_manual(values = color,
                                       labels = cls,
                                       name = "Clusters") +
          ggraph::geom_node_label(ggplot2::aes(filter = deg > deg[order(deg, decreasing = TRUE)][ifelse(length(deg) > 100, 100, length(deg))] & igraph::V(graph) %in% igraph::articulation.points(graph),
                                               label = igraph::V(graph)$genes,
                                               color = igraph::V(graph)$cluster
          ),
                                  repel = TRUE,
                                  hjust = "inward",
                                  size = 7,
                                  show.legend = FALSE) +
          ggraph::scale_edge_width_continuous(range = c(0, 1)) +
          ggplot2::ggtitle(name) +
          #scale_size(range = c(1,6))+
          ggraph::theme_graph(base_family = "sans") +
          ggplot2::theme(legend.position = "left"))

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
#'                               "exp_nils_bm_human.csv",
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
                        threshold = 50) {

  target_type = str_split(target, "\\|")[[1]][[2]]

  if (!is.null(target)) {
    #data <- lrobj_tbl[grepl(target, lrobj_tbl$allpair), ]

    if (target_type == "R") {
      data <- lrobj_tbl %>%
        filter(Receptor == target)
    } else if (target_type == "L") {
      data <- lrobj_tbl %>%
        filter(Ligand == target)
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
            ggalluvial::geom_alluvium(aes(fill = .data$LRScore > 0, color = 'b'),
                                      width = 1 / 12,
                                      discern = FALSE) +
            ggalluvial::geom_stratum(width = 1 / 12) +
            ggplot2::geom_label(stat = ggalluvial::StatStratum,
                                ggplot2::aes(label = ggplot2::after_stat(.data$stratum)),
                                size = 4) +
            ggplot2::scale_x_discrete(limits = tmp_cols, expand = c(.05, .05)) +
            ggplot2::scale_fill_manual(values = colp,
                                       limits = names(colp),
                                       name = "Upregulated") +
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


#'This function selected gene transcription factor interaction related sankey plot
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
#'                               "exp_nils_bm_human.csv",
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
plot_sankey_tf <- function(lrobj_tbl,
                           target = NULL,
                           cluster = NULL,
                           target_type = NULL,
                           plt_name = NULL,
                           threshold = 50) {

  lrobj_tbl = lrobj_tbl %>%
    filter(Cluster.x == cluster)

  if (!is.null(target)) {
    if (target_type == "R") {
      data <- lrobj_tbl %>%
        filter(Receptor == target)
    } else if (target_type == "L") {
      data <- lrobj_tbl %>%
        filter(Ligand == target)
    } else {
      data <- lrobj_tbl %>%
        filter(TF == target)
    }
  }
  else {
    data <- lrobj_tbl
  }

  plt_name = paste0(cluster, " Sankey ", target)
  colp <- c(Blue2DarkOrange18Steps[4], Blue2DarkOrange18Steps[14])
  tmp_cols <- c("Receptor", "TF", "Ligand")
  names(colp) <- c("FALSE", "TRUE")
  if (dim(data)[1] >= 1) {
    data$freq <- 1
    tmp <- dplyr::top_n(data, ifelse(dim(data)[1] > threshold, threshold,
                                     dim(data)[1]), abs(.data$TFScore.x))
    print(ggplot2::ggplot(tmp, aes(y = .data$freq, axis1 = .data$Receptor,
                                   axis2 = stats::reorder(.data$TF, .data$TFScore.x),
                                   axis3 = stats::reorder(.data$Ligand, .data$TFScore.x))) +
            ggalluvial::geom_alluvium(aes(fill = .data$TFScore.x > 0),
                                      width = 1 / 12,
                                      discern = FALSE) +
            ggalluvial::geom_stratum(width = 1 / 12) +
            ggplot2::geom_label(stat = ggalluvial::StatStratum,
                                ggplot2::aes(label = ggplot2::after_stat(.data$stratum)),
                                size = 4) +
            ggplot2::scale_x_discrete(limits = tmp_cols, expand = c(.05, .05)) +
            ggplot2::scale_fill_manual(values = colp,
                                       limits = names(colp),
                                       name = "Upregulated") +
            ggplot2::ggtitle(plt_name) +
            ggplot2::theme(text = element_text(size = 8)) +
            ggplot2::theme_minimal()
    )
  }
  else {
    print(paste0("Gene->", target, "Not Found"))
  }
}


#'This function selected gene transcription factor interaction related sankey plot
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
#'                               "exp_nils_bm_human.csv",
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
plot_sankey_tf_pagerank <- function(lrobj_tbl,
                                    pagerank_table,
                                    target = NULL,
                                    cluster = NULL,
                                    target_type = NULL,
                                    plt_name = NULL,
                                    threshold = 50,
                                    save_path = NULL) {

  lrobj_tbl = lrobj_tbl %>%
    filter(Cluster.x == cluster)

  if (!is.null(target)) {
    if (target_type == "R") {
      data <- lrobj_tbl %>%
        filter(Receptor == target)
    } else if (target_type == "L") {
      data <- lrobj_tbl %>%
        filter(Ligand == target)
    } else {
      data <- lrobj_tbl %>%
        filter(TF == target)
    }
  }
  else {
    data <- lrobj_tbl
  }

  data$combined_gene = paste(data$Cluster.x, data$Receptor, sep = "/")
  receptor_list = data$combined_gene

  receptor_list_unique = unique(receptor_list)

  pagerank_table <- pagerank_table %>%
    filter(nodes %in% receptor_list_unique) %>%
    arrange(desc(abs(Pagerank)))

  sankey_table = data %>%
    merge(pagerank_table, by.x = "combined_gene", by.y = "nodes") %>%
    arrange(desc(abs(Pagerank)))

  if (!is.null(save_path)) {
    write.csv(sankey_table, paste0(dirname(save_path), "/", target, "_", cluster, "_sankey_table.csv"))
  }

  filtered_pagerank <- pagerank_table %>%
    head(15)

  filtered_data = data %>%
    filter(combined_gene %in% filtered_pagerank$nodes)

  data = filtered_data

  plt_name = paste0(cluster, " Sankey ", target)
  colp <- c(Blue2DarkOrange18Steps[4], Blue2DarkOrange18Steps[14])
  tmp_cols <- c("Receptor", "TF", "Ligand")
  names(colp) <- c("FALSE", "TRUE")
  if (dim(data)[1] >= 1) {
    data$freq <- 1
    # tmp <- dplyr::top_n(data, ifelse(dim(data)[1] > threshold, threshold,
    #                                  dim(data)[1]), abs(.data$TFScore.x))
    tmp <- data
    print(ggplot2::ggplot(tmp, aes(y = .data$freq, axis1 = .data$Receptor,
                                   axis2 = stats::reorder(.data$TF, .data$TFScore.x),
                                   axis3 = stats::reorder(.data$Ligand, .data$TFScore.x))) +
            ggalluvial::geom_alluvium(aes(fill = .data$TFScore.x > 0),
                                      width = 1 / 12,
                                      discern = FALSE) +
            ggalluvial::geom_stratum(width = 1 / 12) +
            ggplot2::geom_label(stat = ggalluvial::StatStratum,
                                ggplot2::aes(label = ggplot2::after_stat(.data$stratum)),
                                size = 4) +
            ggplot2::scale_x_discrete(limits = tmp_cols, expand = c(.05, .05)) +
            ggplot2::scale_fill_manual(values = colp,
                                       limits = names(colp),
                                       name = "Upregulated") +
            ggplot2::ggtitle(plt_name) +
            ggplot2::theme(text = element_text(size = 8)) +
            ggplot2::theme_minimal()
    )
  }
  else {
    print(paste0("Gene->", target, "Not Found"))
  }
}


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
  df$y = coords
  if (type == "R") {
    df$x = 5
  } else if (type == "TF") {
    df$x = 15
  } else {
    df$x = 25
  }
  return(df)
}


#'This function selected gene transcription factor interaction related sankey plot
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
#'                               "exp_nils_bm_human.csv",
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


      recptors_coord = gene_list1 %>%
        select(Ligand) %>%
        rename(gene = Ligand) %>%
        unique()

      ligands_coord = gene_list2 %>%
        select(Receptor) %>%
        rename(gene = Receptor) %>%
        unique()

      tf_coord_r = gene_list1 %>%
        select(Receptor) %>%
        rename(gene = Receptor) %>%
        unique()
      tf_coord_l = gene_list2 %>%
        select(Ligand) %>%
        rename(gene = Ligand) %>%
        unique()
      tf_coord = rbind(tf_coord_r, tf_coord_l) %>%
        unique()

      if (dim(ligands_coord)[1] > 0) {
        res_df = set_coords(ligands_coord, "L")
      }

      if (dim(tf_coord) > 0) {
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
                                         mid = colorBlindness::Blue2DarkOrange18Steps[10],
                                         high = colorBlindness::Blue2DarkOrange18Steps[14], midpoint = 0, name = "Pagerank Score") +
              coord_cartesian(clip = "off") +
              theme_void() +
              annotate(geom = "text", x = 5, y = max(res_df$y) + 5, label = "Receptor", size = 5, fontface = "bold") +
              annotate(geom = "text", x = 15, y = max(res_df$y) + 5, label = "Transcription Factor", size = 5, fontface = "bold") +
              annotate(geom = "text", x = 25, y = max(res_df$y) + 5, label = "Ligand", size = 5, fontface = "bold") +
              theme(plot.margin = unit(rep(30, 4), "points"))
      )

    } else {
      print(paste0("Target gene ", target, " not found in selected cluster ", cluster, "!"))
    }

  }else {
    print("Please provide an target gene to filter the interactions!")
  }
}


#'This function selected genes correlation
#'
#'@param data Seurat object
#'@param lrobject LRObj object
#'@param pair pair of interest
#'@param lrslot table from lrobject
#'@param assay Seurat assay
#'@import Seurat ggsci ggnewscale tibble ggridges gtools Rmagic
#'@importFrom tidyr %>%
#'@importFrom stats reorder
#'@return R default plot
#'@export
plotInterInfo <- function(data, lrobject, pair, lrslot, assay = 'RNA') {
  Seurat::DefaultAssay(data) <- assay
  pair <- lrobject@tables[[lrslot]][grepl(pair, datalr@tables$EXP_x_CTR$allpair),]
  sender <- pair$Ligand.Cluster
  receiver <- pair$Receptor.Cluster
  exp <- Seurat::GetAssayData(data)
  lpair <- tibble::tibble(cells = names(exp[pair$Ligand, grepl(sender, Idents(data))]),
                          expr = exp[pair$Ligand, grepl(sender, Idents(data))],
                          cellpop = pair$ligpair,
                          stage = data@meta.data[names(exp[pair$Ligand, grepl(sender, Idents(data))]),]$protocol)
  rpair <- tibble::tibble(cells = names(exp[pair$Receptor, grepl(receiver, Idents(data))]),
                          expr = exp[pair$Receptor, grepl(receiver, Idents(data))],
                          cellpop = pair$recpair,
                          stage = data@meta.data[names(exp[pair$Receptor, grepl(receiver, Idents(data))]),]$protocol)
  pairt <- dplyr::bind_rows(lpair, rpair)
  p1 <- ggplot2::ggplot(pairt, aes(x = expr, y = cellpop, fill = stage)) +
    ggridges::geom_density_ridges() +
    ggplot2::facet_grid(. ~ stage) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = 'Expression density plot')
  tmp <- Seurat::Embeddings(data, 'umap')
  tmp <- tibble::tibble(c1 = tmp[, "UMAP_1"], c2 = tmp[, "UMAP_2"], cells = rownames(tmp))
  joined <- merge(tmp, pairt, by = 'cells', all = T)
  joined$l <- gtools::na.replace(joined$expr, 0)
  joined$l[joined$cellpop == pair$recpair] = 0
  joined$r <- gtools::na.replace(joined$expr, 0)
  joined$r[joined$cellpop == pair$ligpair] = 0
  sel <- unique(joined$cellpop)
  p2 <- ggplot2::ggplot(joined, aes(x = c1, y = c2, l)) +
    ggplot2::geom_point(data = subset(joined, cellpop == pair$ligpair), aes(color = l, alpha = l), size = 2) +
    ggplot2::scale_color_gradient(low = 'gray', high = 'blue') +
    ggplot2::labs(color = 'ligand') +
    ggnewscale::new_scale_color() +
    ggplot2::geom_point(data = subset(joined, cellpop == pair$recpair), aes(color = r, alpha = r), size = 2) +
    ggplot2::scale_color_gradient(low = 'gray', high = 'red') +
    ggplot2::labs(color = 'receptor') +
    ggnewscale::new_scale_color() +
    ggplot2::geom_point(data = subset(joined, !(cellpop %in% c(pair$ligpair, pair$recpair))), aes(color = expr), size = 2) +
    ggplot2::scale_color_gradient(low = 'black', high = 'black') +
    ggplot2::labs(color = 'unselected') +
    ggplot2::theme_minimal() +
    ggplot2::xlab('UMAP_1') +
    ggplot2::ylab('UMAP_2') +
    ggplot2::labs(title = 'LR pair expression')
  impulr <- Rmagic::magic(data, genes = c(pair$Ligand, pair$Receptor))
  Seurat::DefaultAssay(impulr) <- 'MAGIC_RNA'
  scells <- unique(na.omit(joined$cells[joined$expr != 0]))
  p3 <- Seurat::FeatureScatter(impulr, feature1 = pair$Ligand, feature2 = pair$Receptor, cells = scells, cols = ggsci::pal_d3()(2), pt.size = 2)
  p3 <- p3 + ggplot2::labs(title = 'Imputed Matrix Feature Plot')
  return(p1 +
           p2 +
           p3 * theme(plot.title = element_text(face = "plain")))
}
