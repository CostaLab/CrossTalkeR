#'Plot Cell Cell Interaction
#'
#'This function do a CCI plot
#'
#'@param graph Paths of single condition LR data
#'@param col Cell type (Cluster) Colors
#'@param plt_name Plot Name (Title)
#'@param emax Max MeanLR across the all inputs, if its not defined, the method going to consider the max find within a sample
#'@param col Cluster coords according do desired layout default is circle layout
#'@param leg Set color legend
#'@param low Lower threshold: This parameter low and high defines the edges which will be filtered. Edges within the interval [low,high] are filtered.
#'@param high Higher threshould
#'@param efactor edge scale factor
#'@param vfactor edge scale factor
#'@return R default plot
#'@importFrom tidyr %>%
#'@export
plot_cci <- function(graph,
                    col,
                    plt_name,
                    coords,
                    emax = NULL,
                    leg = FALSE,
                    low = 25,
                    high = 75,
                    ignore_alpha = F,
                    log = F,
                    efactor = 8,
                    vfactor = 12) {
  # Check Maximal Weight
  if (is.null(emax)) {
    emax <- max(abs(igraph::E(graph)$mean))
  }
  # Using color pallet to up and down regulation
  col_pallet <- pals::coolwarm(9)
  # Expanding the pallet range
  col_pallet <- colorRampPalette(col_pallet)(201)
  # Checking looops
  edge_start <- igraph::ends(graph,
                             es = igraph::E(graph),
                             names = FALSE)
  # Scale nodes coordinates
  if (nrow(coords) != 1) {
    coords_scale <- scale(coords)
  }else{
    coords_scale <- coords
  }
  # It will make the loops in a correct angle
  loop_angle <- ifelse(coords_scale[igraph::V(graph)$name, 1] > 0,
                       -atan(coords_scale[igraph::V(graph)$name, 2] / coords_scale[igraph::V(graph)$name, 1]),
                       pi-atan(coords_scale[igraph::V(graph)$name, 2] / coords_scale[igraph::V(graph)$name, 1])
                       )
  # Setting node colors
  igraph::V(graph)$color <- col[igraph::V(graph)$name]
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
  max_deg <- max(igraph::degree(subgraph, normalized = T))
  igraph::V(graph)$size <- (igraph::degree(subgraph, normalized = T) / max_deg)
  igraph::V(graph)$size <- igraph::V(graph)$size * vfactor
  if (log) {
        igraph::E(graph)$width <- ifelse(igraph::E(graph)$inter != 0,
                                         log2(1 + igraph::E(graph)$inter),
                                         0) * efactor
  }else{
        igraph::E(graph)$width <- ifelse(igraph::E(graph)$inter != 0,
                                         igraph::E(graph)$inter,
                                         0) * efactor
  }
  igraph::E(graph)$arrow.size <- 0.25
  igraph::E(graph)$arrow.width <- igraph::E(graph)$width + 0.6
  if (sum(edge_start[, 2] == edge_start[, 1]) != 0) {
    igraph::E(graph)$loop.angle[which(edge_start[,2]==edge_start[,1])]<-loop_angle[edge_start[which(edge_start[,2]==edge_start[,1]),1]]
    igraph::E(graph)$loop.angle[which(edge_start[,2]!=edge_start[,1])]<-0
  }
  plot(graph,
       layout = coords_scale,
       xlim = c(-1.5, 1.5),
       ylim = c(-1.5, 1.5),
       edge.curved = 0.5,
       vertex.label = NA,
       vertex.shape = "circle",
       margin = 0.0,
       loop.angle = igraph::E(graph)$loop.angle,
       edge.label = NA,
       main = plt_name
       )
  # Thicknesse legenf
  amin <- min(igraph::E(graph)$inter[igraph::E(graph)$inter != 0])
  amax <- max(igraph::E(graph)$inter)
  e_wid_sp <- c(amin,
                amin + amax / 2,
                amax)
  legend("topleft",
         legend = round(e_wid_sp, 1),
         col = "black",
         title = "Percentage of the interactions",
         pch = NA,
         bty = "n",
         cex = 1,
         lwd = e_wid_sp,
         lty = c(1, 1, 1),
         horiz = FALSE)
  # Node names
  x <- coords_scale[, 1] * 1.2
  y <- coords_scale[, 2] * 1.2
  coord_ratio <- coords_scale[, 1] / coords_scale[, 2]
  angle <- ifelse(
              atan(-coord_ratio) * (180 / pi) < 0,
              90 + atan(-coord_ratio) * (180 / pi),
              270 + atan(-coord_ratio) * (180 / pi))
  for (i in seq_len(length(x))) {
    text(x = x[i],
         y = y[i],
         labels = igraph::V(graph)$name[i],
         adj = NULL,
         pos = NULL,
         cex = 0.8,
         col = "black",
         xpd = T)
  }
  if (leg) {
      # Edge Colormap
      netdiffuseR::drawColorKey(seq(1, 200),
                                tick.marks = c(1, 200),
                                color.palette = col_pallet,
                                labels = c(-round(emax, 3), round(emax, 3)),
                                nlevels = 200,
                                main = "Weights",
                                pos = 2,
                                key.pos = c(0.98, 1.0, 0.0, 0.2),
                                border = "transparent")
  }
}



#'This function do a ggi plot and higest degree nodes
#'
#'@param graph graph
#'@param color cluster color
#'@import ggplot2
#'@import ggraph
#'@import igraph
#'@return R default plot
#'@importFrom tidyr %>%
#'@export
plot_ggi <- function(graph, color) {
  deg <- igraph::degree(graph)
  v_names <- igraph::V(graph)$name
  v_names <- tibble::as_tibble(v_names)
  v_names <- v_names %>%
    tidyr::separate(.data$value, c("genes", "cluster"), "/")
  igraph::V(graph)$genes <- v_names$genes
  igraph::V(graph)$cluster <- v_names$cluster
  lr2col <- round(oce::rescale(igraph::E(graph)$MeanLR,
                               xlow = min(igraph::E(graph)$MeanLR),
                               rlow = 1,
                               rhigh = 25
                               )
                  )
  ecolors <- rev(pals::coolwarm(25))[lr2col]
  igraph::E(graph)$colors <- ecolors
  cls <- names(color)
  names(color) <- NULL
  ewidth <- abs(igraph::E(graph)$MeanLR) - mean(abs(igraph::E(graph)$MeanLR))
  ewidth <- ewidth / sd(abs(igraph::E(graph)$MeanLR))
  ewidth <-  (ewidth - min(ewidth)) / (max(ewidth) - min(ewidth))
  print(ggraph::ggraph(graph, layout = "stress") +
        ggraph::geom_edge_link0(ggplot2::aes(edge_width = ewidth,
                                             color = igraph::E(graph)$MeanLR),
                                             alpha = ewidth) +
        ggraph::scale_edge_color_gradient2(low = pals::coolwarm(25)[1],
                                           high = pals::coolwarm(25)[25]) +
        ggraph::geom_node_point(size = ((deg / max(deg)) * 10),
                                alpha = 1,
                                ggplot2::aes(color = igraph::V(graph)$cluster)
                                ) +
        ggplot2::scale_colour_manual(values = color,
                                     labels = cls,
                                     name = "Clusters") +
        ggraph::geom_node_label(ggplot2::aes(filter = deg > deg[order(deg, decreasing = T)][ifelse(length(deg) > 100, 100, length(deg))] & igraph::V(graph) %in% igraph::articulation.points(graph),
                                             label = igraph::V(graph)$genes,
                                             color = igraph::V(graph)$cluster
                                            ),
                               repel = TRUE,
                               fontface = "italic",
                               hjust = "inward",
                               size = 7,
                               show_guide = F) +
        ggraph::scale_edge_width_continuous(range = c(0, 1)) +
        #scale_size(range = c(1,6))+
        ggraph::theme_graph() +
        ggplot2::theme(legend.position = "left"))

}

### Parei aqui

#'This function do a ggi plot and articulation
#'
#'@param graph graph
#'@param color cluster color
#'@import ggplot2
#'@import ggraph
#'@return R default plot
#'@importFrom tidyr %>%
#'@export
plot_articulation <- function(graph, color) {
  deg <- igraph::degree(graph)
  v_names <- igraph::V(graph)$name
  v_names <- tibble::as_tibble(v_names)
  v_names <- v_names %>%
    tidyr::separate(.data$value, c("genes", "cluster"), "/")
  igraph::V(graph)$genes <- v_names$genes
  igraph::V(graph)$cluster <- v_names$cluster
  ecolors <- rev(pals::coolwarm(25))[round(oce::rescale(igraph::E(graph)$MeanLR,
                                           xlow = min(igraph::E(graph)$MeanLR),
                                           rlow = 1,
                                           rhigh = 25))]
  igraph::E(graph)$colors <- ecolors
  cls <- names(color)
  names(color) <- NULL
  ewidth <- (abs(igraph::E(graph)$MeanLR) - mean(abs(igraph::E(graph)$MeanLR)))
  ewidth <- ewidth / sd(abs(igraph::E(graph)$MeanLR))
  ewidth <-  (ewidth - min(ewidth)) / (max(ewidth) - min(ewidth))
  ap <- igraph::V(graph) %in% igraph::articulation.points(graph)
  print(ggraph::ggraph(graph, layout = "stress") +
        ggraph::geom_edge_link0(ggplot2::aes(edge_width = ewidth,
                                             color = igraph::E(graph)$MeanLR),
                                alpha = ewidth) +
        ggraph::scale_edge_color_gradient2(low = pals::coolwarm(25)[1],
                                           high = pals::coolwarm(25)[25]) +
        ggraph::geom_node_point(size = (deg / max(deg) * 10),
                                alpha = 1,
                                ggplot2::aes(color = igraph::V(graph)$cluster)
                                ) +
        ggplot2::scale_colour_manual(values = color,
                                     labels = cls,
                                     name = "Clusters") +
        ggraph::geom_node_text(ggplot2::aes(filter = ap,
                                            label = igraph::V(graph)$genes)) +
        ggraph::scale_edge_width_continuous(range = c(0, 1)) +
        ggraph::theme_graph() +
        ggplot2::theme(legend.position = "left"))

}



#'This function selected genes sankey plot
#'
#'@param lrobj_tbl LRobject with all data
#'@param target gene
#'@param lcls Ligand Clusters
#'@param rcls Receptor Clusters
#'@import ggplot2
#'@import dplyr
#'@import ggalluvial
#'@importFrom ggalluvial StatStratum
#'@return R default plot
#'@importFrom tidyr %>%
#'@export
plot_sankey <- function(lrobj_tbl,
                        target = NULL,
                        ligand_cluster = NULL,
                        receptor_cluster = NULL,
                        pltname = NULL) {
  data <- lrobj_tbl
  StatStratum <- ggalluvial::StatStratum
  if (!is.null(target)) {
    data <- lrobj_tbl[grepl(target, lrobj_tbl$allpair), ]
  }
  if (!is.null(ligand_cluster)) {
    tmp_sel <- match(data$Ligand.Cluster, ligand_cluster, nomatch = F) != F
    data <- data[tmp_sel, ]
  }
  if (!is.null(receptor_cluster)) {
    tmp_sel <- match(data$Receptor.Cluster, receptor_cluster, nomatch = F) != F
    data <- data[tmp_sel, ]
  }
  data$freq <- 1
  colp <- pals::coolwarm(2)
  tmp_cols <- c("Ligand Cluster", "Ligand", "Receptor", "Receptor Cluster")
  names(colp) <- c("FALSE", "TRUE")
  if (dim(data)[1] >= 1) {
    tmp <- dplyr::top_n(data, ifelse(dim(data)[1] > 50, 50,
                        dim(data)[1]), abs(.data$MeanLR))
    print(ggplot2::ggplot(tmp, aes(y = .data$freq, axis1 = .data$Ligand.Cluster,
                                   axis2 = reorder(.data$Ligand, .data$MeanLR),
                                   axis3 = reorder(.data$Receptor, .data$MeanLR),
                                   axis4 = .data$Receptor.Cluster)) +
          ggalluvial::geom_alluvium(aes(fill = .data$MeanLR > 0),
                                    width = 1 / 12,
                                    discern = F) +
          ggalluvial::geom_stratum(width = 1 / 12) +
          ggplot2::geom_label(stat = "stratum",
                              ggplot2::aes(label = ggplot2::after_stat(.data$stratum)),
                              size = 2) +
          ggplot2::scale_x_discrete(limits = tmp_cols, expand = c(.05, .05)) +
          ggplot2::scale_fill_manual(values = colp,
                                     limits = names(colp),
                                     name = "Upregulated") +
          ggplot2::ggtitle(pltname) +
          ggplot2::theme(text = element_text(size = 8)) +
          ggplot2::theme_minimal()
          )
  }
  else{
      print(paste0("Gene->", target, "Not Found"))

  }
}
