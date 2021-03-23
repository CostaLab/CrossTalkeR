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
#'@importFrom tidyr %>%
#'@import colorBlindness
#'@return R default plot
#'@export
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
                    pg = NULL) {

  # Check Maximal Weight
  if (is.null(emax)) {
    emax <- max(abs(igraph::E(graph)$mean))
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
  }else{
    coords_scale <- coords
  }
  # It will make the loops in a correct angle
  loop_angle <- ifelse(coords_scale[igraph::V(graph)$name, 1] > 0,
                       -atan(coords_scale[igraph::V(graph)$name, 2] / coords_scale[igraph::V(graph)$name, 1]),
                       pi-atan(coords_scale[igraph::V(graph)$name, 2] / coords_scale[igraph::V(graph)$name, 1])
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
  if(is.null(pg)){
    pg <- igraph::page.rank(graph, weights = abs(igraph::E(graph)$weight),directed=T)$vector
    igraph::V(graph)$size <-  oce::rescale(pg,rlow = 0, rhigh = 60)
  }
  else{
    igraph::V(graph)$size <-  oce::rescale(pg,rlow = 0, rhigh = 75)
  }

  if (log) {
        igraph::E(graph)$width <- ifelse(igraph::E(graph)$inter != 0,
                                         log2(1 + igraph::E(graph)$inter),
                                         0) * efactor
  }else{
        igraph::E(graph)$width <- ifelse(igraph::E(graph)$inter != 0,
                                         igraph::E(graph)$inter,
                                         0) * efactor
  }
  igraph::E(graph)$arrow.size <- 0.6
  igraph::E(graph)$arrow.width <- igraph::E(graph)$width + 0.8
  if (sum(edge_start[, 2] == edge_start[, 1]) != 0) {
    igraph::E(graph)$loop.angle[which(edge_start[,2]==edge_start[,1])]<-loop_angle[edge_start[which(edge_start[,2]==edge_start[,1]),1]]
    igraph::E(graph)$loop.angle[which(edge_start[,2]!=edge_start[,1])]<-0
  }
  coords_scale[,1] <- scales::rescale(coords_scale[,1],from=c(-1,1),to=c(-2,2))
  coords_scale[,2] <- scales::rescale(coords_scale[,2],from=c(-1,1),to=c(-2,2))
  plot(graph,
       layout = coords_scale,
       xlim = c(-4, 4),
       ylim = c(-4, 4),
       rescale=F,
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
  a <- graphics::legend('bottomleft',
                        title="Node Pagerank",
                        legend=c("","",""),
                        pt.cex=c(min(v)+1,mean(v),max(v))/12,col='black',
               pch=21, pt.bg='black',box.lwd = 0,y.intersp=2)
  graphics::text(a$rect$left + a$rect$w, a$text$y,
                  c(round(min(pg),2),round(mean(pg),2),round(max(pg),2)), pos = 2)
  #x <- (a$text$x + a$rect$left) / 2
  #y <- a$text$y
  #symbols(x,y,circles=c(min(v),mean(v),max(v))/500,inches=FALSE,add=TRUE,bg='black')
  # Node names
  x <- coords_scale[, 1] * 1.2
  y <- coords_scale[, 2] * 1.2
  coord_ratio <- coords_scale[, 1] / coords_scale[, 2]
  angle <- ifelse(
              atan(-coord_ratio) * (180 / pi) < 0,
              90 + atan(-coord_ratio) * (180 / pi),
              270 + atan(-coord_ratio) * (180 / pi))
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
  if (leg) {
      # Edge Colormap
      if(min(igraph::E(graph)$mean) < 0){
        netdiffuseR::drawColorKey(seq(1, 200),
                                  tick.marks = c(1,101,200),
                                  color.palette = col_pallet,
                                  labels = c(-round(emax, 3),0,round(emax, 3)),
                                  nlevels = 200,
                                  main = "Weights",
                                  pos = 2,
                                  key.pos = c(0.98, 1.0, 0.0, 0.2),
                                  border = "transparent")
      }
      else{
        netdiffuseR::drawColorKey(seq(100, 200),
                                  tick.marks = c(100, 200),
                                  color.palette = col_pallet[100:201],
                                  labels = c(0, round(emax, 3)),
                                  nlevels =100,
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
#'@import ggplot2
#'@import ggraph
#'@import igraph
#'@import colorBlindness
#'@import graphlayouts
#'@return R default plot
#'@importFrom tidyr %>%
#'@export
plot_ggi <- function(graph,color,name) {
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
  ecolors <- rev(grDevices::colorRampPalette(colorBlindness::Blue2DarkOrange18Steps)(25))[lr2col]
  igraph::E(graph)$colors <- ecolors
  cls <- names(color)
  names(color) <- NULL
  ewidth <- abs(igraph::E(graph)$MeanLR) - mean(abs(igraph::E(graph)$MeanLR))
  ewidth <- ewidth / stats::sd(abs(igraph::E(graph)$MeanLR))
  ewidth <-  (ewidth - min(ewidth)) / (max(ewidth) - min(ewidth))
  print(ggraph::ggraph(graph, layout = "stress") +
        ggraph::geom_edge_link0(ggplot2::aes(edge_width = ewidth,
                                             color = igraph::E(graph)$MeanLR),
                                             alpha = ewidth) +
        ggraph::scale_edge_color_gradient2(low = grDevices::colorRampPalette(colorBlindness::Blue2DarkOrange18Steps)(25)[1],
                                           high = grDevices::colorRampPalette(colorBlindness::Blue2DarkOrange18Steps)(25)[25]) +
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
        ggplot2::ggtitle(name)+
        #scale_size(range = c(1,6))+
        ggraph::theme_graph() +
        ggplot2::theme(legend.position = "left"))

}


#'This function do a ggi plot and articulation
#'
#'@param graph graph
#'@param color cluster color
#'@import ggplot2
#'@import ggraph
#'@import graphlayouts
#'@import colorBlindness
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
  ecolors <- rev(grDevices::colorRampPalette(colorBlindness::Blue2DarkOrange18Steps)(25))[round(oce::rescale(igraph::E(graph)$MeanLR,
                                           xlow = min(igraph::E(graph)$MeanLR),
                                           rlow = 1,
                                           rhigh = 25))]
  igraph::E(graph)$colors <- ecolors
  cls <- names(color)
  names(color) <- NULL
  ewidth <- (abs(igraph::E(graph)$MeanLR) - mean(abs(igraph::E(graph)$MeanLR)))
  ewidth <- ewidth / stats::sd(abs(igraph::E(graph)$MeanLR))
  ewidth <-  (ewidth - min(ewidth)) / (max(ewidth) - min(ewidth))
  ap <- igraph::V(graph) %in% igraph::articulation.points(graph)
  print(ggraph::ggraph(graph, layout = "stress") +
        ggraph::geom_edge_link0(ggplot2::aes(edge_width = ewidth,
                                             color = igraph::E(graph)$MeanLR),
                                alpha = ewidth) +
        ggraph::scale_edge_color_gradient2(low = grDevices::colorRampPalette(colorBlindness::Blue2DarkOrange18Steps)(25)[1],
                                           high = grDevices::colorRampPalette(colorBlindness::Blue2DarkOrange18Steps)(25)[25]) +
        ggraph::geom_node_point(size = (deg / max(deg) * 10),
                                alpha = 1,
                                ggplot2::aes(color = igraph::V(graph)$cluster)
                                ) +
        ggplot2::scale_colour_manual(values = color,
                                     labels = cls,
                                     name = "Clusters") +
        ggraph::geom_node_label(ggplot2::aes(filter = ap,
                                            label = igraph::V(graph)$genes,
                                            color = igraph::V(graph)$cluster),
                                 show.legend = FALSE) +
        ggraph::scale_edge_width_continuous(range = c(0, 1)) +
        ggraph::theme_graph() +
        ggplot2::theme(legend.position = "left"))

}



#'This function selected genes sankey plot
#'
#'@param lrobj_tbl LRobject table with all data
#'@param target gene
#'@param ligand_cluster Ligand Clusters
#'@param receptor_cluster Receptor Clusters
#'@param plt_name plot title
#'@import ggplot2
#'@import dplyr
#'@import colorBlindness
#'@import ggalluvial
#'@importFrom tidyr %>%
#'@importFrom stats reorder
#'@return R default plot
#'@export
plot_sankey <- function(lrobj_tbl,
                        target = NULL,
                        ligand_cluster = NULL,
                        receptor_cluster = NULL,
                        plt_name = NULL,
                        threshold=50) {
  data <- lrobj_tbl
  if (!is.null(target)) {
    data <- lrobj_tbl[grepl(target, lrobj_tbl$allpair), ]
  }
  if (!is.null(ligand_cluster)) {
    tmp_sel <- grepl(ligand_cluster, data$Ligand.Cluster)
    data <- data[tmp_sel, ]
  }
  if (!is.null(receptor_cluster)) {
    tmp_sel <- grepl(receptor_cluster, data$Receptor.Cluster)
    data <- data[tmp_sel, ]
  }
  data$freq <- 1
  colp <-c(Blue2DarkOrange18Steps[4],Blue2DarkOrange18Steps[14])
  tmp_cols <- c("Ligand Cluster", "Ligand", "Receptor", "Receptor Cluster")
  names(colp) <- c("FALSE", "TRUE")
  if (dim(data)[1] >= 1) {
    tmp <- dplyr::top_n(data, ifelse(dim(data)[1] > threshold, threshold,
                        dim(data)[1]), abs(.data$MeanLR))
    print(ggplot2::ggplot(tmp, aes(y = .data$freq, axis1 = .data$Ligand.Cluster,
                                   axis2 = stats::reorder(.data$Ligand, .data$MeanLR),
                                   axis3 = stats::reorder(.data$Receptor, .data$MeanLR),
                                   axis4 = .data$Receptor.Cluster)) +
          ggalluvial::geom_alluvium(aes(fill = .data$MeanLR > 0),
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
  else{
      print(paste0("Gene->", target, "Not Found"))

  }
}




#'This function signed sending and receiving barplot (CCI)
#'
#'@param all_data LRobject
#'@param curr condition
#'@import ggplot2
#'@import dplyr
#'@import patchwork
#'@import colorBlindness
#'@importFrom tidyr %>%
#'@importFrom stats reorder
#'@return R default plot
#'@export
plot_signedbar <- function(all_data,curr) {
   curr_net <- all_data@graphs[[curr]]
   in_deg_up <- table(all_data@tables[[curr]]$Ligand.Cluster[all_data@tables[[curr]]$MeanLR > 0])
   in_up <- tibble::tibble(as.data.frame(in_deg_up))
   in_deg_down <- table(all_data@tables[[curr]]$Ligand.Cluster[all_data@tables[[curr]]$MeanLR < 0])
   in_down <- tibble::tibble(as.data.frame(in_deg_down))
   in_down$Freq <- 0-in_down$Freq
   in_all <- dplyr::bind_rows(in_up,in_down)
   in_all$Expression <- ifelse(in_all$Freq<0,'Downregulated','Upregulated')
   in_all$rank <- ifelse(in_all$Freq<0,0,1)
   out_deg_up <- table(all_data@tables[[curr]]$Receptor.Cluster[all_data@tables[[curr]]$MeanLR > 0])
   out_up <- tibble::tibble(as.data.frame(out_deg_up))
   out_deg_down <- table(all_data@tables[[curr]]$Receptor.Cluster[all_data@tables[[curr]]$MeanLR < 0])
   out_down <- tibble::tibble(as.data.frame(out_deg_down))
   out_down$Freq <- 0-out_down$Freq
   out_all <- dplyr::bind_rows(out_up,out_down)
   out_all$Expression <- ifelse(out_all$Freq<0,'Downregulated','Upregulated')
   out_all$rank <- ifelse(out_all$Freq<0,0,1)
   p1 <- ggplot2::ggplot(in_all,ggplot2::aes(x=.data$Freq,y=reorder(.data$Var1,.data$Freq*.data$rank),fill=.data$Expression))+
     ggplot2::geom_bar(stat = 'identity',position = "identity")+
     ggplot2::geom_text(ggplot2::aes(label=.data$Freq),size=3.5)+
     ggplot2::scale_fill_manual(values=grDevices::colorRampPalette(colorBlindness::Green2Magenta16Steps)(2))+
     ggplot2::ggtitle('Ligands')+
     ggplot2::ylab('Cell')+
     ggplot2::xlab('Number of interactions')+
     ggplot2::theme_minimal()
   p2 <- ggplot2::ggplot(out_all,ggplot2::aes(x=.data$Freq,y=reorder(.data$Var1,.data$Freq*.data$rank),fill=.data$Expression))+
     ggplot2::geom_bar(stat = 'identity',position = "identity")+
     ggplot2::geom_text(ggplot2::aes(label=.data$Freq),size=3.5)+
     ggplot2::scale_fill_manual(values=grDevices::colorRampPalette(colorBlindness::Blue2DarkOrange18Steps)(2))+
     ggplot2::ggtitle('Receptors')+
     ggplot2::ylab('Cell')+
     ggplot2::xlab('Number of interactions')+
     ggplot2::theme_minimal()
   print((p1+p2)+patchwork::plot_annotation(title = curr,tag_levels = 'A'))
}



#'This function signed sending and receiving barplot (GGI)
#'
#'@param all_data LRobject
#'@param curr condition
#'@import ggplot2
#'@import dplyr
#'@import patchwork
#'@import colorBlindness
#'@importFrom tidyr %>%
#'@importFrom stats reorder
#'@return R default plot
#'@export
plot_signedbar_ggi <- function(all_data,curr) {
  curr_net <- all_data@graphs_ggi[[curr]]
  up_graph <- igraph::subgraph.edges(curr_net, E(curr_net)[E(curr_net)$MeanLR > 0])
  down_graph <- igraph::subgraph.edges(curr_net, E(curr_net)[E(curr_net)$MeanLR < 0])
  in_deg_up <- igraph::degree(up_graph, mode = 'in')
  in_deg_down <- igraph::degree(down_graph, mode = 'in')
  in_up <- tibble::tibble(genes = paste0(names(in_deg_up),'_up'), values=as.array(in_deg_up))
  in_down <- tibble::tibble(genes = paste0(names(in_deg_down),'_down'), values=as.array(in_deg_down))
  in_deg_data_up <-dplyr::top_n(in_up, 10, .data$values)
  in_deg_data_down <- dplyr::top_n(in_down, 10, .data$values)
  in_deg_data_down$values <- 0 -in_deg_data_down$values
  in_deg_data <- dplyr::bind_rows(in_deg_data_up,in_deg_data_down )
  in_deg_data$Expression <- ifelse(in_deg_data$values <0,'Downregulated','Upregulated')
  p1 <- ggplot2::ggplot(in_deg_data,ggplot2::aes(x=.data$values,y=reorder(.data$genes,.data$values),fill=.data$Expression))+
    ggplot2::geom_bar(stat = 'identity',position = "identity")+
    ggplot2::scale_fill_manual(values=c(Blue2DarkOrange18Steps[4],Blue2DarkOrange18Steps[14]))+
    ggplot2::ggtitle('Receiving')+
    ggplot2::ylab('Cell')+
    ggplot2::xlab('Number of interactions')+
    ggplot2::theme_minimal()
  out_deg_up <- igraph::degree(up_graph, mode = 'out')
  out_deg_down <- igraph::degree(down_graph, mode = 'out')
  out_up <- tibble::tibble(genes = paste0(names(out_deg_up),'_up'), values=as.array(out_deg_up))
  out_down <- tibble::tibble(genes = paste0(names(out_deg_down),'_down'), values=as.array(out_deg_down))
  out_deg_data_up <-dplyr::top_n(out_up, 10, .data$values)
  out_deg_data_down <- dplyr::top_n(out_down, 10, .data$values)
  out_deg_data_down$values <- 0-out_deg_data_down$values
  out_deg_data <- dplyr::bind_rows(out_deg_data_up,out_deg_data_down )
  out_deg_data$Expression <- ifelse(out_deg_data$values <0,'Downregulated','Upregulated')
  p2 <- ggplot2::ggplot(out_deg_data,ggplot2::aes(x=.data$values,y=reorder(.data$genes,.data$values),fill=.data$Expression))+
    ggplot2::geom_bar(stat = 'identity',position = "identity")+
    ggplot2::scale_fill_manual(values=c(Blue2DarkOrange18Steps[4],Blue2DarkOrange18Steps[14]))+
    ggplot2::ggtitle('Sending')+
    ggplot2::ylab('Cell')+
    ggplot2::xlab('Number of interactions')+
    ggplot2::theme_minimal()
  print((p2+p1)+plot_annotation(title = curr,tag_levels = 'A'))
}
