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
    igraph::V(graph)$size <-  60
  }
  else{
    igraph::V(graph)$size <-  scales::rescale(pg,c(1,60))
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
  igraph::E(graph)$arrow.size <- 0.4
  igraph::E(graph)$arrow.width <- igraph::E(graph)$width+0.8
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
  if(!is.null(pg)){
    a <- graphics::legend('bottomleft',
                          title="Node Pagerank",
                          legend=c("","",""),
                          pt.cex=c(min(v)+1,mean(v),max(v))/12,col='black',
                 pch=21, pt.bg='black',box.lwd = 0,y.intersp=2)
    graphics::text(a$rect$left + a$rect$w, a$text$y,
                    c(round(min(pg),2),round(mean(pg),2),round(max(pg),2)), pos = 2)
  }
  x <- coords_scale[, 1] * 1.2
  y <- coords_scale[, 2] * 1.2
  coord_ratio <- coords_scale[, 1] / coords_scale[, 2]
  angle <- ifelse(
              atan(-coord_ratio) * (180 / pi) < 0,
              90 + atan(-coord_ratio) * (180 / pi),
              270 + atan(-coord_ratio) * (180 / pi))
  if(vnames){
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
      if(min(igraph::E(graph)$weight) < 0){
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
plot_ggi <- function(graph,color,name) {
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
  ewidth <-  (ewidth - min(ewidth)) / (max(ewidth) - min(ewidth))
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
        ggraph::scale_edge_width_continuous(range = c(0, 1))  +
        ggplot2::ggtitle(name)+
        #scale_size(range = c(1,6))+
        ggraph::theme_graph(base_family="sans") +
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
                        threshold=50) {

  if (!is.null(target)) {
      data <- lrobj_tbl[grepl(target, lrobj_tbl$allpair), ]
  }
  else{
      data <- lrobj_tbl
  }
  if (!is.null(ligand_cluster)) {
    tmp_sel <- grepl(ligand_cluster, data$ligpair)
    data <- data[tmp_sel, ]
  }
  if (!is.null(receptor_cluster)) {
    tmp_sel <- grepl(receptor_cluster, data$recpair)
    data <- data[tmp_sel, ]
  }
  colp <-c(Blue2DarkOrange18Steps[4],Blue2DarkOrange18Steps[14])
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
          ggalluvial::geom_alluvium(aes(fill = .data$LRScore > 0,color='b'),
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
  else{
      print(paste0("Gene->", target, "Not Found"))
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
plotInterInfo<-function(data,lrobject,pair,lrslot,assay='RNA'){
      Seurat::DefaultAssay(data) <- assay
      pair <- lrobject@tables[[lrslot]][grepl(pair,datalr@tables$EXP_x_CTR$allpair),]
      sender <- pair$Ligand.Cluster
      receiver <- pair$Receptor.Cluster
      exp <- Seurat::GetAssayData(data)
      lpair <- tibble::tibble(cells=names(exp[pair$Ligand,grepl(sender,Idents(data))]),
                      expr=exp[pair$Ligand,grepl(sender,Idents(data))],
                      cellpop=pair$ligpair,
                      stage=data@meta.data[names(exp[pair$Ligand,grepl(sender,Idents(data))]),]$protocol)
      rpair <- tibble::tibble(cells=names(exp[pair$Receptor,grepl(receiver,Idents(data))]),
                      expr=exp[pair$Receptor,grepl(receiver,Idents(data))],
                      cellpop=pair$recpair,
                      stage=data@meta.data[names(exp[pair$Receptor,grepl(receiver,Idents(data))]),]$protocol)
      pairt <- dplyr::bind_rows(lpair,rpair)
      p1<-ggplot2::ggplot(pairt,aes(x=expr,y=cellpop,fill=stage))+
          ggridges::geom_density_ridges()+
          ggplot2::facet_grid(.~stage)+
          ggplot2::theme_minimal()+
          ggplot2::labs(title='Expression density plot')
      tmp <- Seurat::Embeddings(data,'umap')
      tmp <- tibble::tibble(c1 = tmp[,"UMAP_1"],c2 = tmp[,"UMAP_2"],cells=rownames(tmp))
      joined <- merge(tmp,pairt,by='cells',all=T)
      joined$l <- gtools::na.replace(joined$expr,0)
      joined$l[joined$cellpop==pair$recpair] = 0
      joined$r <- gtools::na.replace(joined$expr,0)
      joined$r[joined$cellpop==pair$ligpair] = 0
      sel <- unique(joined$cellpop)
      p2<-ggplot2::ggplot(joined,aes(x=c1,y=c2,l))+
        ggplot2::geom_point(data=subset(joined,cellpop == pair$ligpair),aes(color=l,alpha=l),size=2)+
        ggplot2::scale_color_gradient(low='gray',high='blue')+
        ggplot2::labs(color='ligand')+
        ggnewscale::new_scale_color()+
        ggplot2::geom_point(data=subset(joined,cellpop == pair$recpair),aes(color=r,alpha=r),size=2)+
        ggplot2::scale_color_gradient(low='gray',high='red')+
        ggplot2::labs(color='receptor')+
        ggnewscale::new_scale_color()+
        ggplot2::geom_point(data=subset(joined,!(cellpop %in% c(pair$ligpair,pair$recpair))),aes(color=expr),size=2)+
        ggplot2::scale_color_gradient(low='black',high='black')+
        ggplot2::labs(color='unselected')+
        ggplot2::theme_minimal()+ggplot2::xlab('UMAP_1')+ggplot2::ylab('UMAP_2')+ggplot2::labs(title='LR pair expression')
      impulr<- Rmagic::magic(data,genes=c(pair$Ligand,pair$Receptor))
      Seurat::DefaultAssay(impulr) <- 'MAGIC_RNA'
      scells <- unique(na.omit(joined$cells[joined$expr!=0]))
      p3<-Seurat::FeatureScatter(impulr,feature1 = pair$Ligand,feature2 = pair$Receptor,cells =scells,cols = ggsci::pal_d3()(2),pt.size = 2)
      p3<-p3+ggplot2::labs(title='Imputed Matrix Feature Plot')
      return(p1+p2+p3*theme(plot.title = element_text(face = "plain")))
  }
