#require(igraph)
#require(ggraph)
#library(oce)
#library("RColorBrewer")
#require(pals)
#require(scales)
#require(patchwork)
#'Plot Cell Cell Interaction
#'
#'This function do a CCI plot
#'
#'@param ce Paths of single condition LR data
#'@param col Cell type (Cluster) Colors
#'@param plt_name Plot Name (Title)
#'@param emax Max MeanLR across the all inputs, if its not defined, the method going to consider the max find within a sample
#'@param col Cluster coords according do desired layout default is circle layout
#'@param leg Set color legend
#'@param low Lower threshold: This parameter low and high defines the edges which will be filtered. Edges within the interval [low,high] are filtered.
#'@param high Higher threshould
#'@return R default plot
#'@importFrom tidyr %>%
#'@export
plot_cci <- function(ce,col,plt_name,coords,emax=NULL,leg=FALSE,low=25,high=75,ignore_alpha=F,log=F){
  if(is.null(emax)){ # Check Maximal Weight
    emax <-max(abs(igraph::E(ce)$mean))
  }
  col_pallet <- pals::coolwarm(9) # Using color pallet to up and down regulation
  col_pallet <- colorRampPalette(col_pallet)(201) # Expanding the pallet range
  edge.start <- igraph::ends(ce, es=igraph::E(ce), names=FALSE) # Checking looops
  if(nrow(coords)!=1){ # Scale nodes coordinates
    coords_scale=scale(coords)
  }else{
    coords_scale<-coords
  }
  # It will make the loops in a correct angle
  loop.angle<-ifelse(coords_scale[igraph::V(ce)$name,1]>0,
                     -atan(coords_scale[igraph::V(ce)$name,2]/coords_scale[igraph::V(ce)$name,1]),
                     pi-atan(coords_scale[igraph::V(ce)$name,2]/coords_scale[igraph::V(ce)$name,1]))
  # Setting node colors
  igraph::V(ce)$color<-col[igraph::V(ce)$name]
  # Setting label settings
  #graph::V(ce)$label.color<-'black'
  #igraph::V(ce)$label.cex<-1.5
  # Setting edge settings
  ## Color scheme
  we <- round(oce::rescale(igraph::E(ce)$weight,xlow=(-emax),xhigh=emax,rlow=1,rhigh=200,clip=TRUE),0)
  igraph::E(ce)$color <- col_pallet[we]
  alpha <- ifelse((igraph::E(ce)$inter > low) & (igraph::E(ce)$inter < high), 0, igraph::E(ce)$inter)
  subgraph <- igraph::delete.edges(ce, igraph::E(ce)[alpha==0 | is.na(alpha)])
  if(!ignore_alpha){
    igraph::E(ce)$color <- scales::alpha(igraph::E(ce)$color, alpha)
  }
  ## Thickness and arrow size
  igraph::V(ce)$size<-igraph::degree(subgraph,normalized = T)/max(igraph::degree(subgraph,normalized = T))*10
  if(log){
        igraph::E(ce)$width <- ifelse(igraph::E(ce)$inter!=0,log2(1+igraph::E(ce)$inter),0)*5#abs(E(ce)$weight)
  }else{
        igraph::E(ce)$width <- ifelse(igraph::E(ce)$inter!=0,igraph::E(ce)$inter,0)*5#abs(E(ce)$weight)
  }
  igraph::E(ce)$arrow.size<-0.25
  igraph::E(ce)$arrow.width<-igraph::E(ce)$width+0.6
  if(sum(edge.start[,2]==edge.start[,1])!=0){
      igraph::E(ce)$loop.angle[which(edge.start[,2]==edge.start[,1])]<-loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
      igraph::E(ce)$loop.angle[which(edge.start[,2]!=edge.start[,1])]<-0
  }
  plot(ce,
       layout=coords_scale,
       xlim=c(-1.5,1.5),
       ylim=c(-1.5,1.5),
       edge.curved=0.5,
       vertex.label=NA,
       vertex.shape='circle',
       margin=0.0,
       loop.angle=igraph::E(ce)$loop.angle,
       edge.label=NA,
       main=plt_name
       )
  # Thicknesse legenf
  e_wid_sp <- c(min(igraph::E(ce)$inter[igraph::E(ce)$inter!=0]),(min(igraph::E(ce)$inter[igraph::E(ce)$inter!=0])+max(igraph::E(ce)$inter))/2,max(igraph::E(ce)$inter))
  legend('topleft', legend=round(e_wid_sp,1),col='black',title='Percentage of the interactions',pch = NA, bty = "n", cex =1, lwd=e_wid_sp,lty=c(1,1,1),horiz = FALSE)
  # Node names
  x = coords_scale[,1]*1.2
  y = coords_scale[,2]*1.2
  angle = ifelse(atan(-(coords_scale[,1]/coords_scale[,2]))*(180/pi) < 0,  90 + atan(-(coords_scale[,1]/coords_scale[,2]))*(180/pi), 270+atan(-coords_scale[,1]/coords_scale[,2])*(180/pi))
  for (i in 1:length(x)) {
    text(x=x[i], y=y[i], labels=igraph::V(ce)$name[i], adj=NULL, pos=NULL, cex=0.8, col="black", xpd=T)#,srt=angle[i], xpd=T)
  }
  if(leg){
      # Edge Colormap
      netdiffuseR::drawColorKey(seq(1,200),tick.marks = c(1,200), color.palette=col_pallet,labels=c(-round(emax,3),round(emax,3)),nlevels=200,main="LR Interactions",pos=2,key.pos = c(0.98,1.0, 0.0, 0.2),border='transparent')
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
plot_ggi<-function(graph,color){
  deg <- igraph::degree(graph)
  #root <-names(deg[order(deg,decreasing = T)])[1] # this will be used in the root based node
  v_names <- igraph::V(graph)$name
  v_names <- tibble::as_tibble(v_names)
  v_names <-v_names %>%
    tidyr::separate(.data$value, c("genes", "cluster"), "_")
  igraph::V(graph)$genes <- v_names$genes
  igraph::V(graph)$cluster <- v_names$cluster
  ecolors <- rev(pals::coolwarm(25))[round(oce::rescale(igraph::E(graph)$MeanLR,xlow = min(igraph::E(graph)$MeanLR),rlow = 1, rhigh = 25))]
  igraph::E(graph)$colors <- ecolors
  cls <- names(color)
  names(color) <- NULL
  ewidth <- (abs(igraph::E(graph)$MeanLR) - mean(abs(igraph::E(graph)$MeanLR)))/sd(abs(igraph::E(graph)$MeanLR))
  ewidth <-  (ewidth-min(ewidth))/(max(ewidth)-min(ewidth))
  print(ggraph::ggraph(graph,layout='stress')+
          ggraph::geom_edge_link0(ggplot2::aes(edge_width=ewidth,color=igraph::E(graph)$MeanLR),alpha=ewidth)+
          ggraph::scale_edge_color_gradient2(low=pals::coolwarm(25)[1],high=pals::coolwarm(25)[25])+
          ggraph::geom_node_point(size=(deg/max(deg)*10),alpha=1,ggplot2::aes(color=igraph::V(graph)$cluster))+
          ggplot2::scale_colour_manual(values=color,labels=cls,name='Clusters')+
          ggraph::geom_node_label(ggplot2::aes(filter=deg>deg[order(deg,decreasing = T)][ifelse(length(deg)>100, 100, length(deg))] & V(graph) %in% igraph::articulation.points(graph),
                                               label=igraph::V(graph)$genes,color=igraph::V(graph)$cluster),
                                 repel = TRUE,
                                 fontface = "italic",
                                 hjust = "inward",size = 5)+
          ggraph::scale_edge_width_continuous(range = c(0,1))+
          #scale_size(range = c(1,6))+
          ggraph::theme_graph()+
          ggplot2::theme(legend.position = "left"))

}

#'This function do a ggi plot and articulation
#'
#'@param graph graph
#'@param color cluster color
#'@import ggplot2
#'@import ggraph
#'@return R default plot
#'@importFrom tidyr %>%
#'@export
plot_articulation<-function(graph,color){
  deg <- igraph::degree(graph)
  v_names <- igraph::V(graph)$name
  v_names <- tibble::as_tibble(v_names)
  v_names <-v_names %>%
    tidyr::separate(.data$value, c("genes", "cluster"), "_")
  igraph::V(graph)$genes <- v_names$genes
  igraph::V(graph)$cluster <- v_names$cluster
  ecolors <- rev(pals::coolwarm(25))[round(oce::rescale(igraph::E(graph)$MeanLR,xlow = min(igraph::E(graph)$MeanLR),rlow = 1, rhigh = 25))]
  igraph::E(graph)$colors <- ecolors
  cls <- names(color)
  names(color) <- NULL
  ewidth <- (abs(igraph::E(graph)$MeanLR) - mean(abs(igraph::E(graph)$MeanLR)))/sd(abs(igraph::E(graph)$MeanLR))
  ewidth <-  (ewidth-min(ewidth))/(max(ewidth)-min(ewidth))
  print(ggraph::ggraph(graph,layout='stress')+
        ggraph::geom_edge_link0(ggplot2::aes(edge_width=ewidth,color=igraph::E(graph)$MeanLR),alpha=ewidth)+
        ggraph::scale_edge_color_gradient2(low=pals::coolwarm(25)[1],high=pals::coolwarm(25)[25])+
        ggraph::geom_node_point(size=(deg/max(deg)*10),alpha=1,ggplot2::aes(color=igraph::V(graph)$cluster))+
        ggplot2::scale_colour_manual(values=color,labels=cls,name='Clusters')+
        ggraph::geom_node_text(ggplot2::aes(filter=igraph::V(graph) %in% igraph::articulation.points(graph),label=igraph::V(graph)$genes))+
        ggraph::scale_edge_width_continuous(range = c(0,1))+
        #scale_size(range = c(1,6))+
        ggraph::theme_graph()+
        ggplot2::theme(legend.position = "left"))

}



#'This function selected genes sankey plot
#'
#'@param LRObj_tbl LRobject with all data
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
plot_sankey<-function(LRObj_tbl,target=NULL,Ligand.Cluster=NULL,Receptor.Cluster=NULL,pltname=NULL){
  data <- LRObj_tbl
  if(!is.null(target)){
    data <- LRObj_tbl[grepl(target, LRObj_tbl$allpair), ]
  }
  if(!is.null(Ligand.Cluster)){
    data <- data[match(data$Ligand.Cluster, Ligand.Cluster,nomatch = F)!=F,]
  }
  if(!is.null(Ligand.Cluster)){
    data <- data[match(data$Receptor.Cluster, Receptor.Cluster,nomatch = F)!=F, ]
  }
  data$Freq <- 1
  colp <- pals::coolwarm(2)
  names(colp) <- c('FALSE','TRUE')
  if(dim(data)[1] >= 1){
    tmp <- dplyr::top_n(data, ifelse(dim(data)[1] > 50, 50, dim(data)[1]), abs(.data$MeanLR))
    print(ggplot2::ggplot(tmp, aes(y = .data$Freq, axis1 = .data$Ligand.Cluster, axis2 = reorder(.data$Ligand,.data$MeanLR), axis3 = reorder(.data$Receptor,.data$MeanLR), axis4 = .data$Receptor.Cluster)) +
          ggalluvial::geom_alluvium(aes(fill = .data$MeanLR > 0), width = 1/12,discern=F) +
          ggalluvial::geom_stratum(width = 1/12) +
          # geom_text(stat = "stratum", aes(label = after_stat(.data$stratum)), size = 2) +
          ggplot2::geom_label(stat = "stratum", ggplot2::aes(label = ggplot2::after_stat(.data$stratum)), size = 2) +
          ggplot2::scale_x_discrete(limits = c("Ligand Cluster", "Ligand", "Receptor", "Receptor Cluster"), expand = c(.05, .05)) +
          ggplot2::scale_fill_manual(values=colp,limits=names(colp),name='Upregulated')+
          ggplot2::ggtitle(pltname) +
          ggplot2::theme(text = element_text(size = 8))+
          ggplot2::theme_minimal()
          )
  }
  else{
      print(paste0('Gene->',target,' Not Found'))

  }
}
