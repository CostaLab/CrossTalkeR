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
plotLR_diff1 <- function(ce,col,plt_name, coords,emax=NA,leg=FALSE,low=25,high=75){
  if(is.na(emax)){
    emax <-max(abs(igraph::E(ce)$mean))
    emin <-min(abs(igraph::E(ce)$mean))
  }
  else{
    emin = emax
  }
  col_pallet <- pals::coolwarm(3)
  col_pallet <- colorRampPalette(col_pallet)(201)
  edge.start <- ends(ce, es=E(ce), names=FALSE)
  if(nrow(coords)!=1){
    coords_scale=scale(coords)
  }else{
    coords_scale<-coords
  }
  loop.angle<-ifelse(coords_scale[V(ce),1]>0,-atan(coords_scale[V(ce),2]/coords_scale[V(ce),1]),pi-atan(coords_scale[V(ce),2]/coords_scale[V(ce),1]))
  loop.angle<-0
  igraph::V(ce)$size<-10
  igraph::V(ce)$color<-col[as.numeric(igraph::V(ce))]
  igraph::V(ce)$label.color<-'black'
  igraph::V(ce)$label.cex<-1.5
  we <- round(oce::rescale(igraph::E(ce)$weight,xlow=(-emax),xhigh=emax,rlow=1,rhigh=200,clip=TRUE),0)
  igraph::E(ce)$label.color<-'black'
  igraph::E(ce)$label.cex<-1
  igraph::E(ce)$color <- col_pallet[we]
  aux <- E(ce)$inter
  aux <- ifelse(aux < high & aux > low, 0, aux)
  igraph::E(ce)$arrow.width<- ifelse(aux < high & aux > low, 0, (aux/25)+1.5)
  igraph::E(ce)$color <- scales::alpha(igraph::E(ce)$color, aux)
  igraph::E(ce)$width <- aux#abs(E(ce)$weight)
  if(sum(edge.start[,2]==edge.start[,1])!=0){
    igraph::E(ce)$loop.angle[edge.start[,2]==edge.start[,1]]<-loop.angle[edge.start[edge.start[,2]==edge.start[,1],1]]
  }
  #vertex.label=NA
  plot(ce,edge.curved=0.5,vertex.label=NA,vertex.shape='circle',vertex.lwd=1,layout=coords_scale,margin=0.0,xlim=c(-2,2),ylim=c(-2,2),edge.lty=igraph::E(ce)$lty,loop.angle=igraph::E(ce)$loop.angle,edge.label=NA,main=plt_name)
  e_wid_sp <- c(min(igraph::E(ce)$inter[E(ce)$inter!=0]),(min(igraph::E(ce)$inter[E(ce)$inter!=0])+max(igraph::E(ce)$inter))/2,max(igraph::E(ce)$inter))
  legend('topleft', legend =round(e_wid_sp,1),col='gray',title='Percentage of the interactions',
         pch = NA, bty = "n", cex =1, lwd=e_wid_sp,lty=c(1,1,1),
         text.col = "black", horiz = FALSE)

  x = coords_scale[,1]*1.2
  y = coords_scale[,2]*1.2
  angle = ifelse(atan(-(coords_scale[,1]/coords_scale[,2]))*(180/pi) < 0,  90 + atan(-(coords_scale[,1]/coords_scale[,2]))*(180/pi), 270+atan(-coords_scale[,1]/coords_scale[,2])*(180/pi))
  for (i in 1:length(x)) {
    text(x=x[i], y=y[i], labels=igraph::V(ce)$name[i], adj=NULL, pos=NULL, cex=0.8, col="black", xpd=T)#srt=angle[i], xpd=T)
  }
  if(leg){
    netdiffuseR::drawColorKey(seq(0,200),tick.marks = c(0,200), color.palette=col_pallet,labels=c(-round(emax,4),round(emax,4)),nlevels=200,main="LR Interactions",pos=2,key.pos = c(0.97,1.0, 0.0, 0.3),border='transparent')
  }
}







#'This function do a ggi plot and higest degree nodes
#'
#'@param graph graph
#'@param color cluster color
#'@import ggplot2
#'@import ggraph
#'@return R default plot
#'@importFrom tidyr %>%
#'@export
plot_ggi<-function(graph,color){
  deg <- igraph::degree(graph)
  deg1 <- deg[order(deg,decreasing = T)]
  root <-names(deg1)[1]
  v_names <- igraph::V(graph)$name
  v_names <- tibble::as_tibble(v_names)
  v_names <-v_names %>%
    tidyr::separate(value, c("genes", "cluster"), "_")
  igraph::V(graph)$genes <- v_names$genes
  igraph::V(graph)$cluster <- v_names$cluster
  ecolors <- rev(pals::coolwarm(25))[round(oce::rescale(E(graph)$MeanLR,xlow = min(E(graph)$MeanLR),rlow = 1, rhigh = 25))]
  igraph::E(graph)$colors <- ecolors
  cls <- names(color)
  names(color) <- NULL
  ewidth <- (abs(igraph::E(graph)$MeanLR) - mean(abs(igraph::E(graph)$MeanLR)))/sd(abs(igraph::E(graph)$MeanLR))
  ewidth <-  (ewidth-min(ewidth))/(max(ewidth)-min(ewidth))
  print(ggraph::ggraph(graph,layout='stress')+
          ggraph::geom_edge_link0(ggplot2::aes(edge_width=ewidth,color=igraph::E(graph)$MeanLR),alpha=ewidth)+
          ggraph::scale_edge_color_gradient2(low=pals::coolwarm(25)[1],high=pals::coolwarm(25)[25])+
          ggraph::geom_node_point(size=(deg/max(deg)*10),alpha=1,ggplot2::aes(color=V(graph)$cluster))+
          ggplot2::scale_colour_manual(values=color,labels=cls,name='Clusters')+
          ggraph::geom_node_text(ggplot2::aes(filter=deg>deg1[10],label=genes))+
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
  deg1 <- deg[order(deg,decreasing = T)]
  v_names <- igraph::V(graph)$name
  v_names <- tibble::as_tibble(v_names)
  v_names <-v_names %>%
    tidyr::separate(value, c("genes", "cluster"), "_")
  igraph::V(graph)$genes <- v_names$genes
  igraph::V(graph)$cluster <- v_names$cluster
  igraph::V(graph)$lcolors <-color[V(graph)$cluster]
  ecolors <- rev(pals::coolwarm(25))[round(oce::rescale(E(graph)$MeanLR,xlow = min(E(graph)$MeanLR),rlow = 1, rhigh = 25))]
  igraph::E(graph)$colors <- ecolors
  cls <- names(color)
  names(color) <- NULL
  ewidth <- (abs(igraph::E(graph)$MeanLR) - mean(abs(igraph::E(graph)$MeanLR)))/sd(abs(igraph::E(graph)$MeanLR))
  ewidth <-  (ewidth-min(ewidth))/(max(ewidth)-min(ewidth))
  print(ggraph::ggraph(graph,layout='stress')+
          ggraph::geom_edge_link0(ggplot2::aes(edge_width=ewidth,color=igraph::E(graph)$MeanLR),alpha=ewidth)+
          ggraph::scale_edge_color_gradient2(low=pals::coolwarm(25)[1],high=pals::coolwarm(25)[25])+
          ggraph::geom_node_point(size=(deg/max(deg)*10),alpha=1,ggplot2::aes(color=V(graph)$cluster))+
          ggplot2::scale_colour_manual(values=color,labels=cls,name='Clusters')+
          ggraph::geom_node_text(ggplot2::aes(filter=igraph::V(graph) %in% igraph::articulation.points(graph),label=genes))+
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
#'@import ggalluvial
#'@importFrom ggalluvial StatStratum
#'@return R default plot
#'@importFrom tidyr %>%
#'@export
plot_sankey<-function(LRObj_tbl,target,Ligand.Cluster, Receptor.Cluster,pltname){
  data <- LRObj_tbl[grep(paste(target,"_", sep=""), LRObj_tbl$allpair), ]
  data <- data[grep(paste(Ligand.Cluster, collapse='|'), data$Ligand.Cluster), ]
  data <- data[grep(paste(Receptor.Cluster, collapse='|'), data$Receptor.Cluster), ]
  data$Freq <- 1
  if(dim(data)[1] >= 1){
    print(ggplot2::ggplot(
      top_n(data, ifelse(dim(data)[1] > 50, 50, dim(data)[1]), abs(MeanLR)),
      aes(
        y = Freq,
        axis1 = Ligand.Cluster, axis2 = Ligand, axis3 = Receptor, axis4 = Receptor.Cluster)
    ) +
      ggalluvial::geom_alluvium(aes(fill = MeanLR > 0), width = 1/12) +
      ggalluvial::geom_stratum(width = 1/12) +
      # geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2) +
      ggplot2::geom_label(stat = "stratum", ggplot2::aes(label = ggplot2::after_stat(stratum)), size = 2) +
      ggplot2::scale_x_discrete(limits = c("Ligand Cluster", "Ligand", "Receptor", "Receptor Cluster"), expand = c(.05, .05)) +
      ggplot2::scale_fill_manual(values=if(sum(table(MeanLR > 0)) == 0 ,pals::coolwarm(n=2)[1],pals::coolwarm(n=2)),name='Upregulated')+
      ggplot2::ggtitle(pltname) +
      ggplot2::theme(
        text = element_text(size = 8)
      ))
  }
  else{
      print(paste0('Gene->',target,' Not Found'))

  }
}

