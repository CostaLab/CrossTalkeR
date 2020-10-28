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

  igraph::V(ce)$size<-5
  igraph::V(ce)$color<-col[as.numeric(igraph::V(ce))]
  igraph::V(ce)$label.color<-'black'
  igraph::V(ce)$label.cex<-1.5
  igraph::E(ce)$width <- igraph::E(ce)$inter*100
  we <- oce::rescale(igraph::E(ce)$weight,xlow=(-emax),xhigh=emax,rlow=1,rhigh=200,clip=TRUE)
  we <- round(we,0)
  igraph::E(ce)$arrow.width<- ifelse(igraph::E(ce)$width > 20, 5,igraph::E(ce)$width)
  igraph::E(ce)$label.color<-'black'
  igraph::E(ce)$label.cex<-1
  igraph::E(ce)$color <- col_pallet[we]
  aux <- oce::rescale(abs(igraph::E(ce)$weight), rlow=1,rhigh=100)
  aux <- ifelse(aux < high & aux > low, 0, aux)
  igraph::E(ce)$color <- scales::alpha(igraph::E(ce)$color, aux)
  igraph::E(ce)$width <- ifelse(igraph::E(ce)$width < 150 & igraph::E(ce)$width > 50, 0, igraph::E(ce)$width)#abs(E(ce)$weight)
  if(sum(edge.start[,2]==edge.start[,1])!=0){
    igraph::E(ce)$loop.angle[which(edge.start[,2]==edge.start[,1])]<-loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
  }
  #vertex.label=NA
  plot(ce,edge.curved=0.5,vertex.label=NA,vertex.shape='circle',vertex.lwd=1,layout=coords_scale,margin=0.0,xlim=c(-2,2),ylim=c(-2,2),edge.lty=igraph::E(ce)$lty,loop.angle=igraph::E(ce)$loop.angle,edge.label=NA,main=plt_name)
  e_wid_sp <- c(min(igraph::E(ce)$width),(min(igraph::E(ce)$width)+max(igraph::E(ce)$width))/2,max(igraph::E(ce)$width))
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
    netdiffuseR::drawColorKey(seq(0,200),tick.marks = c(0,200), color.palette=col_pallet,labels=c(-round(emax,4),round(emax,4)),nlevels=200,main="LR Interactions",pos=2,key.pos = c(0.97,1.0, 0.0, 0.4),border='transparent')
  }
}







