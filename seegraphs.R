#!/usr/bin/env Rscript

load_data <- function(filename) {
    return(read.delim(filename, as.is=T))
}

unique_points <- function(graph) {
    pts <- data.frame(Node=c(graph$Node1, graph$Node2),
                      x = c(graph$x1, graph$x2), 
                      y = c(graph$y1, graph$y2), 
                      z = c(graph$z1, graph$z2))
    return( pts[ !duplicated(pts), ] )
}

#seg_params <- function(graph, ...) {
#    list(x=cbind(graph$x1, graph$x2),
#         y=cbind(graph$y1, graph$y2),
#         z=cbind(graph$z1, graph$z2), ...)
#}

seg_params <- function(graph, ...) {
    list(x=as.vector(t(graph[,c("x1","x2")])),
           y=as.vector(t(graph[,c("y1","y2")])),
           z=as.vector(t(graph[,c("z1","z2")])))
}

seegraph <- function(filename, sph.radius=1, #sph_par=list(),
                     seg_par=list(), save=NULL) {
    require(rgl)
    edges <- load_data(filename)
    nodes <- unique_points(edges)
    #sph_par <- c(list(x=nodes$x, y=nodes$y, z=nodes$z), sph_par)
    sph_par <- list(x=nodes$x, y=nodes$y, z=nodes$z, radius=sph.radius)
    do.call(rgl.spheres, sph_par)
    #rgl.texts(nodes$x, nodes$y, nodes$z, text=nodes$Node, adj=c(1,1))
    text3d(nodes$x+sph.radius, nodes$y+sph.radius, nodes$z+sph.radius,
           text=nodes$Node, adj=c(1,1))
    do.call(segments3d, c(seg_params(edges), seg_par))
    #if( !is.null(save) )
    #    rgl.postscript(filename=save, fmt="pdf")
    return(list(edges=edges, nodes=nodes))
}

# rotation taken from the PDBfold ssm structural alignment of 2CPK and 3LCK
rot_mat <- matrix(c(
  0.127,   0.812,   0.569, 
  0.073,  -0.580,   0.811, 
  0.989,  -0.061,  -0.133), nrow=3, byrow=T) 

tr_mat <- c(-73.389, -38.434, -7.648)

rotatetranslate <- function(line) {
    coord <- as.numeric(line[2:4])
    return( rot_mat %*% coord + tr_mat )
}

if( !interactive() ) {
    CPK <- seegraph("2CPK_graphs_01.tsv", list(col='lightblue'), list(col='grey20'),
                save="2CPK.png") 
    LCK <- seegraph("3LCK_graphs_01.tsv", list(col='lightgreen'),
                list(col='grey20'), save="3LCK.png")
}
#LCK_p_adj <- t(apply(LCK_p, 1, rotatetranslate))
#colnames(LCK_p_adj) <- c('x', 'y', 'z')
 
#rgl.spheres(LCK_p_adj[,1], LCK_p_adj[,2], LCK_p_adj[,3], col='lightgreen')
#do.call(segments3d, seg_params(LCK_p_adj, col="lightblue"))
 
