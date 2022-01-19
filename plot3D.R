#!env Rscript

library(rgl)
library(bio3d)


prot3D.init <- function (new.device = FALSE, bg = "gray", width = 640) { 
    if( new.device | rgl.cur() == 0 ) {
        rgl.open()
        par3d(windowRect = 50 + c( 0, 0, width, width ) )
        rgl.bg(color = bg )
    }
    rgl.clear(type = c("shapes", "bboxdeco"))
    rgl.viewpoint(theta = 15, phi = 20, zoom = 0.7)
}

prot3D.loadEdges <- function (fname) {
    # convention : vertex2 x1 y1 z1 vertex1 x2 y2 z2
    edges <- read.table(fname,header=T)
    return(edges)
}

plot3D.uniPoints <- function (edges) {
    pts <- data.frame(Node=c(edges$Node1, edges$Node2),
                         x = c(edges$x1, edges$x2),
                         y = c(edges$y1, edges$y2),
                         z = c(edges$z1, edges$z2))
    return(pts)
}

prot3D.plot <- function (edges,color) {
    pts <- plot3D.uniPoints(edges)
    pts <- pts[ !duplicated(pts$Node), ]
    rgl.spheres(pts$x,pts$y,pts$z,r=0.4,col=color)
    # rgl.texts(pts$x1,pts$y1,pts$z1,pts$Node1)
    segments3d(cbind(edges$x1,edges$x2),cbind(edges$y1,edges$y2),cbind(edges$z1,edges$z2),col=color)
}

# 2CPK_graphs.tsv
# 3LCK_graphs.tsv
plotSub <- function () {
    prot3D.init()
    edges <- prot3D.loadEdges("2CPK_graphs.tsv")
    prot3D.plot(edges,color="blue")
    edges <- prot3D.loadEdges("3LCK_graphs.tsv")
    prot3D.plot(edges,color="green")
}


plotSub()
