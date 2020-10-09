#' analyze network characteristics
#' 
#' @param graph \code{\link{igraph}} object, network graph representing the public transportation network, vertices represent stations, which are linked by an edge if there is a direct transfer between them
#' @return \code{'data.frame':} 1 obs. of 7 variables: 
#'   \itemize{
#'      \item \code{vcount} number of nodes, 
#'      \item \code{ecount} number of edges, 
#'      \item \code{density} network graph density, 
#'      \item \code{av_deg} average degree, 
#'      \item \code{av_cent} average unit betweenness,
#'      \item \code{av_spl} average shortest path length, 
#'      \item \code{diam} diameter, and 
#'      \item \code{trans} transitivity.
#'      \item \code{cluster} number of clusters
#' }
#' 
#' @references Details to the computation and interpretation can be found in:\itemize{
#'  \item Kolaczyk, E. D. (2009). Statistical analysis of network data: methods and models. Springer series in statistics. Springer. <DOI: 10.1007/978-0-387-88146-1>
#'  \item Manitz, J. (2014): Statistical Inference for Propagation Processes on Complex Networks. Ph.D. thesis, Georg-August-University Goettingen. Verlag Dr.~Hut, ISBN 978-3-8439-1668-4. Available online: \url{http://ediss.uni-goettingen.de/handle/11858/00-1735-0000-0022-5F38-B}.
#' }
#'
#' @examples
#' data(ptnAth)
#' analyze_net(ptnAth)
#'
#' data(ptnGoe)
#' analyze_net(ptnGoe)
#' 
#' @import igraph
#' @family network helper
#' @export
analyze_net <- function(graph){
    K <- vcount(graph)
    # compute graph characteristics
    res <-  data.frame(vcount = K,
                       ecount = ecount(graph),
                       # graph density
                       density = edge_density(graph),
                       # average degree
                       av_deg = mean(degree(graph)),
                       # av unit betweenness
                       av_cent = mean(betweenness(graph))/((K-1)*(K-2)/2),
                       # av shortest path length
                       av_spl = mean_distance(graph),
                       # diameter
                       diam = diameter(graph),
                       # transitivity
                       trans = transitivity(graph, type='global'),
                       # number of clusters
		       cluster = clusters(graph)$no)
#    power.law.fit(degree(graph))
    # output 
    return(res)
}

#' A plot method for public transportation networks (PTNs).
#' 
#' @param graph \code{\link{igraph}} object, network graph representing the public transportation network, vetrices represent stations, which are linked by an edge if there is a direct transfer between them
#' @param color.coding numeric vector with length equal to the number of network nodes
#' @param color.scheme character vector of length 5 indicating the \code{vertex.color}, default is \code{rev(sequential_hcl(5))}
#' @param legend logical indicating whether legend for color-coding should be added or not.
#' @param ... further arguments to be passed to \code{\link{plot.igraph}}
#'
#' @examples
#' data(ptnAth)
#' plot_ptn(ptnAth)
#'
#' data(ptnGoe)
#' plot_ptn(ptnGoe)
#' 
#' @import igraph
#' @import colorspace
#' @family network helper
#' @export
plot_ptn <- function(graph, color.coding = NULL, color.scheme = rev(sequential_hcl(5)), legend = FALSE, ...){

    # color.coding
    if(is.null(color.coding)){
        vertex.color <- color.scheme[4]
    }else{
        r <- cut(color.coding, breaks=5)
        vertex.color <- color.scheme[r]
        legend <- TRUE
    }
    # layout
    ltr <- layout_with_kk(graph) 

    # plot
    plot.igraph(graph, xlim=c(-2,2), ylim=c(-1.75,1.75), asp=0.6, margin=-1.05, layout=ltr,
         vertex.size=2, vertex.label = NA, vertex.color = vertex.color, #vertex.shape = 'circle', 
         edge.curved = FALSE, ...) #edge.lty=3, edge.width=2,

    if(legend){
        legend(-1.25,-0.55, legend=levels(r), col=color.scheme, pch=16, bty='n')
    }
    
}


