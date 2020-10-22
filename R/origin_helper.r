origin <- function(events, ...) UseMethod("origin")
#' Origin Estimation for Propagation Processes on Complex Networks
#'
#' This is the main function for origin estimation for propagation processes on complex networks. Different methods are available: effective distance median (\code{'edm'}), recursive backtracking (\code{'backtracking'}), and centrality-based source estimation (\code{'centrality'}).
#' For details on the methodological background, we refer to the corresponding publications.
#'
#' @rdname origin
#' @author Juliane Manitz with contributions by Jonas Harbering
#'
#' @references \itemize{
#'   \item Comin, C. H. and da Fontoura Costa, L. (2011). Identifying the starting point of a spreading process in complex networks. Physical Review E, 84. <DOI: 10.1103/PhysRevE.84.056105>
#'   \item Manitz, J., J. Harbering, M. Schmidt, T. Kneib, and A. Schoebel (2017): Source Estimation for Propagation Processes on Complex Networks with an Application to Delays in Public Transportation Systems. Journal of Royal Statistical Society C (Applied Statistics), 66: 521-536.
#'   \item Manitz, J. (2014). Statistical Inference for Propagation Processes on Complex Networks. Ph.D. thesis, Georg-August-University Goettingen. Verlag Dr.~Hut, ISBN 978-3-8439-1668-4. Available online: \url{http://ediss.uni-goettingen.de/handle/11858/00-1735-0000-0022-5F38-B}.
#'   \item Manitz, J., Kneib, T., Schlather, M., Helbing, D. and Brockmann, D. (2014). Origin detection during food-borne disease outbreaks - a case study of the 2011 EHEC/HUS outbreak in Germany. PLoS Currents Outbreaks, 1. <DOI: 10.1371/currents.outbreaks.f3fdeb08c5b9de7c09ed9cbcef5f01f2>
#' }
#'
#' @param events numeric vector of event counts at a specific time point; if type is 'bayesian', 'events' is a matrix, number of nodes x time points; entries represent number of cases
#' @param type character specifying the method, \code{'edm'}, \code{'backtracking'}, \code{'centrality'} and \code{'bayesian'} are available.
#' @param ... parameters to be passed to origin methods \code{\link{origin_edm}}, \code{\link{origin_backtracking}}, \code{\link{origin_centrality}} or \code{\link{origin_centrality}}
#'
#' @family origin-est
#' @export
origin <- function(events, type=c('edm', 'backtracking', 'centrality', 'bayesian'), ...){
    type <- match.arg(type)
    switch(type,
           edm = origin_edm(events = events, ...),
           backtracking = origin_backtracking(events = events, ...),
           centrality = origin_centrality(events = events, ...),
           bayesian = origin_bayesian(events = events, ...))
}


#################### standard methods for origin objects ######################

# add generic 
#print <- function(x) UseMethod("print")
#' @name origin-methods
# #' @aliases print.origin
# #' @aliases summary
# #' @aliases plot
# #' @aliases performance
#'
#' @title methods for origin estimation objects of class \code{origin}
#' 
#' @description \code{print} produces an output for objects of class \code{origin}.
#'
#' @param x object of class \code{\link{origin}}, origin estimation object from function \code{origin_xxx}
#'
#' @rdname origin-methods
#' @seealso \code{\link{origin}} \code{\link{plot_performance}}
#' @export
print.origin <- function(x, ...){
  if(!is.na(x$est)){ 
     switch(x$type, edm = cat('Effective distance median origin estimation:\n\n'),
		    backtracking = cat('Backtracking origin estimation:\n\n'),
                    centrality = cat('Centrality-based origin estimation:\n\n'))
    cat(paste('estimated node of origin', x$est))    
    if(!is.null(rownames(x[[2]]))) cat(paste(':',rownames(x[[2]])[x$est],'\n'))
    else cat('\n')
  }else{
    cat('source estimation not available\n')
  }
  return(invisible(x))
}

# add generic 
#summary <- function(x, ...) UseMethod("summary")
#' \code{summary} produces an object summary for objects of class \code{origin}.
#'
#' @param object object of class \code{\link{origin}}, origin estimation object from function \code{origin_xxx}; passed to \code{x}
#' 
#' @rdname origin-methods
#' @export
summary.origin <- function(object, x = object, ...){
    print(x)
    cat('\nauxiliary variables:\n')
    print(summary(x$aux))
    return(invisible(x))
}

#' \code{plot} generates an illustration of an origin object using the variable to be optimized.
#'
#' @param y character specifying the variable being plotted at the y-axis; options are \code{'id'} for node identifier (default), \code{'mdist'} for mean distance (only available for \code{\link{origin_edm}}) or \code{'wvar'} for weighted variance (only available for \code{\link{origin_edm}})
#' @param start numeric, giving the node of the true origin
#'
#' @rdname origin-methods
#' 
#' @importFrom graphics abline axis legend par plot points rect text title
#' @export
plot.origin <- function(x, y='id', start, ...){
    # extract estimation result
    k0.hat <- x$est
    res <- x$aux
    K <- nrow(res) # number of nodes
    # extract node names
    node.names <- if( is.null(rownames(res)) ) 1:K else rownames(res)
    # convert start as numeric
    if(is.character(start) && !is.null(rownames(res))){
	start <- match(start,rownames(res))
    }
    # define point size proportional to event magnitude 
    x <- sqrt(res$events)
    px <- x/max(x)*3+0.5 # point size propotional to events observed

    # specify what should be plotted
    y <- match.arg(y, c('id', 'wvar', 'mdist'))
    # plot: aux[,3] (cent,wmean) ~ id scatterplot
    if(y == 'id'){
      xy <- res[,c(3,1)]
    }
    # plot: weighted mean - unweighted mean effective distance scatterplot
    if(y == 'mdist'){
      xy <- res[,c('wmean','mdist')]
    }
    # plot: mean - variance scatterplot
    if(y == 'wvar'){
      xy <- res[,c('wmean','wvar')]
    }

    ### Plot
    plot(xy, las=1, bty='l', col='darkgrey', pch=20, cex=px, ...)
    points(xy[k0.hat,], col='limegreen', pch=13, cex=1.5, lwd=1.5)

    # mark true and estimated origin
    if(is.null(start)){      # if true origin is NOT known
       legend('bottomright', pch=13, lwd=1.5, lty=NA, bty='n',
              legend=paste('estimation:', node.names[k0.hat]), col=c('limegreen'))
    }else{                   # if true origin is known
       name.start <- ifelse(is.na(start), NA, node.names[start])
       points(xy[start,-1], col='indianred1', pch=13, cex=1.5, lwd=1.5)
       legend('bottomright', bty='n',
              col=c('indianred1','limegreen'), pch=13, lwd=1.5, lty=NA, 
              legend=c(paste('true origin:',node.names[start]),
                       paste('estimation:', node.names[k0.hat])))
    }
}

#### evaluation method for origin objects

# add generic for evaluation 
#' generic method for performance evaluation
#' @param x object
#' @param ... further arguments
#' @seealso \code{\link{origin-methods}} \code{\link{plot_performance}}
#' @export
performance <- function(x, ...) UseMethod("performance")


#' \code{performance} evaluates an object of class \code{origin} and returns a \code{data.frame} identifying correct estimation, and computing rank and distance of correct detection.
#' 
#' @param graph \code{\link{igraph}} object specifying the underlying network graph with attribute 'length' on edges for calculation of distance to the correct origin
#' @param ... further arguments to be passed to default \code{plot} function
#'
#' @return \code{performance.origin} returns a \code{data.frame} with variables
#'   \itemize{
#'     \item \code{origin = start} representing the true origin, 
#'     \item \code{est} the estimated node of origin, 
#'     \item \code{hitt} logical indicating whether origin estimation is correct or not, 
#'     \item \code{rank} rank of correct detection, 
#'     \item \code{spj} number of segments from estimated origin to true origin (requires an \code{\link{igraph}} object), 
#'     \item \code{dist} distance along the shortest path from estimated origin to true origin (\code{\link{igraph}} edge attribute \code{length})
#'   }
#'
#' @examples
#' data(ptnGoe)
#' data(delayGoe)
#'
#' res <- origin(events=delayGoe[10,-c(1:2)], type='centrality', graph=ptnGoe)
#' res
#'
#' summary(res)
#' performance(res, start=1, graph=ptnGoe)
#'
#' @import igraph
#' @rdname origin-methods
#' @export
performance.origin <-  function(x, start, graph=NULL, ...){

    # extract estimation result
    k0.hat <- x$est
    aux <- x$aux
    K <- nrow(aux) # number of nodes
    # extract node names
    node.names <- if( is.null(rownames(aux)) ) 1:K else rownames(aux)
    # convert start as numeric
    if(is.character(start) && !is.null(rownames(aux))){
	start <- match(start,rownames(aux))
    }

    ### evaluation measures
    ret <- data.frame(start = node.names[start], est = NA, 
                      hitt = 'missing', rank = NA, spj = NA, dist = NA)
#                      r.err = NA, v.coef = NA)
    # no source detection to evaluate
    if(is.na(k0.hat)) return(ret)
    else ret$est = node.names[k0.hat]
    # correct source detection
    if(start == k0.hat){
       ret$hitt <- 'TRUE'
       ret$rank <- 1
       ret$spj <- ret$dist <- 0
#       ret$r.err <- ret$v.coef <- 0
    }else{
       ret$hitt <- 'FALSE'
       # rank of correct detection
       ret$rank <- switch(x$type,
           # effective distance median: minimize weighted mean
           edm = rank(aux$wmean, na.last=TRUE, ties.method='min')[start],
           # backtracking: maximize backtracking source counts
	   backtracking = rank(-aux$bcount, na.last=TRUE, ties.method='min')[start],
           # centrality-based: maximize centrality
           centrality = rank(-aux$cent, na.last=TRUE, ties.method='min')[start])
       # distance to correct detection
       if(!is.null(graph)){ 
         sp <- shortest_paths(graph, from=node.names[start], to=node.names[k0.hat], output='epath')$epath[[1]]
         ret$spj <- length(sp)
         ret$dist <- sum(sp$length) # grabs edge attribute 'length'
       }else{
         warning("Distance to correct detection cannot be computed, because 'graph' is missing.")
       }
#      # relative error -> not available for centrality
#      ret$r.err <- abs( res$wmean[k0.hat] - res$wmean[start] ) / res$wmean[start]
#      # variation coefficient -> not available for centrality
#      ret$v.coef <- sqrt( res$wvar[k0.hat] ) / res$wmean[k0.hat]
    }
    return(ret)    
}

###################### further evaluation methods

#' A plot method combining a time series of performance results.
#'
#' @param x \code{data.frame} obtained by combined results from \code{\link{performance.origin}} with variables \code{X1} for time point, \code{start} for true origin, \code{est} for estimated origin, and performance variables
#' @param var character, variable to be plotted, \code{\link{performance.origin}} returns \code{rank}, \code{spj}, and \code{dist}, default is \code{'rank'}
#' @param add logical, should be added to another performance plot
#' @param offset \code{POSIXct}, starting time of spreading 
#' @param log logical, should y-axis be logarithmized?
#' @param col numeric or character, color of lines
#' @param ylim numeric vector, range of y axis
#' @param text.padding a numeric value specifying the factor for the text position relative to the y values
#' @param ... further graphical parameters passed to default \code{plot} function
#'
#' @import igraph 
#' @examples
#' \dontrun{ 
#' ### delays on Goettingen bus network
#' # compute effective distance
#' data(ptnGoe)
#' goenet <- igraph::as_adjacency_matrix(ptnGoe, sparse=FALSE)
#' p <- goenet/rowSums(goenet)
#' eff <- eff_dist(p)
#' # apply source estimation
#' data(delayGoe)
#' if (requireNamespace("aplyr", quietly = TRUE)) {
#'    res <- alply(.data=delayGoe[11:20,-c(1:2)], .margins=1, .fun=origin_edm, 
#'                 distance=eff, silent=TRUE, .progress='text')
#'    perfGoe <- ldply(Map(performance, x = res, start = 2, list(graph = ptnGoe)))
#'    # performance plots
#'    plot_performance(perfGoe, var='rank', ylab='rank of correct detection', text.padding=0.5)
#'    plot_performance(perfGoe, var='dist', ylab='distance to correct detection')
#' }
#' 
#' ### delays on Athens metro network
#' # compute effective distance
#' data(ptnAth)
#' athnet <- igraph::as_adjacency_matrix(ptnAth, sparse=FALSE)
#' p <- athnet/rowSums(athnet)
#' eff <- eff_dist(p)
#' # apply source estimation
#' data(delayAth)
#' if (requireNamespace("aplyr", quietly = TRUE)) {
#'    res <- alply(.data=delayAth[11:20,-c(1:2)], .margins=1, .fun=origin_edm, 
#'              distance=eff, silent=TRUE, .progress='text')
#'    perfAth <- ldply(Map(performance, x = res, start = as.list(delayAth$k0),
#'                      list(graph = ptnAth)))
#'    # performance plots
#'    plot_performance(perfAth, var='rank', ylab='rank of correct detection',text.padding=0.5)
#'    plot_performance(perfAth, var='dist', ylab='distance to correct detection')
#' }
#' }
#' @importFrom graphics abline legend lines points text
#' @export
plot_performance <- function(x, var='rank', add=FALSE, offset=NULL, log=FALSE, col=1, ylim=NULL, text.padding = 0.9, ...){
    # specify variable to be plotted
    time <- try(as.POSIXct(x[,1]), silent=TRUE)
    if(inherits(time,'try-error')) time <- x[,1]
    y <- x[,paste(var)]
#print(cbind(time,y))

    # plotting parameter specification
    if(is.null(ylim)) ylim <- range(y, na.rm=TRUE)+c(-0.5,0)
    if(!add){
      if(log) plot(x=time, y=y, ylim=ylim, type='n', las=1, log='y', ...)
      else plot(x=time, y=y, ylim=ylim, type='n', las=1, ...)
    }
    
    # rank ts
    lines(x=time, y=y, lwd=2, col=col, ...)
    # time of start spreading
    if(!is.null(offset)) abline(v=offset, col=2)
    # add source estimates
    est <- levels(factor(x$est))
    pos <- match(est, x$est)   
#print(pos)
    points(x=time[pos], y=y[pos], col=col)
#    abline(v=time[pos], col=col, lty=3)
#    abline(h=1, col='grey')
#    if(log) text(time[pos], y[pos]*0.65, labels=est, cex=0.7, srt=45, col=col)
    text(time[pos], y[pos]+text.padding, labels=est, cex=0.7, srt=45, col=col)
}



