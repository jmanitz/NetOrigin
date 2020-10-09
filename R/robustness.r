
######################### robustness for origin detection methods ###############

# @export
robustness <- function(x, ...) UseMethod("robustness")

#' run robustness analysis for a source estimate by subsampling individual events.
#'
#' @details We create subsamples of individual events and their magnitude using a sampling proportion p in [0, 1]. After aggregating the data, we apply the source estimation approach. Using this result, we deduce the relative frequency of how often the source estimate obtained with the complete data set can be recovered by source estimation based on the subsample. Thus, the estimate robustness is assessed by the proportion of estimate recovery.
#' 
#' @param x \code{data.frame}, dataset with individual events and their magnitude, to be passed to \code{\link{aggr_data}}
#' @param type character, specifying the method, \code{'edm'}, \code{'backtracking'} and \code{'centrality'} are available.
#' @param prop numeric, value between zero and one, proportion of events to be sampled
#' @param n numeric, number of resamplings
#' @param ... parameters to be passed to origin methods \code{\link{origin_edm}}, \code{\link{origin_backtracking}} or \code{\link{origin_centrality}}
#' @return \code{data.frame} with columns
#'   \itemize{
#'         \item \code{est} origin estimated when all data is evaluated
#'         \item \code{rob} estimate uncertainty, computed as the proportion of resamplings when origin estimate was recovered
#'   }
#'
#' @import igraph
#' @examples
#' # generate random delay data
#' data(ptnAth)
#' require(igraph)
#' dat <- data.frame(node  = sample(size = 500, make.names(V(ptnAth)$name), replace = TRUE),
#'                   time  = sample(size = 500, 1:10, replace = TRUE),
#'                   delay = rexp(500, rate=10))
#' # compute effective distance
#' net <- igraph::as_adjacency_matrix(ptnAth, sparse=FALSE)
#' p <- net/rowSums(net)
#' eff <- eff_dist(p)
#' colnames(eff) <- paste('x.',colnames(eff),sep='')
#'
#' # run robustness analysis
#' r5 <- robustness(x=dat, type='edm', prop=0.5, n=10, distance=eff)
#' summary(r5)
#' plot(r5)
#'
#' # compare results
#' r9 <- robustness(x=dat, type='edm', prop=0.9, n=10, distance=eff)
#' plot(r9, add=TRUE, col='gray')
#'
#' @rdname robustness
#' @seealso \code{\link{robustness-methods}}
#' @export
robustness <- function(x, type=c('edm', 'backtracking', 'centrality'), prop, n=100, ...){
    N <- dim(x)[1]
    
    # run source detection on all data
    datA <- aggr_data(x,cumsum=TRUE)#, from='06:00')
    res <- apply(datA, FUN=origin, type=type, MARGIN=1, ...)
    est <- unlist(lapply(res, function(y) rownames(y[[2]])[y$est]))
    nt <- length(est) # number of time points

    # initialize return vector counting the number of estimate replica 
    rob <- integer(nt)
    
    # check robustness of source estimation
    cat('Run robustness analysis for source estimate: \n')
    for(i in 1:n){
        # buffer bar
        cat('.')
        # sample data
        repeat{
            sam <- sample(size=round(prop*N), 1:N, replace=FALSE)
            dat.sam <- aggr_data(x[sam,], cumsum=TRUE)#, from='06:00') 
            if(dim(dat.sam)[1] == nt) break # not neccessarily all times sampled
        }
        # run source detection
        res.sam <- apply(dat.sam, FUN=origin, type=type, MARGIN=1, ...)
        est.sam <- unlist(lapply(res.sam, function(x) rownames(x[[2]])[x$est]))
        # update number of estimate replica
        rob <- rob + as.numeric( est == est.sam )
    }

    cat('\n')
    ret <- data.frame(est=est, rob=rob/n)
    class(ret) <- 'robustness'
    return(ret)
}

#' @name robustness-methods
#' @aliases print.robustness
#' @aliases summary
#' @aliases plot
#' 
#' @title methods for robustness estimation objects of class \code{robustness}
#'
#' @description \code{print} produces an output for objects of class \code{robustness}
#'
#' @param x \code{data.frame} obtained by \code{\link{robustness}}, robustness estimation object for source estimation from function \code{\link{robustness}}
#' @param add logical specifying whether this should be added to another robustness plot
#' @param ... further arguments passed to the default \code{print} method
#'
#' @seealso \code{\link{robustness}}
#' @rdname robustness-methods
#' @export
print.robustness <- function(x, ...){
    ret <- data.frame(time = attr(x,'row.names'), estimate = x$est, robustness = x$rob)
    print(ret, ...)
}

#' \code{summary} produces an object summary for objects of class \code{robustness}
#'
#' @param object object of class \code{\link{origin}}, origin estimation object from function \code{origin_xxx}; passed to \code{x}
#' 
#' @rdname robustness-methods
#' @export
summary.robustness <- function(object, x = object, ...){
    ret <- data.frame(time = attr(x,'row.names'), estimate = x$est, robustness = x$rob)
    summary(ret, ...)
}

#' \code{plot} produces a time series plot of the \code{robustness} estimate object
#' 
#' @param y not used; default \code{NULL}
#' 
#' @rdname robustness-methods
#' @importFrom graphics abline legend lines points text
#' @export
plot.robustness <- function(x, y=NULL, add = FALSE, ...){
    # lines plot
    if(add){ 
       lines(x$rob, ...) 
    }else{
    plot(x$rob, type='l', ylim=c(-0.15,1.0), ylab='proportion of estimate recovery', ...)   
    # add estimates
    est <- levels(x$est)
    pos <- match(est, x$est)
    abline(v=pos, col='grey', lty=3)
    text(pos, -0.1, labels=est, cex=0.7, srt=45)
    abline(h=0)
    legend('bottomright',legend='source estimate', bty='n')
}}

