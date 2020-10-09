# #' generic method for hpd
# #' @param object object of class \code{\link{origin}}
# #' @param ... further arguments
# #' @seealso \code{\link{hpd.origin}} \code{\link{origin-methods}} \code{\link{hpd-methods}}
# #' @export
hpd <- function(object, ...) UseMethod("hpd")

#' \code{hpd} computes the highest posterior density region for origin estimator of type \code{'gaussian'}
#' 
#' @param alpha numeric, type-I-error for HPD region, probability between 0 and 1
#' @param connected logical, indicating whether HPD region should be connected
# #' @param graph \code{\link{igraph}} object specifying the underlying network graph (required for computation of a connected HPD region)
#' 
#' @return \code{hpd.origin} returns an object of class \code{hpd}, which is a \code{list} with elements
#'   \itemize{
#'     \item \code{hpd} \code{list} with elements 
#'       \itemize{
#'          \item \code{nodes} nodes in HPD region
#'          \item \code{emp_prob} empirical probability covered
#'       }
#'     \item \code{hpd_connected} \code{list} with elements 
#'       \itemize{
#'          \item \code{nodes} nodes in HPD region and connecting nodes on paths between HPD region nodes
#'          \item \code{emp_prob} empirical probability covered
#'       }
#'     \item \code{origin} object of class \code{\link{origin}}
#'   } 
#'
#' @rdname origin-methods
#' @aliases hpd
#' @seealso \code{\link{hpd-methods}}
#' @export
hpd.origin <- function(object, alpha=0.05, connected=FALSE, graph=NULL){

    if(object$type!='gaussian') stop('HPD intervals only available for Gaussian source estimator')
    
    # posterior probability
    post <- object$aux$post

    # select HPD region nodes
    sel <- order(post)[!cumsum(sort(post)) < alpha]
    # compute empirical coverage probability
    emp.prob <- sum(post[sel])

    # construct connected set from HPD nodes
    if(connected){
        if(is.null(graph)){
#            warning("Connected HPD region can be only computed if 'graph' is given.")
            csel <- NA
            emp.prob.connect <- NA            
        }else{
            csel <- unique(unlist(shortest_paths(graph, from=sel, to=sel)$vpath))
            emp.prob.connect <- sum(post[csel])
        }
    }else{
        csel <- NULL
        emp.prob.connect <- emp.prob
    }

    # return object
    ret <- list(hpd = list(nodes = sel, emp_prob = emp.prob),
                hpd_connected = list(nodes = csel, emp_prob = emp.prob.connect),
                origin = object)
    class(ret) <- 'hpd'
    return(ret)
}

#' @name hpd-methods
#' @aliases print.hpd
#' @aliases summary.hpd
#' @aliases plot.hpd
#'
#' @title methods for highest posterior density region objects of class \code{hpd}
#' 
#' @description \code{print} produces an output for objects of class \code{hpd}.
#'
#' @param x object of class \code{\link{hpd}}
#'
#' @rdname hpd-methods
#' @seealso \code{\link{hpd}} \code{\link{hpd.origin}} \code{\link{origin-methods}}
#' @export
print.hpd <- function(x, ...){
    print(x$origin)
    
    cat('\nhighest posterior density region covers the nodes:', x$hpd$nodes, '\nwith empirical coverage ', round(x$hpd$emp_prob,3)*100, '%\n')

    if(!is.null(x$hpd_connected$nodes)){
        cat("Connected HPD region covers the nodes:", x$hpd$nodes, '\nwith empirical coverage ', round(x$hpd_connected$emp_prob,3)*100, '%\n')
        
    }
}

#' \code{summary} produces an object summary for objects of class \code{hpd}.
#'
#' @param object object of class \code{\link{hpd}}; passed to \code{x}
#'
#' @rdname hpd-methods
#' @export
summary.hpd <- function(object, x = object, ...){
    print(x)
    cat('\nauxiliary variables:\n')
    print(summary(x$origin$aux))
    return(invisible(x))
}

#' \code{plot} generates an illustration of the highest posterior density region for an origin object of type 'gaussian'
#'
#' @param graph \code{\link{igraph}} object specifying the underlying network graph 
#' @param ... further parameters to be passed to \code{\link{plot.igraph}}
#' @param y \code{NULL}
#'
#' @rdname hpd-methods
#' @export
plot.hpd <- function(x, y=NULL, graph, ...){

    # set alayout
    graph <- set_graph_attr(graph, "layout", layout_with_kk(graph))
    V(graph)$size <- 9; E(graph)$width <- 2
    V(graph)$frame.color <- NA

    # mark node set
    V(graph)$color <- 2
    V(graph)$color[x$hpd$nodes] <- 1
    V(graph)$color[x$hpd_connected$nodes] <- 4
    V(graph)$color[x$origin$est] <- 'brown2'

    # node size proportional to posterior 
  #  V(graph)$size <- sqrt(x$origin$aux$post)*50

    # plot
    plot(graph, mark.groups=x$hpd$nodes, mark.col="grey90", mark.border='grey20', ...)
    legend(-1,0.3, bty='n', pch=20, cex=1, col=c('brown2', categorical_pal(4)[c(1,4,2)]), legend=c('origin estimate','HPD region','connecting nodes', 'other nodes'))
}


