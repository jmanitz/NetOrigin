#' @include origin_helper.r

######################### origin detection methods #############################
################################################################################

#origin_edm <- function(x) UseMethod("origin_edm")
#' @title Origin Estimation for Propagation Processes on Complex Networks
#'
#' @description 
#' \code{origin_edm} for effective distance-median origin estimation (Manitz et al., 2016)
#'
#' @param distance numeric matrix specifying the distance matrix (for \code{type='edm'})
#' @return \code{origin_edm} returns an object of class \code{origin}, list with 
#' \itemize{
#'  \item \code{est} origin estimate
#'  \item \code{aux} \code{data.frame} with auxiliary variables 
#'    \itemize{
#'       \item \code{id} as node identifier, 
#'       \item \code{events} for event magnitude, 
#'       \item \code{wmean} for weighted mean, 
#'       \item \code{wvar} for weighted variance, and 
#'       \item \code{mdist} mean distance from a node to all other nodes.
#'    }
#'  \item \code{type = 'edm'} effective distance median origin estimation
#'  }
#' 
#' @examples
#' data(delayGoe)
#'
#' # compute effective distance
#' data(ptnGoe)
#' goenet <- igraph::as_adjacency_matrix(ptnGoe, sparse=FALSE)
#' p <- goenet/rowSums(goenet)
#' eff <- eff_dist(p)
#' # apply effective distance median source estimation
#' om <- origin(events=delayGoe[10,-c(1:2)], type='edm', distance=eff)
#' summary(om)
#' plot(om, 'mdist',start=1)
#' plot(om, 'wvar',start=1)
#' performance(om, start=1, graph=ptnGoe)
#' 
#' @rdname origin
#' @import Hmisc igraph
#' @export
origin_edm <- function(events, distance, silent=TRUE){    
    ### error handling
    if(!is.vector(events)) events <- as.vector(events)
    # NA handling in events
    nas <- which(is.na(events))
    if(length(nas)>0){
        events <- events[-nas]
#        distance <- distance[-nas,-nas]
        if(!silent) message(cat('missing values in data events at nodes:\n',nas))
    }
    # Inf handling in distance
    if(any(distance == Inf, na.rm=TRUE)){
        distance[which(distance == Inf)] <- max(distance[which(distance < Inf)]) * 1000
    }
    # check events/distance matching 
    if(is.null(names(events))){
        if(!silent) message('Note: events and distance cannot be matched, ensure they are in the same order')
    }else{
        order.events <- match(names(events),colnames(distance))
        # remove nodes with events that are not in the network
        nas <- which(is.na(order.events))
        if(length(nas) > 0){
          if(!silent) message(cat('Note: data events for nodes given that are not in the network:\n',names(events)[nas]))
          events <- events[-nas]
          order.events <- order.events[-nas]
	}
        # sort distance + remove unneeded information
        distance <- distance[order.events, order.events]
    }   

    ### convert input
    # define dimension
    K <- nrow(distance)    
    # convert event vector into dat matrix object
    dat <- matrix(as.numeric(events),ncol=1, nrow=K)

    ### compute source detection measures
    # weighted mean and variance
    wmean <- (distance %*% dat) / sum(dat)
    wvar <- apply(distance, MARGIN=2, FUN=var_wtd_mean_cochran, w=dat)
    
    # source detection 
    k0.hat <- which.min(wmean)
    if(length(k0.hat)==0) k0.hat <- NA

    # add mdist - mean distance from a node to all other nodes
#    mdist <- colMeans(distance)

    ### return
    ret <- list(est = k0.hat, 
                aux = data.frame(id=1:length(events), events=dat, wmean = wmean, wvar = wvar, mdist = colMeans(distance), row.names=names(events)),
                type = 'edm')
    class(ret) <- 'origin'
    return(ret)
}

# Computes the variance of a weighted mean
#' Computes the variance of a weighted mean following the definition by Cochran (1977; see Gatz and Smith, 1995)
#'
#' This is a helper method for weighted variance computation in \code{\link{origin_edm}}, which is the closest to the bootstrap.
#'
#' @param x numeric vector of values
#' @param w numeric vector of weights
#' @return numeric value of weighted variance
#' 
#' @references \itemize{
#'   \item Gatz, D. F., and Smith, L. (1995). The standard error of a weighted mean concentration-I. Bootstrapping vs other methods. Atmospheric Environment, 29(11), 1185-1193. <DOI: 10.1016/1352-2310(94)00210-C>
#'   \item Gatz, D. F., and Smith, L. (1995). The standard error of a weighted mean concentration-II. Estimating confidence intervals. Atmospheric Environment, 29(11), 1195-1200. <DOI: 10.1016/1352-2310(94)00209-4>
#'   \item \url{http://r.789695.n4.nabble.com/Problem-with-Weighted-Variance-in-Hmisc-td826437.html}
#' }
#'
#' @export
var_wtd_mean_cochran <- function(x,w){
  n = length(w)
  xWbar = wtd.mean(x,w, na.rm=TRUE)
  wbar = mean(w, na.rm=TRUE)
  out = n/((n-1)*sum(w)^2)*(sum((w*x-wbar*xWbar)^2) -
        2*xWbar*sum((w-wbar)*(w*x-wbar*xWbar))+xWbar^2*sum((w-wbar)^2))
  return(out)
#  if(out < 0){
#    var <- NA
#    if(!silent) warning('Negative value in var.wtd.mean.cochran:sqrt() NAs produced')
#  }else{
#    var <- sqrt(out*n)
#  }
#  return(var)
}

#' \code{origin_backtracking} for recursive backtracking origin estimation (Manitz et al., 2016)
#'
#' @param start_with_event_node logical specifying whether backtracking only starts from nodes that experienced events (for \code{type='backtracking'})
#' @return \code{origin_backtracking} returns an object of class \code{origin}, list with 
#' \itemize{
#'  \item \code{est} origin estimate
#'  \item \code{aux} \code{data.frame} with auxiliary variables 
#'    \itemize{
#'       \item \code{id} as node identifier, 
#'       \item \code{events} for event magnitude, and
#'       \item \code{bcount} for backtracking counts, how often backtracking identifies this source node.
#'    }
#'  \item \code{type = 'backtracking'} backtracking origin estimation
#'  }
#'
#' @examples
#' # backtracking origin estimation (Manitz et al., 2016)
#' ob <- origin(events=delayGoe[10,-c(1:2)], type='backtracking', graph=ptnGoe)
#' summary(ob)
#' plot(ob, start=1)
#' performance(ob, start=1, graph=ptnGoe)
#' 
#' @rdname origin
#' @import igraph
#' @export
origin_backtracking <- function(events, graph, start_with_event_node = TRUE, silent = TRUE){
    
   # input errors
    if(!is.vector(events)) events <- as.vector(events)
#    if(is.null(names(events))) warning('\nWarning: events and node names cannot be matched - we assume that events and graph nodes have the same order')

   # match events and graph names
    idm <- match(names(events), V(graph)$name) #<FIXME> correct order?

   # remove events that not have nodes in the network
    nas <- which(is.na(idm))
    if(length(nas) > 0){
       if(!silent) message(cat('\nNote: data events for nodes given that are not in the network are removed:\n',names(events)[nas]))
       events <- events[-nas]
       idm <- idm[-nas]
    }

   # add events to graph
    if(length(idm)){
       V(graph)$events <- events[idm]
    }else{ 
       V(graph)$events <- events 
       if(!silent) warning('\nWarning: events and node names cannot be matched - we assume that events and graph nodes have the same order')
    }

   # initialize result object
    result <- integer(length(events))

    K <- vcount(graph)

    for(k in 1:K){
        current_node <- k
        # NA handling
        if(is.na(V(graph)[current_node]$events)){ next }
        # check starting node   
        if(V(graph)[current_node]$events > 0 | !start_with_event_node){
            changed <- TRUE
            # as long there is a connected node with larger events
            while(changed){
                changed <- FALSE
                # check neighbors events
                nb <- neighbors(graph, current_node)
                current <- V(graph)[current_node]$events
                nb_delay <- unlist(V(graph)[nb]$events)
                # compare neighbours events with events at current node
                if(current < max(c(0, nb_delay), na.rm=TRUE)){
                    changed <- TRUE
                    current_node <- nb[which.max(nb_delay)]
                }else{
                    break # for secruity if all NA
                }
            }
            # update result object
            if( V(graph)[current_node]$events > 0){          
                result[current_node] <- result[current_node] + 1
            }
        }
    }

    # return argument
    aux <- data.frame(id=1:length(events), events = t(events), bcount=result)
    colnames(aux) <- c('id', 'events','bcount')
    ret <- list(est = which.max(result), aux = aux, type='backtracking')
    class(ret) <- 'origin'
    return(ret)
}

#origin_centrality <- function(x) UseMethod("origin_centrality")
#' \code{origin_centrality} for centrality-based origin estimation (Comin et al., 2011)
#'
#' @param graph igraph object specifying the underlying network graph (for \code{type='backtracking'} and \code{type='centrality'})
#' @param silent locigal, should the messages be suppressed?
#' @return \code{origin_centrality} returns an object of class \code{origin}, list with 
#' \itemize{
#'  \item \code{est} origin estimate
#'  \item \code{aux} \code{data.frame} with auxiliary variables 
#'    \itemize{
#'       \item \code{id} as node identifier, 
#'       \item \code{events} for event magnitude, and
#'       \item \code{cent} for node centrality (betweenness divided degree).
#'    }
#'  \item \code{type = 'centrality'} centrality-based origin estimation
#'  }
#'
#' @examples
#' # centrality-based origin estimation (Comin et al., 2011)
#' oc <- origin(events=delayGoe[10,-c(1:2)], type='centrality', graph=ptnGoe)
#' summary(oc)
#' plot(oc, start=1)
#' performance(oc, start=1, graph=ptnGoe)
#' 
#' @rdname origin
#' @import igraph
#' @export
origin_centrality <- function(events, graph, silent=TRUE){

    # input errors
    if(is.null(names(events))) stop('\nError: events and node names cannot be matched')
    if(!is.vector(events)) events <- as.vector(events)
    idm <- match(names(events), V(graph)$name)
    # remove events that not have nodes in the network
    nas <- which(is.na(idm))
    if(length(nas) > 0){
       if(!silent) message(cat('\nNote: data events for nodes given that are not in the network are removed:\n',names(events)[nas]))
       events <- events[-nas]
       idm <- idm[-nas]
    }

    # prepare auxiliary information of return argument
    aux <- data.frame(id=1:length(events), events=t(events), cent=NA)

    # add event info to network graph
    sub <- subset(events, select=events>0) #
    idv <- match(names(sub), V(graph)$name) #<FIXME>

    # no delays generated
    if(length(idv)==0){
       k0.hat <- NA
    }else{ # more than two nodes are infected
    # trivial solution
    if(length(idv)==1){
       k0.hat <- idv
    }else{ # for two-node networks betweenness cannot be computed
    if(length(idv)==2){
       k0.hat <- which.max(events)
    }else{ # more than two nodes are infected
#       if(length(sub) < 5){ # minimal spreading pattern
#         if(!silent) message(cat('\nNote: events for less than four network nodes, so that no origin detection can be performed\n'))
#	 aux <- data.frame(events=t(events), cent=NA, id=1:length(events))
#	 ret <- list(est = NA, aux = aux, type = 'centrality')
#         class(ret) <- 'origin'
#         return(ret)
#    }}
    
      # induced subgraph
      gsub <- induced.subgraph(graph, idv)

      # compute centrality
      cent <- betweenness(gsub)/degree(gsub)
      # estimate source
      k0.hat <- match(names(which.max(cent)), names(events))

      # return argument
#      aux <- data.frame(id=1:length(events), events=t(events), cent=NA)
      aux$cent[match(names(cent),names(events))] <- cent
    }}}
    # finalize return argument
    colnames(aux) <- c('id', 'events','cent')
    ret <- list(est = k0.hat, aux = aux, type='centrality')
    class(ret) <- 'origin'
    return(ret)
}


origin_multiple <- function(events, ...) UseMethod("origin_multiple")
#' Multiple origin estimation using community partitioning
#' 
#' @param events numeric vector of event counts at specific time point
#' @param type character specifying the method, \code{'edm'}, \code{'backtracking'} and \code{'centrality'} are available.
#' @param fast logical specifying community partitioning algorithm, default is \code{'TRUE'} that uses \code{\link{fastgreedy.community}}, \code{'FALSE'} refers to \code{\link{leading.eigenvector.community}} 
#' @param graph igraph object specifying the underlying network graph
#' @param no numeric specifying the number of supposed origins
#' @param distance numeric matrix specifying the distance matrix
#' @param ... parameters to be passed to origin methods \code{\link{origin_edm}}, \code{\link{origin_backtracking}} or \code{\link{origin_centrality}}
#' @return \code{origin_multiple} returns a list object with objects of class \code{\link{origin}} of length \code{no}
#'
#' @family origin-est
#'
#' @references Zang, W., Zhang, P., Zhou, C. and Guo, L. (2014) Discovering Multiple Diffusion Source Nodes in Social Networks. Procedia Computer Science, 29, 443-452. <DOI: 10.1016/j.procs.2014.05.040>
#'
#' @examples
#' data(ptnAth)
#' # backtracking
#' origin_multiple(events=delayAth[10,-c(1:2)], type='backtracking', graph=ptnAth, no=2)
#' # edm
#' athnet <- igraph::as_adjacency_matrix(ptnAth, sparse=FALSE)
#' p <- athnet/rowSums(athnet)
#' eff <- eff_dist(p)
#' origin_multiple(events=delayAth[10,-c(1:2)], type='edm', graph=ptnAth, no=2, distance=eff)
#' 
#' @import igraph
#' @export
origin_multiple <- function(events, type=c('edm', 'backtracking', 'centrality'), graph, no=2, distance, fast=TRUE, ...){
    # prepare graph
    K <- length(events)
    V(graph)$delay <- events
    V(graph)$station <- 1:K
   
    # define subgraph
    gsub <- induced.subgraph(graph, which(V(graph)$delay>0))

    # define communities 
    com <- if(fast) fastgreedy.community(gsub) else leading.eigenvector.community(gsub)
    comC <- cutat(com, no=no)

    # run source detection for each community
    res <- list()
    for(i in 1:no){
        if(sum(comC==i)<2){
	  res[i] <- unlist(V(gcom)$station)
	  next
	}
        # subset graph and distance
        gcom <- induced.subgraph(gsub, which(comC==i))
        # apply source detection
        if(type == 'edm'){
          distC <- distance[unlist(V(gcom)$station),unlist(V(gcom)$station)]
	  res[[i]] <- origin_edm(events = unlist(V(gcom)$delay), 
                                  distance = distC, ...)
#        res[i] <- unlist(V(gcom)$station)[fo.con(unlist(V(gcom)$delay), dist=distC)]
        }
        if(type == 'backtracking'){
	  res[[i]] <- origin_backtracking(events = unlist(V(gcom)$delay), 
                                 graph = gcom, ...) 
	}
        if(type == 'centrality'){
	  res[[i]] <- origin_centrality(events = unlist(V(gcom)$delay), 
                                 graph = gcom, ...) 
        }
    }
    return(res)
}


