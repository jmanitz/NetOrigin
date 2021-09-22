#' Computation of effective path distance 
#'
#' \code{eff_dist} computes the effective distance between all nodes in the network
#' 
#' @rdname eff_dist
#'
#' @param p numeric matrix, representing the transition probability matrix for the network graph
#' @return A numeric matrix, representing the effective distance between all nodes in the network graph.
#'
#' @references \itemize{
#'  \item Dijkstra, E. W. (1959): A note on two problems in connexion with graphs. Numerische Mathematik, 1, 269-271. <DOI: 10.1007/BF01386390>
#'  \item Brockmann, D. and Helbing, D. (2013): The hidden geometry of complex, network-driven contagion phenomena. Science, 342, 1337-1342. <DOI: 10.1126/science.1245200>
#'  \item Manitz, J. (2014): Statistical Inference for Propagation Processes on Complex Networks. Ph.D. thesis, Georg-August-University Goettingen. Verlag Dr. Hut, ISBN 978-3-8439-1668-4. Available online: \url{https://ediss.uni-goettingen.de/handle/11858/00-1735-0000-0022-5F38-B}.
#' }
#' 
#' @import igraph
#' @family distance
#'
#' @examples 
#' # compute effective shortest path distance
#' data(ptnAth)
#' require(igraph)
#' net <- igraph::as_adjacency_matrix(ptnAth, sparse=FALSE)
#' p <- net/rowSums(net)
#' eff <- eff_dist(p)
#'
#' @export
eff_dist <- function(p){
  K <- dim(p)[1]
  ppd <- NULL
  cat('Computing the effective distance between', K, 'nodes:\n 1')
  for(i in 1:K){
    cat(ifelse(i %% 100 == 0, paste('\n',i), '.'))
    ppd <- cbind(ppd,eff_dijkstra(p=p,start=i)[[1]])
  }
  cat('done\n')
  colnames(ppd) <- rownames(ppd) <- make.names(colnames(p))
  return(ppd)
}

#' \code{eff_dijkstra} computes the shortest effective paths using the dijkstra algorithm 
#'
#' @param start start of path
#'
#' @rdname eff_dist
#' @family distance
#' @export
eff_dijkstra <- function(p, start){
  ### initialize graph information
  N <- dim(p)[2]
  
  # Q..Set of unoptimzed nodes: All nodes in the graph are unoptimized 
  Q <- 1:N
  # is start a valid node
  if(!start %in% Q) stop('start is not a valid node')
  # Unknown distance function from source to all v
  distance <- rep(Inf, times=N)
  # pistance from source to source is zero
  distance[start] <- 0
  # Previous node in optimal path from source
  previous <- rep(NA, times=N)

  ### main loop
  while(length(Q)>0){ # as long there are unoptimized nodes
    # define shortest path for u: nodes in Q with smallest value in distance
    u <- Q[which.min(distance[Q])]
#    if(!is.finite(u)) break
    # delete u from Q
    Q <- Q[-match(u,Q)]
    # for all neighbours v of u 
    for(v in which(p[,u]>0)){
      # where v has not yet been removed from Q.
      if(v %in% Q){
        # update of distance
        alt <- distance[u] +(1 - log(p[v,u]) )
        if(alt < distance[v]){
          distance[v] <- alt
	  previous[v] <- u
        }
      }
		}
  }
  return(list(distance, previous))
}

# spd_dijkstra and compare with shortest paths from igraph
#' \code{spd_dijkstra} computes the shortest paths using the dijkstra algorithm 
#'
#' @examples
#' # compute shortest path distance
#' data(ptnAth)
#' athnet <- as_adj(ptnAth, sparse=FALSE)
#' spd <- spd_dijkstra(athnet, start=1)
#' 
#' # compare calculations with the one from igraph
#' spd_igraph <- igraph::distances(ptnAth, v=1, algorithm='dijkstra')
#' all(spd[[1]] == spd_igraph)
#'
#' @rdname eff_dist
#' @family distance
#' @export
spd_dijkstra <- function(p, start){
  ### initialize graph information
  N <- dim(p)[2]
  
  # Q..Set of unoptimzed nodes: All nodes in the graph are unoptimized 
  Q <- 1:N
  # is start a valid node
  if(!start %in% Q) stop('start is not a valid node')
  # Unknown distance function from source to all v
  distance <- rep(Inf, times=N)
  # pistance from source to source is zero
  distance[start] <- 0
  # Previous node in optimal path from source
  previous <- rep(NA, times=N)

  ### main loop
  while(length(Q)>0){ # as long there are unoptimized nodes
    # define shortest path for u: nodes in Q with smallest value in distance
    u <- Q[which.min(distance[Q])]
#    if(!is.finite(u)) break
    # delete u from Q
    Q <- Q[-match(u,Q)]
    # for all neighbours v of u 
    for(v in which(p[,u]>0)){
      # where v has not yet been removed from Q.
      if(v %in% Q){
        # update of distance
        alt <- distance[u] + p[v,u]
        if(alt < distance[v]){
          distance[v] <- alt
	  previous[v] <- u
        }
      }
		}
  }
  return(list(distance, previous))
}
