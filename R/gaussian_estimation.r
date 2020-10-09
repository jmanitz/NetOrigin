# GAUSSIAN SOURCE ESTIMATION JM 16/10/14, modified 17/01/02
#' \code{origin_gaussian} for Gaussian source estimation
#' 
#' @param mu expected propagation delay, default \code{mu=1}. If \code{length(mu)==1}, iid Gaussian distribution is assumed. 
#' @param sigma standard deviation of propagation delay, default \code{sigma=1}. If \code{length(sigma)==1}, iid Gaussian distribution is assumed. 
#' @param sigma_delay numeric, standard deviation of the measurement delay, default \code{sigma_delay=1}
#' @param prior numeric vector specifying the prior probability of each network node to be the source, default \code{NULL} refers to equal priors (Pinto et al., 2012)
#' @param post logical whether posterior probability should be computed?  Default is \code{post=TRUE}.
#'
#' @return \code{origin_gaussian} returns an object of class \code{origin}, list with 
#' \itemize{
#'  \item \code{est} origin estimate
#'  \item \code{prob} posterior detection probability
#'  \item \code{aux} \code{data.frame} with auxiliary variables 
#'    \itemize{
#'       \item \code{id} as node identifier, 
#'       \item \code{events} for event magnitude, 
#'       \item \code{delta} for discriminant functon value, 
#'       \item \code{prior} prior,  
#'       \item \code{post} posterior probability.
#'    }
#'  \item \code{type = 'gaussian'} Gaussian source estimation
#'  }
#' 
#' @rdname origin
#' @import igraph 
#' @importFrom mvtnorm dmvnorm
#' @export
origin_gaussian <- function(events, graph, mu=1, sigma=1, sigma_delay=1, prior=NULL, post=TRUE, silent=TRUE){

#-----------------------------------------------------------------------------------------
    # input errors
    if(!is.vector(events)) events <- as.vector(events)
    if(is.null(V(graph)$set)) stop("Estimates cannot be computed: Define observer nodes and partition network graph using 'partition_network()'.")

    # initialization
    N <- vcount(graph)
    obs <- which(V(graph)$obs)
    ref <- which(V(graph)$ref)
    K <- length(obs)
    
    # delay observation
    delay <- events - events[ref]
    del <- delay[setdiff(obs,ref)]

    if(any(is.na(del))) stop('Missing values in sensor nodes.') #<TODO> this has also information <TODO>

    # add prior if NULL
    if(is.null(prior)) prior <- rep(1/N, N)

#-----------------------------------------------------------------------------------------
    ### iid Gaussian 
    if(length(mu)==1 && length(sigma)==1){
        ### network topology estimates
        est <- compute_est(graph)
        mu <- est$mu
        lam1 <- solve(est$Lambda)

        ### compute discriminant function
        # discriminant function
        delta <- integer(N)
        for(i in 1:N){
            delta[i] <- del %*% lam1 %*% mu[i,] - 1/2 * t(mu[i,]) %*% lam1 %*% mu[i,] 
        }
        # add prior
        delta  <- delta + log(prior)
        
        ### compute posterior probabilities - computationally expensive
        if(post){    
            dm <- integer(N)
            for(i in 1:N) dm[i] <- dmvnorm(del, mean = mu[i,], sigma = est$Lambda)
            post <- prior*dm/sum(prior*dm, na.rm=TRUE)
        }else post <- NA

#-----------------------------------------------------------------------------------------        
    ### jointly Gaussian    
    }else{
        
        # reporting delay
        repd <- sigma_delay^2 * (matrix(nrow=K-1, ncol=K-1, data=1) + diag(nrow=K-1, ncol=K-1))

        delta <- integer(N)
        dm <- integer(N)

        for(s in 1:N){
            # routing matrix
            sp_obs <- shortest_paths(graph, from=s, to=setdiff(obs, ref), output='epath')$epath
            sp_ref <- shortest_paths(graph, from=s, to=ref, output='epath')$epath[[1]]

            Cs <- matrix(nrow=K-1, ecount(graph), data=0)
            for(k in 1:(K-1)){
                Cs[k,match(sp_obs[[k]],E(graph))] <- 1
                Cs[k,match(sp_ref,E(graph))] <- -1
            }

            # compute estimates
            mus <- Cs %*% mu # mean
            Lambdas <- Cs %*% diag(sigma^2) %*%t(Cs) + repd  # covariance

            ### discriminant function
            delta[s] <- (-1/2)*log(det(Lambdas)) * t(del - mus) %*% solve(Lambdas) %*% (del - mus) + log(prior[s])

            ### posterior
            if(post){ dm[s] <- dmvnorm(del, mean = mus, sigma = Lambdas) }
 
        }

        ### compute posterior probabilities - computationally expensive
        if(post){ post <- prior*dm/sum(prior*dm, na.rm=TRUE)
        }else post <- NA   

    }
#-----------------------------------------------------------------------------------------        

    ### derive estimator
    est <- which.max(delta)  
 
    ### return argument
    aux <- data.frame(id=1:N, events = events, delta=delta, prior=prior, post=post)
    ret <- list(est = est, prob = post[est], aux = aux, type='gaussian')
    class(ret) <- 'origin'
    return(ret)
}

prior <- function(x, ...) UseMethod("prior")
#' compute prior distribution for Gaussian source estimator 
#'
#' @param graph igraph object specifying the underlying network graph
#' @param type character, specifying the prior type to be computed (see details) 
#' @param zeta numeric, scaling parameter for informative priors
#' @param center integer, identifying the center of informative priors 
#'
#' @details Different prior specifications are implemented: 
#' \itemize{
#'    \item \code{type='equal'} assigns equal priors to all nodes, which equals in the Gaussian source estimator without priors (Pinto et al., 2012)
#'    \item \code{type='random'} is an non-informative prior with \deqn{\pi \sim N(\mu=1/N, \sigma=1/N^4)}
#'    \item \code{type='ring'} is an informative prior with exponential decay according to the shortest path distance from the \code{center} node and scaling parameter \code{zeta}, i.e. \deqn{\pi_i(s) = \exp(-\zeta \cdot d(s,i))}
#' } All priors are normalized, so that \deqn{\sum_{i=1}^N \pi_i = 1}
#' 
#' @references \itemize{
#'   \item Pinto, P. C., P. Thiran, and M. Vetterli (2012). Locating the source of diffusion in large-scale networks. Physical Review Letters 109 (6).
#' }
#' 
#' @return numeric vector of length \code{vcount(graph)}.
#'
#' @examples
#' # define network
#' set.seed(160716)
#' g <- igraph::sample_pa(30, directed=FALSE)
#' # compute prior
#' prior(g, type='equal')
#' prior(g, type='random')
#' prior(g, type='ring', zeta=1, center=8)
#' 
#' @importFrom stats rnorm
#' @import igraph
#' @seealso \code{\link{origin_gaussian}}
#' @export
prior <- function(graph, type=c('equal', 'random', 'ring'), zeta=1, center=NULL){
    N <- vcount(graph)
    type <- match.arg(type)                  
    if(type == 'equal'){
      pi <- rep(1/N,N) 
    }
    if(type == 'random'){ 
      pi <- rnorm(N, mean=1/N, 1/N^2)     
    }
    if(type == 'ring'){
 	x <- distances(graph, v=center)
	pi <- as.vector(exp(-zeta*x))
    }

    # normalize
    pi_norm <- pi/sum(pi)
    return(pi_norm)
}


### estimate network topology
#' @title compute network topology estimates for Gaussian source estimate
#' @description compute network topology estimates of deterministic delay and delay covariance for Gaussian source estimator
#'
#' @param graph igraph object specifying the underlying network graph
#' @param mean integer specifying normalizing constant for deterministic delay, default \code{mean=1}
#' @param sigma integer specifying normalizing constant for delay covariance, default \code{sigma=1}
#'
#' @return \code{compute_est()} returns a list of 
#' \itemize{
#'  \item \code{mu} matrix, estimate for deterministic delay, dimension (K-1)xN
#'  \item \code{Lambda} matrix, estimate for delay covariance, dimension (K-1)x(K-1),
#' } where K is the number of sensors and N the number of network nodes.
#' 
#' @examples
#' # define network
#' set.seed(160716)
#' g <- igraph::sample_pa(30, directed=FALSE)
#' # partition network
#' g <- partition_network(graph = g, obs = c(1,12, 27))
#' # compute estimates
#' compute_est(graph = g)
#' 
#' @import igraph
#' @export
compute_est <- function(graph, mean=1, sigma=1){ 

    # error handling 
    if(is.null(V(graph)$set)) stop("Estimates cannot be computed: Define observer nodes and partition network graph using 'partition_network()'.")
    
    obs <- which(V(graph)$obs)
    ref <- which(V(graph)$ref)
    K <- length(obs)

    # shortest path distances
    spd_obs <- distances(graph, to=obs) 
    spd_ref <- distances(graph, to=ref) 
#    # disconnected nodes
#    spd_obs[is.infinite(spd_ref),] <- spd_ref[is.infinite(spd_ref),] <- NA 
    # compute mu: 
    mu <- matrix(nrow=vcount(graph), ncol=K, data=NA)
    for(i in 1:K){ 
        mu[,i] <-  spd_obs[,i] - spd_ref
    }
    # compute Lamda
    sp_obs <- shortest_paths(graph, from=ref, to=obs, output='epath')$epath
    Lambda <- matrix(nrow=K, ncol=K, data=NA)
    #diag(Lambda) <- spd_obs[1,obs]
    for(i in 1:K){
        for(j in i:K){
            Lambda[i,j] <- Lambda[j,i] <- length(intersect(sp_obs[[i]], 
                                                           sp_obs[[j]]))
        }
    }
    # return object
    est <- list()
    rem <- which(obs==ref)
    est$mu <- mean*mu[,-rem]
    est$Lambda <- sigma*Lambda[-rem,-rem]
    return(est)
}

### as function: partition_network()
#' @title partitioning network for Gaussian source estimation
#' @description partitioning network graph into sets of observers, resolvable sets, and unresolvable sets L and U. 
#' @param graph igraph object, network graph
#' @param obs integer, observer nodes, while first node will be the reference node 
#' @return igraph object, with logical graph attributes \code{'obs'}, \code{'ref'}, \code{'Rset'}, \code{'Uset'}, and \code{'Lset'}. The node attribute \code{'set'} summarizes to partioning set assignment: 
#' \itemize{
#' \item \code{'O'} observer node
#' \item \code{'R'} resolvable node
#' \item \code{'L'} unresolvable node, adjacent to one observer
#' \item \code{'U'} unresolvable node, adjacent to a resolvable node from set R
#' }
#' 
#' @examples
#' # define network
#' set.seed(160716)
#' g <- igraph::sample_pa(30, directed=FALSE)
#' #partition network
#' g <- partition_network(graph = g, obs = c(1,12, 27))
#' plot_partition(g, vertex.size=15, edge.width=2)
#' 
#' @import igraph
#' @seealso plot_partition
#' @export
partition_network <- function(graph, obs){

    # initialize set attribute
    V(graph)$set <- NA
    V(graph)$obs <- V(graph)$ref <- V(graph)$Rset <- V(graph)$Uset <- V(graph)$Lset <- FALSE

    # assign observers
    V(graph)$set[obs] <- 'Observer'
    V(graph)$obs[obs] <- TRUE
    V(graph)$ref[obs[1]] <- TRUE 

    # assign resolvable nodes = union of paths between every pair of observers
    obs_paths <- shortest_paths(graph, from=obs, to=obs)
    paths_nodes <- unique(unlist(obs_paths))
    res <- setdiff(paths_nodes, obs)
    res <- V(graph)[res]
    V(graph)$set[res] <- 'Resolvable'
    V(graph)$Rset[res] <- TRUE 

    # assign L sets = adjacent to one observer
    x <- union(obs,res) # set of assigned nodes
    nodes <- obs        # neighbor search set
    while(length(nodes)>0){
        can <- unique(unlist(adjacent_vertices(graph, v=nodes)))
        nodes <- setdiff(can,x)
        V(graph)$set[nodes] <- 'Unresolvable L'
        V(graph)$Lset[nodes] <- TRUE 
        x <- union(x,V(graph)[nodes])
    }

    # assign U sets = unresolved adjacent to resolvable node
    u_nodes <- which(is.na(V(graph)$set))
    V(graph)$set[u_nodes] <- 'Unresolvable U'
    V(graph)$Uset[u_nodes] <- TRUE

    # return network 
    return(graph)
}

#' @title plot partitioned network
#' @description plot partitioned network with coloring according to sets of observers, resolvable sets, and unresolvable set L and U. 
#' @param graph igraph object specifying the network graph with attribute \code{'set'} defining the partion
#' @param ... further parameters to be passed to \code{\link{plot.igraph}}
#' 
#' @examples
#' # define network
#' set.seed(160716)
#' g <- igraph::sample_pa(30, directed=FALSE)
#' #partition network
#' g <- partition_network(graph = g, obs = c(1,12, 27))
#' plot_partition(g, vertex.size=15, edge.width=2)
#' 
#' @import igraph
#' @seealso partition_network plot_ptn
#' @export
plot_partition <- function(graph, ...){

    if(is.null(V(graph)$set)) stop("Partitioned network cannot be plotted: Define observer nodes and partition network graph using 'partition_network()'.")
    obs <- which(V(graph)$obs)
    
    # set alayout
    graph <- set_graph_attr(graph, "layout", layout_with_kk(graph))
    V(graph)$size <- 9; E(graph)$width <- 2

    # mark node set
    V(graph)$color <- as.numeric(factor(V(graph)$set))
    # highlight reference observer
    V(graph)$frame.color <- NA
    V(graph)$frame.color[V(graph)$ref] <- 'red'
    # highlight true source
    if(!is.null(V(graph)$origin)) V(graph)$color[V(graph)$origin] <- 'red'
    # plot network graph
    plot(graph, ...)
    if(!is.null(V(graph)$origin)){
       legend('bottomright', pch=20, cex=1.1, bty='n',
              col=c(categorical_pal(4),'red'), #[c(2:3,1,4)]
              legend=c(levels(factor(V(graph)$set)),'origin s*'))
    }else{
       legend(-1,0.4, pch=20, cex=1, bty='n',
              col=categorical_pal(4), legend=levels(factor(V(graph)$set)))
    }
    return(invisible(graph))
}

#' procedure for probability estimates for multi-class comparisons by pairwise coupling (Wu et al, 2004) 
#' @details The implementation is an adapted version of \code{kernlab::minpair} (see references).
#'
#' @param p2 symmetric matrix of pairwise probabilities
#' @return numeric vector with estimated probability P(y=i|x) for each class 
#' 
#' @references \itemize{
#'   \item Wu, Ting-Fan, Chih-Jen Lin, and Ruby C. Weng (2004). Probability estimates for multi-class classification by pairwise coupling. Journal of Machine Learning Research. 975-1005.
#'   \item Alexandros Karatzoglou, Alex Smola, Kurt Hornik, Achim Zeileis (2004). kernlab - An S4 Package for Kernel Methods in R. Journal of Statistical Software 11(9), 1-20. URL http://www.jstatsoft.org/v11/i09/
#' }
#' @export
couple_pair <- function(p2){

    # input handling
    nclass <- ncol(p2)
#    if(!isSymmetric(p2)) warning("Matrix is not symmetric.")

    ### define Q matrix (see formula Eq.20, p.981)
    # off-diagonal elements
    Q <- -p2*t(p2)
    # diagonal elements 
    diag(p2) <- 0
    diag(Q) <- colSums(p2^2)

    ### set up linear equation system (see formula Eq.21, p.981)
    # left side
    e <- rep(1, nclass)
    SQ <- cbind(rbind(Q, e), c(e, 0))
    # right side
    rhs <- rep(0, nclass + 1)
    rhs[nclass + 1] <- 1

    ### solve linear equation system  
    p <- solve(SQ, rhs)
    # normalize
    p <- p[-(nclass + 1)]/sum(p[-(nclass + 1)])

    return(p)
}



