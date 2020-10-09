TimeMin <- function(num.cases, thres = NA){
  # Compute the time when a node is infected
  ifelse(any(cumsum(num.cases) > thres),
         min(which(cumsum(num.cases) > thres)), 
         NA)
}


#' Inference Source via Gaussian source estimation with prior information
#'
#' \code{infer_source_bayesian} Compute posterior probabilities for all nodes and rank them
#'
#' @rdname infer_source_bayesian
#'
#' @param cases.node.day matrix, number of nodes x time points; entries represent number of cases
#' @param thres.vec vector, length represents number of cities/nodes, representing
#'        thresholds for cities/nodes that they are infected
#' @param obs.vec list of cities ids used as observers
#' @param mu.mat matrix- number of cities/nodes x number of observers, each row represents- 
#'        if this node is the source, the mean of arrival time vector
#' @param lambda.list a length-number of cities/nodes list, each element is a number of observers x number of observers matrix-
#'        if a node is the source, the covariance matrix for arrival time vector
#' @param poss.candidate.vec a boolean vector indicating if a node has the potential to be the source
#' @param prior vector, length - number of cities/nodes, prior for cities
#' @param use.prior boolean, TRUE or FALSE, if use prior, default TRUE
#'
#' @return a dataframe with columns 'nodes' and 'probab', indicating nodes indices and their posteriors
#'
#' @examples
#' library(NetOrigin)
#' # load training data
#' load('.../train_data.rda')
#' obs.vec <- (1:9)
#' candidate.thres <- 0.3
#' mu.lambda.list <- compute_mu_lambda(train.data, obs.vec, candidate.thres)
#' # load matrix representing number of cases per node per day, size-number of nodes x number of day
#' load('.../cases_node_day.rda')
#' # number of nodes
#' nnodes <- dim(cases.node.day)[1] 
#' # fixed threshold for all nodes - 10 infected people
#' thres.vec <- rep(10, nnodes)
#' # flat/non-informative prior
#' prior <- rep(1, nnodes) 
#' result2.df <- infer.source.bayesian(cases.node.day, 
#'                                     thres.vec,
#'                                     obs.vec,
#'                                     mu.lambda.list$mu.mat, mu.lambda.list$lambda.list, 
#'                                     mu.lambda.list$poss.candidate.vec,
#'                                     prior, TRUE)
#' @export


infer_source_bayesian <- function(cases.node.day, 
                                  thres.vec,
                                  obs.vec,
                                  mu.mat, lambda.list, 
                                  poss.candidate.vec,
                                  prior, use.prior = TRUE){
  ## constant parameters
  nnodes <- dim(cases.node.day)[1]
  
  ## extract 'observed.time.all'
  observed.time.all <- rep(NA, 851)
  for (node.num in 1:851) {
    if (length(which(cases.node.day[node.num, ] > thres.vec[node.num])) > 0) {
      observed.time.all[node.num] <- which(cases.node.day[node.num, ] > thres.vec[node.num])[1]
    } else {
      observed.time.all[node.num] <- NA
    }
  }
  ## extract info from 'poss.candidate.vec'
  num.source.candidate = sum(poss.candidate.vec)
  nodes.canbe.detected <- (1:851)[poss.candidate.vec]
  ## prepare 'prob.vec'
  logprob.vec <- rep(NA, num.source.candidate)
  for (j in 1:num.source.candidate) {
    if (is.na(lambda.list[[nodes.canbe.detected[j]]][1])) {
      logprob.vec[j] <- (-Inf)
    } else {
      logprob.vec[j] <-
        dmvnorm(observed.time.all[obs.vec],
                mean = mu.mat[
                  nodes.canbe.detected[j], ],
                sigma = lambda.list[[
                  nodes.canbe.detected[j]]],
                log = TRUE)
    }
  }
  if (!use.prior) {
    prob.vec <- exp(logprob.vec)
  } else {
    prob.vec <- exp(logprob.vec) * prior[poss.candidate.vec]
  }
  ## prepare final 'df'
  probnorm.vec <- prob.vec / sum(prob.vec)
  df <- data.frame(probab = probnorm.vec,
                   nodes = nodes.canbe.detected)
  df <- df[rev(order(df$probab)), ]
  
  ## return
  return(df)
}





