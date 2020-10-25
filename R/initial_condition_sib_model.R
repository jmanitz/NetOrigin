
#' Provide Initial Condition for Function SIB_SS
#'
#' \code{initial_condition_sib_model} Compute Initial Condition for Function SIB_SS
#'
#' @rdname initial_condition_sib_model
#'
#' @param POP_node vector, length represents number of cities/nodes; vector represents
#'        population at each node
#' @param sigma  symptomatic ratio, i.e., fraction of infected people that develop symptoms and are infective. 
#'               (The remaining fraction enters directly the recovered compartment.)
#' @param mu_B death rate of V.cholerae in the aquatic environment (day^-1)
#' @param theta contamination rate 
#' @param node_in index/indices for initial infected node(s)
#' @param in_prevalence initial prevalence of symptomatic infected in a node, default is 0.1\%
#'
#' @return a 5 x number of nodes matrix, each row represents the following for all the nodes:
#'         Row 1: number of suspectible people, i.e., population excpect infected and recovered for each node;
#'         Row 2: number of infected people;
#'         Row 3: number of recovered people;
#'         Row 4: bacteria concentration in equilibrium with infected individuals;
#'         Row 2: number of infected people, but representing cumulative cases
#'
#' @examples
#' data(envirparaList)
#' y0 <- initial_condition_sib_model(envirparaList$popu, envirparaList$sigma, 
#'   envirparaList$mu_B, envirparaList$theta, c(428))
#' @export


initial_condition_sib_model <- function(POP_node, sigma, mu_B, theta, node_in, in_prevalence=0.001) {
  # calculate nnodes
  nnodes = length(POP_node)
  #initial condition row: state variables, column: node
  #state variables 1:Susceptible, 2:Infected, 3:recovered, 4:Baceria conc., 5:cumulative cases
  initial_infected=round(POP_node[node_in]*in_prevalence)
  initial_recovered=round(((1-sigma)/sigma)*initial_infected)
  y0=matrix(0, 5, nnodes)
  y0[1,]=POP_node  #susceptible
  y0[1,node_in]=y0[1,node_in]-(initial_infected+initial_recovered)
  y0[2,node_in]=initial_infected
  y0[3,node_in]=initial_recovered
  y0[4,node_in]=theta[node_in]*in_prevalence/mu_B  #bacteria concentration in equilibrium with infected individuals
  y0[5,node_in]=initial_infected
  # return y0
  return(y0)
}





