
#' Stochastic SIB model for infected cases simulation
#'
#' \code{SIB_SS} Stochastic SIB model for infected cases simulation
#'
#' @rdname SIB_SS
#'
#' @param mu population natality and mortality rate (day^-1)
#' @param beta contact rate 
#' @param rho immunity loss rate (day^-1)
#' @param sigma symptomatic ratio, i.e., fraction of infected people that develop symptoms and are infective. 
#'        (The remaining fraction enters directly the recovered compartment.) 
#' @param gamma rate at which people recover from cholera (day^-1)
#' @param alpha cholera induced mortality rate (day^-1)
#' @param mu_B death rate of V.cholerae in the aquatic environment (day^-1) 
#' @param m parameter for infection force, default value is 0.3
#' @param theta contamination rate 
#' @param nnodes number of nodes/cities
#' @param POP_nodes vector, length represents number of cities/nodes; vector represents
#'        population at each node
#' @param fluxes matrix, number of nodes x number of nodes 
#'        where each row contains the probabilities a person travels from the given city (by Row Index) to another city (by Column Index).
#' @param time_sim time steps for simulation, e.g., seq(0, 100, 0.1)
#' @param y0 initial condition for SIB_SS, output of 'initial_condition_SIB_SS'
#'
#' @return a matrix, nnodes x number of time steps, representing number of new cases at each node, each time step
#'
#' @examples
#' library(NetOrigin)
#' data(envirPara)
#' y0 <- initial_condition_SIB_SS(popu, sigma, mu_B, theta, c(428))
#' time_sim=seq(0, 100, by=0.1)
#' simu.list = SIB_SS(mu = mu, beta = beta, rho = rho, sigma = sigma, gamma = gamma,
#'                    alpha = alpha, mu_B = mu_B, theta = theta, nnodes = length(POP_node), POP_node = popu,
#'                    fluxes = humanmob.mass, time_sim = time_sim, y0 = y0)
#' @export


SIB_SS <- function(mu, beta, rho, sigma, gamma,
                   alpha, mu_B, m = 0.3, theta, nnodes, POP_node, fluxes, 
                   time_sim, y0) {
  # m is gravity parameter
  POP_node_SS=POP_node    #population for the stochastic simulator
  POP=sum(POP_node_SS);   #total population
  S=y0[1, ] 
  sumS=sum(S) # susceptibles
  I=y0[2, ]
  sumI=sum(I) # Infected
  R=y0[3, ] 
  sumR=sum(R) # recovered
  cumcase=y0[5, ]   # cumulative cases
  B=y0[4, ]         # Bacteria concentration
  ## initial values
  FORCE=((1-m)*beta*B/(1+B)+fluxes%*%(beta*B/(B+1))*m)  #force of infection for each node
  INFECTION=FORCE*S;                                    #infection rate for each node
  sumINF=sum(INFECTION) # sum of infection rates
  
  t=time_sim[1] # time 
  t_prev=t # previous time at which B was updated
  index_t=1
  
  I_node_t_SS=matrix(0, nnodes, length(time_sim)) # infected for each node and time point
  cumcases_node_t_SS=matrix(0, nnodes, length(time_sim)) # cumulative cases for each node and time point
  I_node_t_SS[, index_t]=I
  cumcases_node_t_SS[, index_t]=cumcase # initial value
  
  while (t<time_sim[length(time_sim)]) {
    ## below should be included in the while loop, as SIB_SS.m
    #  1 asym infection   2 sym infection    3 recovery  4immunity loss  5 birth  6 death     7 cholera_death
    event_rate=c(sumINF*(1-sigma), sumINF*(sigma), sumI*gamma, sumR*rho, 
                 POP*mu, POP*mu, sumI*alpha)  #total rate for the seven types of events
    deltat=-log(runif(1, 0, 1))/sum(event_rate);  #timestep to the next event
    event=sample(1:7, size=1, prob = event_rate)
    if (event==1) {
      # asym infection
      node=sample(1:length(INFECTION), size = 1, prob = INFECTION)
      # update variables
      S[node]=S[node]-1
      R[node]=R[node]+1
      sumS=sumS-1
      sumR=sumR+1
      INFECTION[node]=INFECTION[node]-FORCE[node]
      sumINF=sumINF-FORCE[node]
    } else if (event==2) {
      # sym infection
      node=sample(1:length(INFECTION), size = 1, prob = INFECTION)
      # update variables
      S[node]=S[node]-1
      I[node]=I[node]+1
      cumcase[node] = cumcase[node]+1
      sumS=sumS-1
      sumI=sumI+1
      INFECTION[node]=INFECTION[node]-FORCE[node]
      sumINF=sumINF-FORCE[node]
    } else if (event==3) {
      # recovery
      node=sample(1:length(INFECTION), size = 1, prob = I)
      I[node]=I[node]-1
      R[node]=R[node]+1
      sumI=sumI-1
      sumR=sumR+1
    } else if (event==4) {
      # immunity loss
      node=sample(1:length(INFECTION), size = 1, prob = R)
      R[node]=R[node]-1
      S[node]=S[node]+1
      sumR=sumR-1
      sumS=sumS+1  
      INFECTION[node]=INFECTION[node]+FORCE[node]
      sumINF=sumINF+FORCE[node]
    } else if (event==5) {
      # birth
      node=sample(1:length(INFECTION), size = 1, prob = POP_node_SS)
      S[node]=S[node]+1
      POP_node_SS[node]=POP_node_SS[node]+1
      sumS=sumS+1
      POP=POP+1
      INFECTION[node]=INFECTION[node]+FORCE[node]
      sumINF=sumINF+FORCE[node]
    } else if (event==6) {
      # death
      node=sample(1:length(INFECTION), size = 1, prob = POP_node_SS)
      comp=sample(1:3, size=1, prob=c(S[node], I[node], R[node]))
      if (comp == 1) {
        S[node]=S[node]-1
        POP_node_SS[node]=POP_node_SS[node]-1
        sumS=sumS-1
        POP=POP-1
        INFECTION[node]=INFECTION[node]-FORCE[node]
        sumINF=sumINF-FORCE[node]
      } else if (comp == 2) {
        I[node]=I[node]-1
        POP_node_SS[node]=POP_node_SS[node]-1
        sumI=sumI-1
        POP=POP-1
      } else {
        R[node]=R[node]-1
        POP_node_SS[node]=POP_node_SS[node]-1
        sumR=sumR-1
        POP=POP-1
      }
    } else {
      # cholera death
      node=sample(1:length(INFECTION), size = 1, prob = I)
      I[node]=I[node]-1
      sumI=sumI-1
      POP_node_SS[node]=POP_node_SS[node]-1
      POP=POP-1
    }
    
    t=t+deltat
    if (t >= time_sim[index_t+1]) {
      ff=(theta / POP_node_SS) * I/mu_B
      B=(B-ff) * exp(-(t-t_prev)*mu_B) + ff
      FORCE=((1-m)*beta*B/(1+B)+fluxes%*%(beta*B/(B+1))*m)
      INFECTION=FORCE*S;
      sumINF=sum(INFECTION)
      
      t_prev=t
      index_t=index_t+1
      cumcases_node_t_SS[,index_t]=cumcase
      I_node_t_SS[,index_t]=I
    }
  }
  
  cumcases_node_t_SS_base = matrix(0, dim(cumcases_node_t_SS)[1], dim(cumcases_node_t_SS)[2])
  cumcases_node_t_SS_base[, 2:dim(cumcases_node_t_SS_base)[2]] = cumcases_node_t_SS[, 1:(dim(cumcases_node_t_SS)[2]-1)]
  cases_node_t_SS = cumcases_node_t_SS - cumcases_node_t_SS_base
  ## return
  return(cases_node_t_SS)
}





