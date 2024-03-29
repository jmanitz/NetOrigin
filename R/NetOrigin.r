#' Origin Estimation for Propagation Processes on Complex Networks
#'
#' Performs different approaches for network-based source estimation: effective distance median, recursive backtracking, and centrality-based source estimation. Additionally, we provide public transportation network data as well as methods for data preparation, source estimation performance analysis and visualization.
#' 
#' The main function for origin estimation of propagation processes on complex network is \code{\link{origin}}. Different methods are available: effective distance median (\code{'edm'}), recursive backtracking (\code{'backtracking'}), and centrality-based source estimation (\code{'centrality'}).
#' For more details on the methodological background, we refer to the corresponding publications.
#'
#' @references \itemize{
#'   \item Manitz, J., J. Harbering, M. Schmidt, T. Kneib, and A. Schoebel (2017): Source Estimation for Propagation Processes on Complex Networks with an Application to Delays in Public Transportation Systems. Journal of Royal Statistical Society C (Applied Statistics), 66: 521-536.
#'   \item Manitz, J., Kneib, T., Schlather, M., Helbing, D. and Brockmann, D. (2014) Origin detection during food-borne disease outbreaks - a case study of the 2011 EHEC/HUS outbreak in Germany. PLoS Currents Outbreaks, 1. <DOI: 10.1371/currents.outbreaks.f3fdeb08c5b9de7c09ed9cbcef5f01f2>
#'   \item Comin, C. H. and da Fontoura Costa, L. (2011) Identifying the starting point of a spreading process in complex networks. Physical Review E, 84. <DOI: 10.1103/PhysRevE.84.056105>
#' }
#'
#' @author Juliane Manitz with contributions by Jonas Harbering
#' 
#' @docType package
#' @name NetOrigin
#' @aliases NetOrigin-package
NULL


