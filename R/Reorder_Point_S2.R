#' @title Reorder point computation for given service level type 1
#'
#' @author Athanasios Migdalas, \email{athmig59@gmail.com}
#'
#' @description
#' Stochastic, single echelon system is assumed with normally distributed demand.
#' An exact value and an approximative (heuristic) of the reorder point is returned 
#' in accordance with Chapter 5 in Axsäter's book.
#' See sections 5.6 and 5.7.2
#'
#' @references
#'  - Axsäter, S. (2015) Inventory Control, 3d Edition, Springer
#'
#' @seealso Plot_R_vs_S1, Reorder_Point_S1
#' 
#' @param {S1} {numeric, service level of type 1}
#' @param {sigma} {numeric, standard deviation of lead-time demand}
#' @param {mu} {numeric, mean of lead-time demand}
#' @param {Q} {numeri, order quantity} 
#' 
#' @return A list containing:
#'   \itemize{
#'     \item{\code{R}} {the "exact" value of the reorder point}
#'     \item{\code{Rheu}} {the heuristic value of the reorder point}
#'     \item{\code{Precision}} {the precision by which the "exact" value has been computed}
#'     \item{\code{Precision2}} {the precision by which the heuristic value has been computed
#'   }} 
#'   
#' @details 
#'   The exact value of R is computed by solving numerically the equation (5.52) in Axsäter's book.
#'   The function \code{uniroot} from package \code{stats} is used for this purpose and \code{Precision}
#'   refers to this numerical approach. Similarly, \code{uniroot} is used to solve equation (5.53) in
#'   Axsäter's book for the heuristic value of \code{R}. \code{Precision} and \code{Precision2} provideReorder_Point_S2(Q=5, S2=0.679, sigma=sqrt(10), mu=10)
#'   a means of judging the quality of the computed values of \code{R} and \code{Rheu}.
#'   The comparison of \code{R} and \code{Rheu} lets one determine how well 
#'   equation (5.53) approximates equation (5.52)
#'   
#' @examples 
#' ## ** Example 1 based on Example 5.4 in Axsäter **
#' ## For order quantity Q=5 and service level S2 = 0.679, lead-time demand with μ=10 and σ=sqrt(10),
#' ## the corresponding reorder point, R, is sought.
#' ##
#' Reorder_Point_S2(Q=5, S2=0.679, sigma=sqrt(10), mu=10)
#'
#' ## ** Example 2 **
#' ## For order quantity and lead-time demand as above but for sevice level S2 of 85%.
#' ## 
#' Reorder_Point_S2(Q=5, S2=0.85, sigma=sqrt(10), mu=10)
#'    
#' 
#' @export
Reorder_Point_S2 <- function(Q, S2, sigma, mu){
  
  G <- function(x){
    return( stats::dnorm(x) - x * (1-stats::pnorm(x)) )
  }
  
  S2fun <- function(R){
    x1 <- (R-mu)/sigma
    x2 <- (R+Q-mu)/sigma
    return( S2 - 1 + (sigma/Q) * ( G(x1) - G(x2) ) )
  }
  
  S2heu <- function(R){
    x1 <- (R-mu)/sigma
    return( S2 - 1 + (sigma/Q) * G(x1) )
  }
  
  # determina a bracketing interval
  
  a <- -mu; b <- mu; h <- Q
  
  while(TRUE){
    if( sign( S2fun(a) )  != sign( S2fun(b) ) ){
      rint <- c(a, b) 
      break
    }else{
      a <- a - h
      b <- b + h
    }
  }
  
  # curve( S2fun, from=a, to=b)
  # curve( S2heu, from=a, to=b)
  # abline(h=0)
  
  sol1 <- stats::uniroot( S2fun, rint, tol=0.5e-10 )
  
  sol2 <- stats::uniroot( S2heu, rint, tol=0.5e-10 )
  
  return( list(R=sol1$root, Rheu=sol2$root, 
               Precision=sol1$estim.prec,  Precision2=sol2$estim.prec) )
}
