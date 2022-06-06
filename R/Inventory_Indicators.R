#' @name Inventory_Indicators
#' @rdname Inventory_Indicators
#'
#' @title Functions to evaluate inventory indicators such as level and position
#'
#' @description
#' Functions to compute inventory indicators such as level and position
#'
#' @param {mu} {mean value of demand}
#' @param {R} {reorder point}
#' @param {Q} {order quantity}
#' @param {L} {lead time}
#'
#' @return A single numeric value. NULL if wrong input.
#'
#' @references
#'  - Axsäter, S. (2015) Inventory Control, 3d Edition, Springer
#'
#'
#'@author Athanasios Migdalas, \email{athmig59@gmail.com}
#'
#'
NULL
#'
#' @rdname Inventory_Indicators
#' @details
#' Function IL(): Computets the average inventory level
#' @examples
#' IL(mu=100, Q=85, R=90)
#' @export
IL <- function(mu, Q, R){
  return( R + 0.5 * Q - mu )
}
#'
#' @rdname Inventory_Indicators
#' @details
#' Function IP(): Computes the average  inventory position
#' @examples
#' IP(Q=100, R=90)
#' @export
IP <- function(Q,R){
  return( R + 0.5 * Q)
}
#'
#' @rdname Inventory_Indicators
#' @details
#' Function EILm(): Computes the expected amount of backorders
#' @examples
#' EILm(mu=150,sigma=40,L=4,Q=100, R=95)
#' @export
EILm <- function(mu,sigma,L,Q,R){
  if(Q == 0) { return(NULL) }
  Lsigma <- sigma * sqrt(L)
  if(Lsigma == 0) { return(NULL) }
  Lmu <- L * mu
  H <- function(x){
    return( 0.5 * ( (x**2+1) * (1-stats::pnorm(x)) - x * stats::dnorm(x) ) )
  }
  z1 <- (R - Lmu)/Lsigma
  z2 <- (R + Q - Lmu)/Lsigma
  return( Lsigma**2 * ( H(z1) - H(z2) )/ Q )
}
