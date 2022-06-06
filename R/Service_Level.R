#' @name Service_Level
#' @rdname Service_Level
#' 
#' @title Service level functions for single echelon stochastic demand inventory model
#' 
#' @description 
#' Functions to compute service level of type I (S1) and type II (S2) as 
#' well as shortage 
#' costs (b1 and b2), for given service level,
#'  in the case of normal distribution of the demand
#'  
#' @param {sigma} {standard deviation of demand}
#' @param {mu} {mean value of demand}
#' @param {R} {reorder point}
#' @param {Q} {order quantity}
#' @param {h:} {holding cost}
#' @param {b1} {shortage cost per unit and time unit}
#' @param {b2} {shortage cost per unit}
#' @param {S1} {service level of type I}
#' @param {S2} {service level of type II}
#'  
#' @return A single numeric value or NULL, if wrong input
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
#' @rdname Service_Level
#' @details
#' Function S1(): Service level S1
#' probability of no stockout per order cycle
#' @examples 
#' S1(sigma, mu, R)
#' @export
S1 <- function(sigma, mu, R){
  if(sigma <= 0)return(NULL)
  return( pnorm( (R-mu)/sigma ))
}
#'
#'
#' @rdname Service_Level
#' @details
#' Function S1b2(): Service level S1, i.e., probability of no stockout 
#' per order cycle in terms of shortage cost per unit, b2
#' (Approximation)
#' @examples 
#' S1b2(h,b2,mu,Q)
#' @export
S1b2 <- function(h, b2, mu, Q){
  x1 <- b2*mu
  if(x1 == 0)return(NULL)
  x2 <- 1 - Q*h/x1
  if(x2 <= 0)return(NULL)
  return(x2)
}
#'
#'
#' @rdname Service_Level
#' @details
#' Function S2(): S2 service level, i.e., "fill rate" fraction, also
#' S3 service level, i.e., "ready rate" fraction
#' @examples 
#' S2(sigma,mu,Q,R)
#' @export
S2 <- function(sigma, mu, Q, R){
  if( Q <= 0 ) return(NULL)
  if( sigma <= 0 ) return(NULL)
  x1 <- R - mu
  x2 <- x1 + Q
  return( 1 - (sigma/Q) * ( G(x1/sigma) - G(x2/sigma) ) )
}
#'
#'
#' @rdname Service_Level
#' @details
#' Function S2b1(): Fill-rate, ready-rate Service Level S2, S3
#' w.r.t. dhortage cost per unit and time unit, b1
#' @examples 
#' S2b1(h,b1)
#' @export
S2b1 <- function(h,b1){
  if(h+b1 <= 0)return(NULL)
  return(b1/(h+b1))
}
#'
#'
#' @rdname Service_Level
#' @details
#' Function b1S2(): Shortage cost per unit per time, b1,
#' in terms of fill-rate, ready-rate service level S2, S3
#' @examples 
#' b1S2(h,S2)
#' @export
b1S2 <- function(h, S2){
  if(S2 >= 1 ) return(NULL)
  return(h * S2 / (1-S2))
}
#'
#'
#' @rdname Service_Level
#' @details
#' Function b2S1(): Shortage cost per unit, b2, in terms of 
#' Service level S1, i.e., probability of no stockout 
#' per order cycle (Approximation)
#' @examples 
#' b2S1(h,mu,Q,S1)
#' @export
b2S1 <-function(h, mu, Q, S1){
  x1 <- 1-S1
  if(x1 <= 0)return(NULL)
  x2 <- mu * x1
  if (x2 <= 0)return(NULL)
  return(Q*h/x2)
}
