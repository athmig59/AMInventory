#' @name Normal_Functions
#' @rdname Normal_Functions
#'
#' @title Standardized normal distribution functions in Axsäter's notation
#'
#' @author Athanasios Migdalas, \email{athmig59@gmail.com}
#'
#' @param {x} {true function argument}
#' @param {Q} {parameter representing order quantity}
#' @param {R} {parameter representing reorder point}
#' @param {mu} {parameter representing demand mean}
#' @param {sigma} {parameter representing demand's standard deviation}
#' @param {L} {parameter representing lead-time}
#'
#' @return a numeric value or NULL in case of erroneous input
#'
#' @details
#' The following functions are made available:
#'
#' 1) Standardized Loss Function $G(x)=\\int_x^\\infty (v-x)\\phi(v)dv$.
#' Named G(x) after Axsäter's book, Chapter 5
#'
#' 2) The integral
#' $H(x) = \\int_x^\\infty G(v)dv$
#' of the standrardized loss function G(x)
#' as defined in Axsäter's book, Chapter 5.
#'
#' 3) Normal distribution functions in Axsäter's notation of Chapter 5:
#' - Standardized normal density function phi(x), i.e., φ(x).
#' - Standardized normal distribution function Phi(x), i.e.,  Φ(x).
#' - Standardized normal inverse distribution function PhInv(x), i.e., $Phi^\{-1\}()$.
#'
#' 4) Distribution and density functions of the inventory level in steady state (Axsäter's notation in Chapter 5)
#' - Distribution function F(x) with additional input parameters Q, R, μ, σ, L.
#' - Density function f(x) with additional input parameters Q, R, μ, σ, L.
#'
#' @references
#'  - Axsäter, S. (2015) Inventory Control, 3d Edition, Springer
#'
NULL
#'
#' @rdname Normal_Functions
#' @examples
#' G(0.45)
#' @export
G <- function(x){
  return( stats::dnorm(x) - x * (1-stats::pnorm(x)) )
}
#'
#' @rdname Normal_Functions
#' @examples
#' H(0.55)
#' @export
H <- function(x){
  return( 0.5 * ( (x**2+1) * (1-stats::pnorm(x)) - x * stats::dnorm(x) ) )
}
#'
#' @rdname Normal_Functions
#' @examples
#' phi(0.58)
#'@export
phi <- function(x){
  return( stats::dnorm(x) )
}
#'
#' @rdname Normal_Functions
#' @examples
#' Phi(0.55)
#'@export
Phi <- function(x){
  return( stats::pnorm(x) )
}
#'
#' @rdname Normal_Functions
#' @examples
#' PhInv(0.56)
#'@export
PhInv <- function(x){
  return( stats::qnorm(x) )
}
#'
#' @rdname Normal_Functions
#' @examples
#' F(x=0.15,Q=50,R=65,mu=68,sigma=25,L=2)
#'@export
F <- function(x,Q,R,mu,sigma,L){
  if(Q == 0) { return(NULL) }
  Lsigma <- sqrt(L)*sigma
  if(Lsigma == 0) { return(NULL) }
  Lmu <- L*mu
  z1 <- (R-x-Lmu)/Lsigma
  z2 <- (R+Q-x-Lmu)/Lsigma
  g1 <- G(z1)
  g2 <- G(z2)
  return( Lsigma * (g1 - g2)/Q )
}
#'
#' @rdname Normal_Functions
#' @examples
#' f(x=0.45, Q=50, R=65, mu=68, sigma=25, L=2)
#'@export
f <- function(x,Q,R,mu,sigma,L){
  if(Q == 0) { return(NULL) }
  Lsigma <- sqrt(L)*sigma
  if(Lsigma == 0) { return(NULL) }
  Lmu <- L*mu
  z1 <- (R-x-Lmu)/Lsigma
  z2 <- (R+Q-x-Lmu)/Lsigma
  g1 <- stats::pnorm(z1)
  g2 <- stats::pnorm(z2)
  return( (g2 - g1)/Q )
}
