#' @title Stochastic model for joint optimal determination of R and Q
#'
#' @author Athanasios Migdalas, \email{athmig59@gmail.com}
#'
#' @details
#' It implements a stochastic model for joint determination of R and Q
#' with uniform distribution of lead-time demand. The model is similar but not equivalent to that of
#' section 6.1.2 in Axsäter's book (see references below).
#'
#' The iterative method of alternatingly update order quantity and reorder point
#'  is applied in order to minimize the expected total cost.
#'
#' @references
#'  - Axsäter, S. (2015) Inventory Control, 3d Edition, Springer
#'
#' @seealso   JointRQb1, JointRQb1iter, JointRQss, JointRQhw, JointRQb1S2
#'
#' @param mu:      mean value of year demand
#' @param h:       holding cost
#' @param b1:      stockout cost
#' @param A:       fix order cost
#' @param a,b:    interval of uniform distribution for lead-time demand, default (0,100)
#' @param maxiter: maximum number of iterations, default 20
#' @param eps:  size of accepted relative difference between two consecutive values
#'
#' @return  Returns a list containing:
#' \itemize{
#'    \item{Q:}  {order quantity}
#'    \item{R:}   {reorder point}
#'    \item{ExpCost:}  {expected total cost}
#'    \item{Iterations:} {the number of iterations executed. If 0 iterations returned,  check solution for feasibility.}
#'  }
#'
#' @examples
#' \dontrun{
#' #Example #
#' #Suppose that the mean value of year demand is μ = 50, and that the lead-time demand
#' #is uniformly distrubted between a lower limit a=10 and an upper limit b=65.
#' #The costs are A = 100, h = 2, and b1 = 20.
#' #The order quantity Q and the reorder point
#' #R are sought.
#' #}
#'
#' JointRQu(mu=50, h=2, b1=20, A=100, a=10, b=65, maxiter=10)
#'
#' @export
JointRQu <- function(A, h, b1, mu, a=0, b=100, maxiter=20, eps=0.5e-5){

  G <- function(x,a,b){
    return( x**2 / (2*(b-a)) - ( b*x ) / (b-a) + b**2 / (2*(b-a)) )
  }

  Lmu <- (b+a)/2
  # Lsigma <- sqrt( (b-a)**2 / 12)

  Qhat <- sqrt( (2 * mu * ( A + b1 * Lmu)) / h)
  Qbar <- b1 * mu / h

  Q0 <- sqrt( 2 * A * mu / h )

  if(Qbar < Qhat) {
    n = 0
    R0 <- qunif( 1 - h * Q0 / (b1 * mu), min=a, max=b )
    EC <-  mu * A / Q0 + h* (0.5 * Q0 + R0 - Lmu) + b1 * mu * G( R0, a, b ) / Q0
    return( list( Q=Q1, R=R1, ExpCost = EC, Iterations=n) )
  }

  R0 <- 0

  n <- 1

  while(TRUE){
    # message("Q: ", Q0, " R: ", R0)
    R1 <- qunif( 1 - h * Q0 / (b1 * mu), min=a, max=b )
    Q1 <- sqrt( 2* mu * ( A + b1 * G( R1, a, b ) ) / h )
    if( abs(R1 - R0) / R1 <= eps ){ break }
    if( abs(Q1 - Q0) / Q1 <= eps ){ break }
    if( n >= maxiter ) { break }
    n <- n + 1
    R0 <- R1
    Q0 <- Q1
  }

  EC <-  mu * A / Q1 + h* (0.5 * Q1 + R1 - Lmu) + b1 * mu * G( R1, a, b ) / Q1

  return( list( Q=Q1, R=R1, ExpCost = EC, Iterations=n) )
}
