#' @title Stochastic model for joint optimal determination of R and Q
#'
#' @author Athanasios Migdalas, \email{athmig59@gmail.com}
#'
#' @details
#' It implements a stochastic model for joint determination of R and Q
#' with normal distribution of demand and stock out cost b1 per unit.
#'
#' Model Assumptions:
#'
#'1) Demand with yearly mean mu (μ) and standard deviation sigma (σ)
#' 2) Stockout cost b1 per unit
#' 3) Lead-time demand normally distributed with mean mu (μ') and standard deviation sigma (σ')
#'
#' The model  (see section 4.3.2 in Snyder and Shen's book) is similar to that of
#' section 6.1.2 in Axsäter's book.
#'
#' The iterative method of alternatingly update order quantity and reorder point
#'  is applied in order to minimize the expected total cost.
#'
#' @references
#'  - Axsäter, S. (2015) Inventory Control, 3d Edition, Springer
#'  - Snyder, L. V. and Shen, Z.-J. M. (2011) Fundamentals of Supply Chain Theory, Wiley
#'
#' @seealso JointRQb1, JointRQb1iter, JointRQu, JointRQhw, JointRQb1S2
#'
#' @param mu:      mean value of demand
#' @param sigma:   standard deviation of demand
#' @param L:       lead-time
#' @param h:       holding cost
#' @param b1:      stockout cost
#' @param A:       fix order cost
#' @param maxiter: maximum number of iterations, default 20
#' @param eps:  size of accepted relative difference between two consecutive values
#'
#' @return  Returns a list containing:
#' \itemize{
#'    \item{Q:}  {order quantity}
#'    \item{R:}   {reorder point}
#'    \item{ExpCost:}  {expected total cost}
#'    \item{Iterations:} {the number of iterations executed. If 0 iterations returned,  check solution for feasibility.}
#' }
#'
#' @examples
#' \dontrun{
#' #Example #
#' #Suppose that the mean value of demand is normally distributed with μ = 50 and σ=20.
#' #The lead-time is 4 weeks. The costs are A = 100, h = 2, and b1 = 20.
#' #The order quantity Q and the reorder point R are sought.
#' #}
#'
#' JointRQss(mu=50, sigma=20, h=2, b1=20, A=100, L=4, maxiter=10)
#'
#' @export
JointRQss <- function(A, h, b1, mu, sigma, L, maxiter=20, eps=0.5e-7){

  G <- function(x){
    return( dnorm(x) - x * (1-pnorm(x)) )
  }

  Lmu <- mu * L
  Lsigma <- sigma * sqrt(L)

  Qhat <- sqrt( (2 * mu * ( A + b1 * Lmu)) / h)
  Qbar <- b1 * mu / h

  Q0 <- sqrt( 2 * A * mu / h )

  #message("Qhat: ", Qhat, " Qbar: ", Qbar)

  if(Qbar < Qhat) {
    n = 0
    R0 <- qnorm( 1 - h * Q0 / (b1 * mu), sd = Lsigma, mean = Lmu )
    EC <-  mu * A / Q0 + h* (0.5 * Q0 + R0 - Lmu) + b1 * mu * Lsigma * G( (R0 - Lmu ) / Lsigma ) / Q0
    return( list( Q=Q0, R=R0, ExpCost = EC, Iterations=n))
  }

  R0 <- 0

  n <- 1

  while(TRUE){
    # message("Q: ", Q0, " R: ", R0)
    R1 <- qnorm( 1 - h * Q0 / (b1 * mu), sd = Lsigma, mean = Lmu )
    Q1 <- sqrt( 2* mu * ( A + b1 * Lsigma * G( ( R1-Lmu ) / Lsigma ) ) / h )

    if( abs(R1 - R0) / R1 <= eps ){ break }
    if( abs(Q1 - Q0) / Q1 <= eps ){ break }
    if( n >= maxiter ) { break }
    n <- n + 1
    R0 <- R1
    Q0 <- Q1
  }

  EC <-  mu * A / Q1 + h* (0.5 * Q1 + R1 - Lmu) + b1 * mu * Lsigma * G( (R1 - Lmu ) / Lsigma ) / Q1

  return( list( Q=Q1, R=R1, ExpCost = EC, Iterations=n) )
}
