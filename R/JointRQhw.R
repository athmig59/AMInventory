#' @title Stochastic model for joint optimal determination of R and Q
#'
#' @author Athanasios Migdalas, \email{athmig59@gmail.com}
#'
#' @details
#' It implements a stochastic model for joint determination of R and Q
#' with normal distribution of demand and stock out cost b2 per unit.
#'
#' Model Assumptions:
#'
#'1) Demand with yearly mean mu (μ) and standard deviation sigma (σ)
#' 2) Stockout cost b2 per unit
#' 3) Holding cost incurs at a rate h*IL per year, IL being the inventory level

#' The model is similar to that of section 6.1.2 in Axsäter's
#' book (see Chapter 4 in Hadley and Whitin and Chapter 4 in Snyder and Shen).
#'
#' The iterative method of alternatingly update order quantity and reorder point
#'  is applied in order to minimize the expected total cost.
#'
#' @references
#'  - Axsäter, S. (2015) Inventory Control, 3d Edition, Springer
#'  - Hadley, G. and Whitin, T. M. (1963) Analysis of Inventory Systems, Prentice-Hall
#'  - Snyder, L. V. and Shen, Z.-J. M. (2011) Fundamentals of Supply Chain Theory, Wiley
#'
#' @seealso JointRQb1, JointRQb1iter, JointRQss, JointRQb1S2, JointRQu
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
#' #Suppose that the mean value of year demand is normally distributed wit μ = 50 and σ=20.
#' #The lead-time is 4 weeks. The costs are A = 100, h = 2, and b1 = 20.
#' #The order quantity Q and the reorder point R are sought.
#' #}
#'
#' JointRQhw(mu=50, sigma=20, h=2, b2=20, A=100, L=4/52, maxiter=10)
#'
#' @export
JointRQhw <- function(A, h, b2, mu, sigma, L, maxiter=20, eps=0.5e-7){

  Lmu <- mu * L
  Lsigma <- sigma * sqrt(L)

  G <- function(x){
    return( dnorm(x) - x * (1-pnorm(x)) )
  }

  Q0 <- sqrt( (2*A*mu) / h)
  R0 <- 0

  n <- 1

  while(TRUE){
    R <- qnorm( 1 - h * Q0 / (b2 * mu), mean = Lmu, sd = Lsigma )

    #message("Q: ", Q0, " R: ", R)

    Q1 <- sqrt( 2 * mu  * ( A + b2 * Lsigma  *  G( ( R - Lmu ) / Lsigma ) )/ h  )

    n <- n + 1
    if( abs( Q0 - Q1 ) / Q1 <= eps ){ break }
    if( abs( R0 - R ) / R <= eps ){ break }
    if( n >= maxiter ){ break }
    Q0 <- Q1
    R0 <- R
  }

  EC <- h * ( R - Lmu + 0.5 * Q1 ) + A * mu / Q1 + b2 * mu * Lsigma * G( ( R - Lmu )/Lsigma )  / Q1

  return( list( Q=Q1, R=R, ExpCost=EC, Iterations=n))
}
