#' @title Stochastic model for joint optimal determination of R and Q
#'
#' @author Athanasios Migdalas, \email{athmig59@gmail.com}
#'
#' @details
#' It implements a stochastic model for joint determination of R and Q
#' with normal distribution of demand and stock out cost b1 per unit
#' and time unit. The model is that of section 6.1.2 in Axsäter's book (see
#' references below)
#'
#' The Nelder-Mead method from the package "optimx" is applied for
#' direct minimization of the expected total cost.
#'
#' @references
#'  - Axsäter, S. (2015) Inventory Control, 3d Edition, Springer
#'
#' @seealso JointRQu, JointRQb1iter, JointRQss, JointRQhw, JointRQb1S2
#'
#' @param mu:      mean value of demand
#' @param sigma:   standard deviation of demand
#' @param L:       lead-time
#' @param h:       holding cost
#' @param b1:      stockout cost
#' @param A:       fix order cost
#' @param Qo:      initial guess of order quantity
#' @param Ro:      initial guess of reorder point
#'
#' @return A list containing:
#' \itemize{
#'    \item{Startpoint:}  {a vector with elements Qo, Ro as above}
#'    \item{Startcost:}   {the cost Co for these guesses}
#'    \item{Finalpoint:}  {a vector with elements Q*, R*}
#'    \item{Finalcost:}   {the cost C* for these final values}
#'    \item{Optimum:}     {TRUE if optimality of final point is verified, FALSE otherwise
#'                 (in case of FALSE, the final point can be used as a new start point)}
#' }
#'
#' @examples
#' \dontrun{
#' #Example 6.1 in Axsäter's book, page 111#
#' #Suppose that the  demand per time unit is normally
#' #distributed with μ = 50 and σ = 20. The costs are A = 100, h = 2, and b1 = 20.
#' #Suppose that the lead-time is L = 4. The order quantity Q and the reorder point
#' #R are sought.
#' #
#' # We need to provide initial guesses for R and Q. Thus, let Ro=224.76 and Qo= 70.71.
#' }
#'
#' JointRQb1(mu=50, sigma=20, L=4, h=2, b1=20, A=100, Qo=70.71, Ro=224.76)
#'
#' @export
JointRQb1 <- function(mu, sigma, L, h, b1, A, Qo, Ro){

  Lmu <- L * mu
  Lsigma <- sigma * sqrt(L)

  CRQ <- function(RQ){

    H <- function(x){
      return( 0.5 * ( (x**2+1) * (1-pnorm(x)) - x * dnorm(x) ) )
    }

    R <- RQ[1]
    Q <- RQ[2]

    z1 <- ( R - Lmu ) / Lsigma
    z2 <- ( R + Q - Lmu ) / Lsigma

    H1 <- H(z1)
    H2 <- H(z2)

    val <-  h * ( R + 0.5 * Q - Lmu ) + ( h + b1 ) * Lsigma**2 * ( H1 - H2 ) / Q + A * mu / Q

    return( val )
  }

  startpoint <- c(Ro=Ro,Qo=Qo)
  startcost  <- CRQ(startpoint)
  names(startcost) <- "Co"

  v <- optimx::optimx( par=startpoint, fn=CRQ, method="Nelder-Mead", control=list(trace=0))

  finalpoint <- c(v$R, v$Q)
  names(finalpoint) <- c("R*", "Q*")
  finalcost  <- v$value
  names(finalcost) <- "C*"

  opt <-  v$kkt1 && v$kkt2

  res <- list( Startpoint = startpoint, Startcost = startcost,
               Finalpoint = finalpoint, Finalcost = finalcost,
               Optimum = opt)

  return( res )
}
