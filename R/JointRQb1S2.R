#' @title Stochastic model for joint optimal determination of R and Q for given service level S2
#'
#' @author Athanasios Migdalas, \email{athmig59@gmail.com}
#'
#' @details
#' It implements a stochastic model for joint determination of R and Q
#' with normal distribution of demand and stock out cost b1 per unit
#' and time unit under constraint for service level of type II (S2).
#' The model is that of section 6.1.3 in Axsäter's book (see
#' references below)
#'
#' A Nonlinear Augmented Lagrangian SQP Based Approach, "solnp()" from the
#' package "Rsolnp" is applied for direct minimization of the expected total cost.
#'
#' @references
#'  - Axsäter, S. (2015) Inventory Control, 3d Edition, Springer
#'
#' @seealso  JointRQb1, JointRQb1iter, JointRQss, JointRQhw, JointRQu
#'
#' @param mu:      mean value of demand
#' @param sigma:   standard deviation of demand
#' @param L:       lead-time
#' @param h:       holding cost
#' @param b1:      stockout cost
#' @param A:       fix order cost
#' @param S2:      service level of type II
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
#' #Example 6.2 in Axsäter's book, page 113#
#' #Suppose that the  demand per time unit is normally
#' #distributed with μ = 50 and σ = 20. The costs are A = 100, h = 2, and b1 = 0.
#' #Suppose that the lead-time is L = 4 and the service level of type II requested at 90%.
#' #The order quantity Q and the reorder point R are sought.
#' #
#' # We need to provide initial guesses for R and Q. Thus, let Ro=224.76 and Qo= 70.71.
#' #}
#'
#' JointRQb1S2(mu=50, sigma=20, L=4, h=2, b1=0, A=100, S2=0.9, Qo=70.71, Ro=224.76)
#'
#' @export
JointRQb1S2 <- function(mu, sigma, L, h, b1, A, S2, Qo, Ro){

  Lmu <- L * mu
  Lsigma <- sigma * sqrt(L)

  startpoint <- c(Ro=Ro,Qo=Qo)


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

  S2con <- function(RQ){

    G <- function(x){
      return( dnorm(x) - x * (1-pnorm(x)) )
    }

    R <- RQ[1]
    Q <- RQ[2]

    z1 <- ( R - Lmu ) / Lsigma
    z2 <- ( R + Q - Lmu ) / Lsigma

    G1 <- G(z1)
    G2 <- G(z2)

    val <- 1 - Lsigma * ( G1 - G2) / Q - S2

    return(val)
  }

  startcost  <- CRQ(startpoint)
  names(startcost) <- "Co"
  feasible0 <- ( abs( S2con(startpoint) ) <= 0.5E-5 )
  names(feasible0) <- "S2"

  res <- Rsolnp::solnp(pars=startpoint, fun=CRQ, eqfun=S2con, eqB=c(0), control=list(trace=0))

  finalpoint <- res$pars
  names(finalpoint) <- c("R*", "Q*")
  finalcost  <- res$values[length(res$values)]
  names(finalcost) <- "C*"
  feasible1 <- ( abs( S2con(finalpoint) ) <= 0.5E-5 )
  names(feasible1) <- "S2"

  opt <-  (res$converge == 0)

  sol <- list( Startpoint = startpoint, StartFeasibility=feasible0, Startcost = startcost,
               Finalpoint = finalpoint, FinalFeasibility=feasible1, Finalcost = finalcost,
               Optimum = opt)

  return(sol)
}
