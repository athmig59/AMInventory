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
#' The iterative method of alternatingly solving the first order optimality
#' equations is applied in order to minimize the expected total cost.
#' The function "uniroot" of the package "stats" is used in order to solve
#' the equation for R.
#'
#' For more details on the iterative approach see section 6.1.2 in Axsäter's book.
#'
#' @references
#'  - Axsäter, S. (2015) Inventory Control, 3d Edition, Springer
#'
#' @seealso JointRQb1, JointRQu, JointRQss, JointRQhw, JointRQb1S2
#'
#' @param mu:      mean value of demand
#' @param sigma:   standard deviation of demand
#' @param L:       lead-time
#' @param h:       holding cost
#' @param b1:      stockout cost
#' @param A:       fix order cost
#'
#' @return A data frame containing the results of the iterations in three columns:
#' \itemize{
#'    \item{Q:}  {the successive order quantities generated}
#'    \item{R:}   {the successive reorder points generated}
#'    \item{C:}  {the corresponding total cost for each successive pair of Q and R}
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
#' # No initial guesses are required since the algorithm will start by utilizing
#' # the EOQ produced order quantity.}
#'
#' JointRQb1iter(mu=50, sigma=20, L=4, h=2, b1=20, A=100)
#'
#' @export
JointRQb1iter <- function(mu, sigma, L, h, b1, A){

  # Tolerances

  tol0 <- 0.5e-4 # to stop iterations
  tol1 <- 0.5e-7 # for the equation solver

  # Lead-time parameters

  Lmu <- L * mu
  Lsigma <- sigma * sqrt(L)

  # G function (Loss function)

  G <- function(x){
    return( dnorm(x) - x * (1-pnorm(x)) )
  }

  # H function

  H <- function(x){
    return( 0.5 * ( (x**2+1) * (1-pnorm(x)) - x * dnorm(x) ) )
  }

  # objective function

  CRQ <- function(R,Q){

    z1 <- ( R - Lmu ) / Lsigma
    z2 <- ( R + Q - Lmu ) / Lsigma

    H1 <- H(z1)
    H2 <- H(z2)

    val <-  h * ( R + 0.5 * Q - Lmu ) + ( h + b1 ) * Lsigma**2 * ( H1 - H2 ) / Q + A * mu / Q

    return( val )
  }

  # partial derivative w.r.t. R

  dCR <- function(R){

    z2 <- ( R - Lmu ) / Lsigma
    z1 <- ( R + Qo - Lmu ) / Lsigma

    G1 <- G(z1)
    G2 <- G(z2)

    val <- h + ( h + b1 ) * Lsigma * ( G1 - G2 ) / Qo

    return(val)

  }

  # Economic order quantity

  EOQ <- function(A, mu, h){
    return( sqrt( 2 * A * mu / h ) )
  }

  # Equation to update Q

  Qnew <- function(Qold, Rold){

    z1 <- ( Rold - Lmu ) / Lsigma
    z2 <- ( Rold + Qold - Lmu ) / Lsigma

    H1 <- H(z1)
    H2 <- H(z2)
    G2 <- G(z2)

    q <- sqrt( 2*A*mu/h + 2*(h+b1)*Lsigma**2 * ( H1-H2 - Qold*G2/Lsigma )/h )

    return(q)
  }

  # determine a bracketing interval

  bracket <- function(a,b,stp){

    while(TRUE){
      if( sign( dCR(a) )  != sign( dCR(b) ) ){
        rint <- c(a, b)
        break
      }else{
        a <- a - stp
        b <- b + stp
      }
    }

    return(rint)
  }

  # initial Q from EOQ

  Qo <- EOQ(A, mu, h)
  Ro <- 0

  a <- -mu; b <- mu; stp <- Qo

  Qs <- c(Qo)
  Rs <- c()
  Cs <- c()

  while(TRUE){

    rint <- bracket(a,b,stp)

    sol <- stats::uniroot( dCR, rint, tol=tol1 )
    R1  <- sol$root
    Rs  <- append(Rs, R1)

    # curve(dCR, from=rint[1], to=rint[2])
    # abline(h=0)
    # abline(v=R1)

    C1 <- CRQ(R1,Qo)
    Cs <- append(Cs, C1)

    if( abs(R1 - Ro) <= tol0 ) { break }

    Q1 <- Qnew(Qo,R1)

    if( abs(Q1 - Qo) <= tol0 ) { break }
    Qs <- append(Qs, Q1)

    Ro   <- R1
    Qo   <- Q1
    stp  <- Q1
  }

  df <- data.frame(Q=Qs, R=Rs, C=Cs)

  return(df)

}
