#' @title Common Cycle Heuristic for the Economic Lot Scheduling Problem (ELSP)
#'
#' @author Athanasios Migdalas, \email{athmig59@gmail.com}
#'
#' @description
#' The Economic Lot Scheduling Problem (ELSP) is concerned with the determination of
#' cyclic schedules for a number of items with constant demands attempting to achieve
#' smooth production. Backorders are not allowed. The production rate is finite.
#' The goal is to minimize standard holding and ordering costs. See section 7.2 in
#' Axsäter's book for further information.
#'
#' The common cycle heuristic generates a feasible solution and an upper bound to
#' the problem, provided it exists. The algorithm is discussed in section 7.2.1.3 in
#' Axsäter's book.
#'
#' @param{d} {vector of demand rates}
#' @param{A} {vector of fixed order/set up costs}
#' @param{h} {vector of holding cost per item per time unit}
#' @param{p} {vector of production rates}
#' @param{s} {vector of set up times}
#'
#' @return A vector with elments:
#'   \itemize{
#'    \item{Tmax} {the approximation (corresponds to Topt
#'     in Axsäter's book)}
#'    \item{Tmin} {the lower bound on the cycle time}
#'    \item{That} {corresponds to equation (7.18) in Axsäter's book}
#'     \item{LBD}  {a lower bound on the sought optimal cost (computed
#'     by the independent solution procedure \code{ELSPihs}}
#'    \item{UBD} {an upper bound on the sought cost, corresponding
#'     to Tmax}
#'   }
#'   or NULL, if inappropriate indata is supplied.
#'
#' @details
#'  Requires AMInventory::ELSPish()
#'
#' @references
#'  - Axsäter, S. (2015) Inventory Control, 3d Edition, Springer
#'
#' @seealso ELSPihs, ELSPsh
#'
#' @examples
#'  ## ** Example 1 **
#'  ## Data from Problem 7.3 in Axsäter
#'  ## A=c(800,500,1000)
#'  ## d=c(48, 20, 32)
#'  ## h=c(0.060, 0.040, 0.048)
#'  ## p=c(200, 100, 100)
#'  ## s=c(0.50, 0.25, 1.00)
#'  ##
#'  ELSPcch(d,A,h,p,s)
#'
#'  ## ** Example 2 **
#'  ## ** Bomberger's example (see Example 7.1 in Axsäter) **
#'  ## d=c(400,400,800,1600,80,80,24,340,340,400)
#'  ## A=c(15,20,30,10,110,50,310,130,200,5)
#'  ## h=c(0.2708, 7.396, 5.313, 4.167, 116, 11.15, 62.50, 245.8, 37.5, 1.667)/10**5
#'  ## p=c(30000, 8000, 9500, 7500, 2000, 6000, 2400, 1300, 2000, 15000)
#'  ## s=c(0.125, 0.125, 0.25, 0.125, 0.5, 0.25, 1, 0.5, 0.75, 0.125)
#'  ##
#'  ELSPcch(d,A,h,p,s)
#'
#' @export
ELSPcch <- function(d,A,h,p,s){

  N <- length(d)

  if( sum(p > 0, na.rm=TRUE) != N){return(NULL)}
  if( sum(p > d, na.rm=TRUE) != N){return(NULL)}

  totSigma <- 0
  totRho <- 0
  totA <- 0
  totHDR <- 0

  for(i in 1:N){
    rho <- d[i]/p[i]
    totSigma <- totSigma + s[i]
    totRho <- totRho + rho
    totA <- totA + A[i]
    totHDR <- totHDR + h[i]*d[i]*(1-rho)
  }

  Tmin <- NULL
  if (totRho < 1){
    Tmin <- totSigma/(1-totRho)
  }

  That <- sqrt(2*totA/totHDR)

  Tmax <- max(Tmin,That,na.rm=TRUE)

  cost <- function(cT){
    if(cT <= 0){return(NULL)}
    return( sum(A/cT+0.5*h*d*(1-d/p)*cT) )
  }

  Cubd <- cost(Tmax) # Upper bound

  sol <- AMInventory::ELSPish(d,A,h,p,s)

  Clbd <- sol$TotalCost

  return( c(Tmax=Tmax, Tmin=Tmin, That=That, LBD=Clbd,
            UBD=Cubd))
}

