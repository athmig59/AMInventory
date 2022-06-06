#' @title Independent Solution Heuristic for the Economic Lot Scheduling Problem (ELSP)
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
#' The independent solution heuristic generates a lower bound to
#' the problem. The algorithm is discussed in section 7.2.1.2 in
#' Axsäter's book.
#'
#' NOTE: See ELSPcsh for the generation of an upper bound.
#'
#' CAUTION: The produced solution is not, in general, feasible
#'
#' @param{d} {vector of demand rates}
#' @param{A} {vector of fixed order/set up costs}
#' @param{h} {vector of holding cost per item per time unit}
#' @param{p} {vector of production rates}
#' @param{s} {vector of set up times}
#'
#' @return A list with elements:
#' \itemize{
#' \item{Solution} {is a data frame with three columns:
#'   \itemize{
#'    \item{CycleTime} {listing individual cycle times}
#'    \item{Cost} {listing individual costs}
#'    \item{ProductionTime} {listing individual production times}
#'  }}
#'  \item{TotalCost}  {is a lower bound on the
#'     sought cost}
#'  }
#'
#' @details
#'  Used by AMInventory::ELSPcch() and AMInventory::ELSPsh()
#'
#' @references
#'  - Axsäter, S. (2015) Inventory Control, 3d Edition, Springer
#'
#' @seealso ELSPcch, ELSPsh
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
#'  ELSPish(d,A,h,p,s)
#'
#'  ## ** Example 2 **
#'  ## ** Bomberger's example (see Example 7.1 in Axsäter) **
#'  ## d=c(400,400,800,1600,80,80,24,340,340,400)
#'  ## A=c(15,20,30,10,110,50,310,130,200,5)
#'  ## h=c(0.2708, 7.396, 5.313, 4.167, 116, 11.15, 62.50, 245.8, 37.5, 1.667)/10**5
#'  ## p=c(30000, 8000, 9500, 7500, 2000, 6000, 2400, 1300, 2000, 15000)
#'  ## s=c(0.125, 0.125, 0.25, 0.125, 0.5, 0.25, 1, 0.5, 0.75, 0.125)
#'  ##
#'  ELSPish(d,A,h,p,s)
#'
#' @export
ELSPish <- function(d,A,h,p,s){

  N <- length(d)

  if( sum(p > 0, na.rm=TRUE) != N){return(NULL)}
  if( sum(p > d, na.rm=TRUE) != N){return(NULL)}

  cT <- vector(mode="numeric",length=N)
  C  <- as.vector(cT)
  sigma <- as.vector(cT)

  totalC <- 0
  for(i in 1:N){
    rho <- d[i]/p[i]
    cT[i] <- sqrt( 2 * A[i] / ( h[i] * d[i] * ( 1 - rho ) ) ) # time
    C[i]  <- sqrt( 2 * A[i] * h[i] * d[i] * ( 1 - rho ) )     # cost
    sigma[i] <- s[i] + rho * cT[i]                            # time
    totalC <- totalC + C[i]
  }

  solution <- data.frame(CycleTime=cT, Cost=C, ProductionTime=sigma)

  return(  list(Solution=solution, TotalCost=totalC) )
}
