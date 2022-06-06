#' @title Simple Heuristic approach for the Economic Lot Scheduling Problem (ELSP)
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
#' The simple heuristic, or Doll-Whybark heuristic, generates successively improved
#' multipliers of a basic period. The algorithm is discussed in section 7.2.1.5 in
#' Axsäter's book.
#'
#' CAUTION: The solution is not necessarily feasible and re-entrance with modified
#' multipliers may be necessary.
#'
#' @seealso ELSPish, ELSPcch
#'
#' @details
#'  Requires AMInventory::ELSPihs()
#'
#' @param{d} {vector of demand rates}
#' @param{A} {vector of fixed order/set up costs}
#' @param{h} {vector of holding cost per item per time unit}
#' @param{p} {vector of production rates}
#' @param{s} {vector of set up times}
#' @param{n} {if given it should be vector of multipliers, otherwise it is initialized using the
#'            independent solution procedure ELSPihs}
#' @param{W} {if given it should be a basic period, otherwise it is initialized using the
#'            independent solution procedure ELSPihs}
#'
#'
#' @return A list with elements:
#' \itemize{
#' \item{Table} {is a data frame with two columns:
#'   \itemize{
#'    \item{Multiplier} {listing individual multipliers of the basic period}
#'    \item{ProductionTime} {listing individual production times}
#'  }}
#'  \item{BasicPeriod}  {the computed basic period}
#'  \item{TotalCost} {the total cost for the obtained multipliers and basic period}
#'  \item{Iterations} {the number of iterations required to convergence}
#'  }
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
#'  ELSPsh(d,A,h,p,s)
#'
#'  ## ** Example 2 **
#'  ## ** Bomberger's example (see Example 7.1 in Axsäter) **
#'  ## d=c(400,400,800,1600,80,80,24,340,340,400)
#'  ## A=c(15,20,30,10,110,50,310,130,200,5)
#'  ## h=c(0.2708, 7.396, 5.313, 4.167, 116, 11.15, 62.50, 245.8, 37.5, 1.667)/10**5
#'  ## p=c(30000, 8000, 9500, 7500, 2000, 6000, 2400, 1300, 2000, 15000)
#'  ## s=c(0.125, 0.125, 0.25, 0.125, 0.5, 0.25, 1, 0.5, 0.75, 0.125)
#'  ##
#'  res <- ELSPsh(d,A,h,p,s)
#'  ##
#'  ## Rerun with
#'  ## n <- res$Table$Multiplier
#'  ## n[9] <- 2
#'  ## W <- res$BasicPeriod
#'  ##
#'  ELSPsh(d,A,h,p,s, n=n, W=W)
#'  ##
#'  ## To get Table 7.3 in Axsäter
#'  ##
#' @export
ELSPsh <- function(d,A,h,p,s,n=NULL, W=NULL){

  totcostfn <- function(W,n){
    l <- length(n)
    tc <- 0.0
    for(i in 1:l){
      tc <- tc + A[i]/(n[i]*W) + h[i]*d[i]*(1-d[i]/p[i])*n[i]*0.5*W
    }
    return(tc)
  }

  if(is.null(n) ){

    sol <- AMInventory::ELSPish(d,A,h,p,s)
    if(is.null(sol)){return(NULL)}

    L <- length(sol$Solution$CycleTime)

    n <- vector(mode="numeric", length=L)

    W <- min(sol$Solution$CycleTime,na.rm=TRUE)

    for(i in 1:L){
      m <- round( log2( sol$Solution$CycleTime[i]/W ) )
      n[i] <- 2**m
    }
  }else{
    L <- length(n)
  }

  nn <- vector(mode="numeric", length=L)
  totA <- sum(A/n, na.rm=TRUE)

  iter <- 1

  repeat{
    cT <- W*n
    totHDRN <- sum(h*d*(1-d/p)*n)
    W <- sqrt(2*totA/totHDRN)

    for(i in 1:L){
      #m <- max(0, round( log2( sol$Solution$CycleTime[i]/W ) ) )
      m <- round( log2( cT[i]/W ) )
      nn[i] <- 2**m
    }

    if( sum(nn==n,na.rm=TRUE) == L){
      cT <- W*n
      break
    }else
    {
      n <- nn
      totA <- sum(A/n, na.rm=TRUE)
      iter <- iter + 1
    }
  }

  sigma <- s + cT*d/p

  table <- data.frame(Multiplier=n, ProductionTime=sigma)
  return( list(Table=table, BasicPeriod=W, TotalCost=totcostfn(W,n), Iterations=iter) )
}
