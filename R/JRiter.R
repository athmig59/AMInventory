#' @title Joint Replenishment - Iterative Approach
#'
#' @author Athanasios Migdalas, \email{athmig59@gmail.com}
#'
#' @description
#' A deterministic joint replenishment model is considered in which
#' there is a joint set up cost independent of the number of items ordered 
#' olus individual set up costs for each item. Cycle times (and, consequently, 
#' batch quantities) are assumed constant. No backorders are allowed and the
#' production time is being disregarded. The problem is to determine the cycle 
#' times that minimize the sum of holding and setup costs.
#' See section 7.3  in Axsäter's book for more information.
#' 
#' The iterative approach (see section 7.3.1.1 in Axsäter) is a heuristic algorithm
#' which attempts to determine cycle times as integer multiples of the smallest 
#' cycle time, initially approximated with the EOQ lot size formula.
#' 
#' @param{d} {demand per time unit for each item}
#' @param{A} {set up cost for the group}
#' @param{h} {holding cost per item unit and time unit}
#' @param{a} {set up cost per item}
#' 
#' @return A list containing
#'  \itemize{
#'  \item{BasicPeriod} {corresponds to cycle time for the "first" item}
#'  \item{Multipliers} {of the basic period to express cycle 
#'    times for the rest of the items}
#'  \item{Ordered} {item indeces ordered}
#'  \item{LBD} {cost lower bound}
#'  \item{UBD} {cost upper bound}
#'  }
#' 
#' @references
#'  - Axsäter, S. (2015) Inventory Control, 3d Edition, Springer
#'
#' @seealso JRroundy
#' 
#' @examples 
#' ## ** Example based on Example 7.2 in Axsäter **
#' ## Four item with joint set up cost $300 and equal individual set up
#' ## and holding costs of $50 and $10 respectively. The demand for the 
#' ## four items are 5000, 1000, 700 and 100 respectively.
#' ##
#' JRiter(d=c(5000,1000,700,100),A=300,h=c(10,10,10,10),a=c(50,50,50,50))
#' 
#' 
#' @export
JRiter <- function(d,A,h,a){
  
  N <- length(d)
  
  Aa <- A+a
  
  if( sum(Aa > 0 ) != N){return(NULL)}
  
  eta <- h*d
  
  if( sum(eta > 0) != N){return(NULL)}
  
  
  aeta <- a/eta
  
  ordered <- order(aeta)
  
  n  <- vector(mode="numeric", length=N)
  nn <- as.vector(n)
  
  one <- ordered[1]
  
  n[one] <- 1
  
  for(i in ordered[2:N]){
    n[i] = round( sqrt(aeta[i] * eta[one]/Aa[one] ) )
  } 
  
  # print(n[ordered])
  
  repeat{
    
    T1 <- sqrt( 2*(A+a[one] + sum(a[ordered[2:N]]/n[ordered[2:N]]))/(eta[one]+sum(eta[ordered[2:N]]*n[ordered[2:N]])) )
    
    C <- sqrt( 2*(A+a[one] + sum(a[ordered[2:N]]/n[ordered[2:N]]))*(eta[one]+sum(eta[ordered[2:N]]*n[ordered[2:N]])))
    
    # print(T1)
    # print(C)
    
    nn[one] <- 1
    
    for(i in ordered[2:N]){
      g <- 2*a[i]/(eta[i]*T1**2)
      nn[i] <-ceiling((-1 + sqrt(1+4*g))/2) 
    }
    
    # print(nn)
    
    if( sum(nn == n) == N ){
      break
    }
    else
    {
      n <- nn
    }
  }
  
  Clbd <- sqrt(2*(A+a[one])*eta[one])+sum(sqrt(2*a[ordered[2:N]]*eta[ordered[2:N]]))
  
  #print(Clbd)
  
  return( list(BasicPeriod=T1, Multipliers=n, Ordered=ordered, LBD=Clbd, UBD=C))
}
