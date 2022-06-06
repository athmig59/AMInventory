#' @title Joint Replenishment - Roundy's 98% Approximation
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
#' Roundy's 98% approximation algorithm is a Lagrangian relaxation approach which
#' produces cycle times that are power of two multiples of a basic cycle time. The
#' approach guarantees that the obtained solution has a cost which is at most 2% higher
#' than the optimal cost. See section 7.3.1.2 in Axsäter for more information.
#'
#' @param{d} {demand per time unit for each item}
#' @param{A} {set up cost for the group}
#' @param{h} {holding cost per item unit and time unit}
#' @param{a} {set up cost per item}
#' @param{q}  {the basic period to be used in power two cycle times, default 1}
#'
#' @return A list containing
#' \itemize{
#' \item{Table} {a data frame containing
#'  \itemize{
#'  \item{Product} {numbers the products}
#'  \item{Order} {indices the products in increasing order of η[i]}
#'  \item{Group} {assigns products to groups, group 0 being the aggregate group}
#'  \item{CycleTimes} {are NOT of power of two}
#'  \item{Power} {lists the exponents of 2 for powe of two cycle times}
#'  \item{PowerTwoTimes} {lists the power of two cycle times}
#'  }}
#'  \item{LBD} {a Lagrangian lower bound on the optimal cost}
#'  \item{UBD} {an upper bound on the optimal cost}
#'  }
#'
#' @seealso JRiter
#'
#'
#' @examples
#' ## ** Example based on Example 7.3 in Axsäter **
#' ## Four item with joint set up cost $300 and equal individual set up
#' ## and holding costs of $50 and $10 respectively. The demand for the
#' ## four items are 5000, 1000, 700 and 100 respectively.
#' ##
#' JRroundy(d=c(5000,1000,700,100),A=300,h=c(10,10,10,10),a=c(50,50,50,50))
#' ##
#' ##
#' JRroundy(d=c(5000,1000,700,100),A=300,h=c(10,10,10,10),a=c(50,50,50,50), q=1.88)
#'
#' @export
JRroundy <- function(d,A,h,a,q=1){

  totcost <- function(T){
    tc <- A/T[1]
    for(i in 1:N) {
      tc <- tc + 0.5 * eta[i] * T[i] + a[i] / T[i]
    }
    return(tc)
  }

  if(q<=0){return(NULL)}

  N <- length(d)

  eta <- h*d
  if(sum(eta > 0) != N){return(NULL)}

  aeta <- a/eta

  ordered <- order(aeta)

  group <- c(1:N)

  k <- 1
  s1 <- A + a[ordered[k]]
  s2 <- eta[ordered[k]]
  group[ordered[k]] <- 0
  while(s1/s2 > a[ordered[k+1]]/eta[ordered[k+1]]){
    s1 <- s1 + a[ordered[k+1]]
    s2 <- s2 + eta[ordered[k+1]]
    group[ordered[k+1]] <- 0
    k <- k+1
    if(k == N){break}
  }

  cT <- vector(mode="numeric", length=N)
  cT[ordered[1:k]] <- sqrt(2*s1/s2)
  for(i in k+1:N){
    cT[ordered[i]] <- sqrt(2*a[ordered[i]]/eta[ordered[i]])
  }

  LBD <-totcost(cT)

  pot <- vector(mode="numeric", length=N)
  for(i in 1:N){
    lg <- log2(cT[i]*(1/q))
    if(lg > 0) {
      pot[i] <- ceiling(lg)
    }
    else
    {
      pot[i] <- round(lg)
    }
  }

  rT <- 2**pot*q

  UBD <- totcost(rT)

  Solution <-data.frame(Product=1:N, Order=ordered, Group=group, CycleTime=cT, Power=pot, PowerTwoTime=rT)
  return(list(Table=Solution, LBD=LBD, UBD=UBD))
}
