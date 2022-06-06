#' @title EOQ model with order quantity (incrimental) discounts
#'
#' @author Athanasios Migdalas, \email{athmig59@gmail.com}
#'
#' @description
#' This model is based on the assumption that the supplier offers
#' discounts based on the quantity ordered. The incremental discount
#' case is considered here, that is, given the quantity break points for
#' the different prices, the units within an interval with end points defined by
#' consecutive breakpoints incur the purchase price for that interval.
#'
#' @references
#'  - Axsäter, S. (2015) Inventory Control, 3d Edition, Springer
#'  - Snyder, L. V. and Shen, Z.-J., M. (2011) Fundamentals of Supply Chain Theory, Wiley
#'
#' @seealso   EOQsensitivity, EOQ_visualize_costs, EOQdisc, EOQ
#'
#' @param {d} {Yearly demand rate}
#' @param {A} {Fixed ordering cost}
#' @param {Qb} {A vector of quantity break points}
#' @param {v} {A vector of unit purchase price corresponding to quantity break points}
#' @param {r} {Interest rate for holding a unit as a percentage of purchase price}
#'
#' @returns A list consisting of two elements:
#'     \itemize{
#'     \item{Sol} {a vector containing the optimal solution in two elements:
#'              \itemize{
#'              \item{optQ} {the optimal quantity that minimizes total cost}
#'              \item{optC} {the optimal total cost}}}
#'      \item{Table} {a data frame with three columns containing intermediate result:
#'              \itemize{
#'              \item{Q} {quantities}
#'              \item{C} {corresponding costs}
#'              \item{Feas} {provides indication of feasibility}}} }
#'
#' @examples
#' ## ** Example 1 (based on Axsäter's Example 4.2) **
#' ## Suppose a supplier offers an item for a unit price of $100 unless the order quantity is 100
#' ## or more units in which case the unit price is $95. Assume that the interest rate is 0.2.
#' ## The yearly demand for the item is constant 300 units. The order cost is $200.
#'
#' EOQincr(d=300, A=200, Qb=c(100), r=0.2, v=c(100,95) )
#'
#' ## ** Example 2 ( based on Snyder-Shen's Example 3.7) **
#' ## A candy with yearly demand of 1300 units is ordered for a fix order cost of $8
#' ## from a supplier who offers it for a unit price of $0.75. However,  if 400 or more up
#' ## to 799 are purchased, the unit price is discounted to 0.72; and if 800 or more are
#' ## purchased, the price is further discounted to $0.68 per unit. The interest rate is 0.3
#' ## of the purchased price per unit per year.
#'
#' EOQincr(d=1300, A=8, Qb=c(400,800),r=0.3, v=c(0.75,0.72,0.68) )
#'
#' @export
EOQincr <- function(d, A, Qb, v, r){

  # Calculate holding costs
  h <- vector(mode='numeric',length=length(v))
  for( i in seq(1,length(v)) ){
    h[i] = r * v[i]
  }

  # Breakpoints
  Qs <- Qb
  Qs <- c(0,Qs)
  Qs <- c(Qs,Inf)

  # costs based on Incremental discounts

  idc <- vector(mode='numeric', length=length(v))
  idc[1] <- 0
  for( i in 2:length(idc)){
    idc[i] <- idc[i-1] + v[i-1]*(Qs[i]-Qs[i-1])
  }
  for ( i in 1:length(idc)){
    idc[i] <- idc[i] - v[i] * Qs[i]
  }

  # holding costs based on Incremental discounts

  idh <- r * idc

  # Trial order quantities based on Incremental discounts

  Q <- vector(mode='numeric', length=length(idc))
  for( i in seq(1,length(idc))){
    Q[i] <- sqrt( (2*(A + idc[i])*d) / h[i] )
  }

  # Compute corresponding costs

  C <- vector(mode='logical', length=length(idc))
  for( i in seq(1,length(idc))){
    C[i] <- v[i]*d + 0.5 * idh[i] + sqrt( 2*(A + idc[i])*d*h[i])
  }

  # Determine feasibility and optimality

  f <-  vector(mode='numeric', length=length(idc))
  f[1:length(f)] <- FALSE
  fcount <- 0
  for( i in seq(1,length(Q))){
    if( Qs[i]<= Q[i] && Q[i] < Qs[i+1] ) {
      f[i] <- TRUE
      fcount <- fcount + 1
    }
  }

  optC <- NULL
  optQ <- NULL

  if( fcount > 0){
    optC <- min(C[f!=0])
    optQ <- Q[match( optC,C)]
  }

  # Prepare output

  Table <- data.frame(Q=Q, C=C, Feas=f)
  sol <- list(Sol=c(optQ=optQ, optC=optC), Table=Table)

  return(sol)
}
