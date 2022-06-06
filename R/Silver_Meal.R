#' @title Siver-Meal heuristic for lot-sizing with time-varying demand
#'
#' @author Athanasios Migdalas, \email{athmig59@gmail.com}
#'
#' @description
#' This model assumes that there is a finite number of periods and the demand for
#' each such period is known and realizable at the beginning of the period. It is
#' further assumed that there is no initial stock. The delivery of an entire order,
#' when it happens, deliverss the entire order quantity. Holding and ordering costs
#' are constant over time. No back orders are allowed.
#'
#' The Silver-Meal algorithm is a heuristic alternative to the optimizing dynamic programming
#' algorithm by Wagner-Whitin. Although it can miss the optimal solution, it does often
#' produces it, however, since it does not verify optimality, you cannot be asured of it.
#' See sections 4.5, 4.6, and 4.7 in Axsäter's book.
#'
#' @references
#'  - Axsäter, S. (2015) Inventory Control, 3d Edition, Springer
#'  - Snyder, L. V. and Shen, Z.-J., M. (2011) Fundamentals of Supply Chain Theory, Wiley
#'
#' @seealso   Wagner_Whitin, Balance_Costs
#'
#' @param {d} {vector specifying the demands of the periods}
#' @param {A} {Fixed ordering cost}
#' @param {v} {Unit purchase price}
#' @param {r} {Interest rate for holding a unit as a percentage of purchase price. Default value is 0}
#'
#' @returns A list consisting of
#'    \itemize{
#'    \item{\code{Costs}} {vector containing three elements:
#'       \itemize{
#'          \item{\code{TotalCost}} {the total policy cost}
#'          \item{\code{TotalOrderCost}} {additive component of total cost}
#'          \item{\code{TotalHoldingCost}} {additive component of total cost}}}
#'     \item{\code{OrderPeriod}} {policy vector indicating for each period demand the order period}}
#'
#' @examples
#' #' ## ** Example 1 (based on Example 4.5 in Axsäter) **
#' ## The ordering cost is $300 and the holding cost $1 per unit and period.
#' ## The horizon is 10 periods long and the demands for the 10 periods are:
#' ## 50, 60, 90, 70, 30, 100, 60, 40, 80, 20
#' ##
#'
#' Silver_Meal(d=c(50, 60, 90, 70, 30, 100, 60, 40, 80, 20), A=300, h=1)
#'
#' ## ** Example 2 (based on Example 3.9 in Snyder-Shen) **
#' ## Fixed ordering cost is $500, the holding cost $2 per unit and period.
#' ## There are four periods with corresponding demand 90, 120, 80 and 70.
#' ##
#'
#' Silver_Meal(d=c(90, 120, 80, 70), A=500, h=2)
#'
#' ## ** Example 3 **
#' ## As in Example 2 but the unit purchase price is $20 and the
#' ## interest rate defining the holding cost per unit and period is 0.1.
#' ##
#'
#' Silver_Meal(d=c(90, 120, 80, 70), A=500, v=20, r=0.1)
#'
#'
#' @export
Silver_Meal <- function(d,A,h=0,v=0, r=0) {
  L  <- length(d)

  if( r > 0 ){
    h <- h + r*v
  }

  in_period <- 1:L
  cost <- A
  op <- 1    # Order period
  p <- 2     # Trial period
  k <- 2     # Holding periods counter

  total_cost <- 0
  total_order_cost <- A
  total_holding_cost <- 0

  repeat{
    holding_cost <- (k-1)*h*d[p]
    trial <- cost + holding_cost
    if( trial/k <= cost/(k-1) ){
      cost <- trial
      total_holding_cost <- total_holding_cost + holding_cost
      in_period[p] <- op
      k <- k + 1
    }
    else
    {
      total_cost <- total_cost + cost
      cost <- A
      op <- p
      k <- 2
      total_order_cost <- total_order_cost + A
    }
    p <- p + 1
    if(p > L) {
      total_cost <- total_cost + cost
      break}
  }
  costs <- c(total_cost,total_order_cost,total_holding_cost)
  names(costs) <- c("TotalCost", "TotalOrderCost", "TotalHoldingCost")
  names(in_period) <- 1:L
  sol <- list(costs,in_period)
  names(sol) <-c("Costs","OrderPeriod")
  return(sol)
}
