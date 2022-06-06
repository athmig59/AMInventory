#' @title Deterministic Economic Order Quantity Model
#'
#' @author Athanasios Migdalas, \email{athmig59@gmail.com}
#'
#' @description
#' It implements three single period, single echelon, economic order quantity (EOQ)
#' model variants:
#' the classic, with production,and with planned backorders.
#' See chapter 4 in Axsäter's book for the details
#'
#' @references
#'  - Axsäter, S. (2015) Inventory Control, 3d Edition, Springer
#'
#' @seealso   EOQsensitivity, DEOQ, IDEOQ, EOQ_visualize_costs
#'
#' @param {d} {demand rate}
#' @param {A} {fixed order cost}
#' @param {h} {holding cost}
#' @param {p} {production rate, default value 0}
#' @param {b1} {backorder cost, default value 0}
#'
#' @details
#' The appropriate model is automatically selected based on the user input.
#' The default model is the classic EOQ and corresponds to the absence (i.e., zero values)
#' of production rate and backorder cost.
#'
#' @returns A list containing:
#' - {Q:} {the optimal order quantity}
#' - {C:} {the optimal cost}
#' - {x:} {is returned only in the case of the backorder model, the fraction of the demand being backorderd }
#' - or NULL values, in case of inadequate input.
#'
#' @examples
#' \dontrun{** Example (based on Examples 4.1 and 4.3 in Axsäter) **
#' An item costs $100 per unit and has a holding cost rate per year corresponding to 20% of
#' this value. The constant demand is 300 units per year. The ordering cost is $200. Clearly,
#' h= 0.20(100) = $20 per year.}
#'
#' EOQ(d=300,A=200,h=20)
#'
#' \dontrun{Assume that bacorders are allowed with a backorder cost $100 per unit and year.
#'  What would then be the order quantity? What the fraction of bacorded demand?}
#'
#'  EOQ(d=300,A=200,h=20, b1=100)
#'
#' \dontrun{Assume instead that the item is produced at a rate of 330 units per year.
#'  What would then be the order quantity?}
#'
#' EOQ(d=300,A=200,h=20, p=330)
#'
#' @export
EOQ <- function(d, A, h, p=0, b1=0){

  # EOQ

  if(p == 0 && b1 == 0){
    if(h <= 0){
      Q <- NULL
      C <- 0
    } else{
      Q <- sqrt(2*A*d/h)
      C <- sqrt(2*A*d*h)
    }

    eoq.sol <- list(Q=Q,C=C)
    return(eoq.sol)
  }

  # BEOQ

  if(b1 > 0){
    x <- h/(h+b1)
    hb <- h*b1
    Q <- NULL
    if (hb > 0){
      Q <- sqrt(2*A*d*(h+b1)/hb)
    }
    C <- sqrt((2*A*d*h*b1)/(h+b1))

    eoq.sol <- list(Q=Q,C=C,x=x)
    return(eoq.sol)
  }

  # PEOQ

  if(p > 0) {
    Q <- NULL
    C <- NULL
    if (p > d){
      hdp <- h*(1-d/p)
      if( hdp > 0 ){
        Q <- sqrt( 2*A*d/hdp)
      }
      C <-0.5*Q*hdp + d*A/Q
    }

    eoq.sol <- list(Q=Q,C=C)
    return(eoq.sol)
  }
}

