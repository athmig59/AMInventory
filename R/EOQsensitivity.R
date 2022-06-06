#' @title Sensitivity of the Economic Order Quantity
#'
#' @author Athanasios Migdalas, \email{athmig59@gmail.com}
#'
#' @description
#' Computes the quotient (C/C*) of the cost (C) for a non-optimal versus
#' the cost (C*) for the optimal order quantity.
#'
#' @references
#'  - Axsäter, S. (2015) Inventory Control, 3d Edition, Springer
#'
#' @seealso   EOQ, DEOQ, IDEOQ, EOQ_visualize_costs
#'
#' @param {Qopt} {the optimal order quantity}
#' @param {Qother} {a tentative order quantity}
#'
#' @returns A numerical value corresponding to the cost quotient, or NULL, in
#' case of inadequate input.
#'
#' @examples
#' EOQsensitivity(Qopt=200,Qother=300)
#'


EOQsensitivity <- function(Qopt,Qother){
  quotient <- NULL
  if(Qopt > 0 && Qother > 0){
    quotient <- 0.5*(Qother/Qopt + Qopt/Qother)
  }
  return(quotient)
}
