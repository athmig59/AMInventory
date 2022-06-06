#' @title Reorder point computation for given service level type 1
#'
#' @author Athanasios Migdalas, \email{athmig59@gmail.com}
#'
#' @description
#' Stochastic, single echelon system is assumed with normally distributed demand.
#' An exact value of the reorder point is returned in accordance with Chapter 5 in Axsäter's book.
#' See sections 5.6 and 5.7.2
#'
#' @references
#'  - Axsäter, S. (2015) Inventory Control, 3d Edition, Springer
#'
#' @seealso Plot_R_vs_S1, Reorder_Point_S2
#' 
#' @param {S1} {numeric, service level of type 1}
#' @param {sigma} {numeric, standard deviation of lead-time demand}
#' @param {mu} {numeric, mean of lead-time demand}
#' 
#' @return A list containing:
#'   \itemize{
#'   \item{\code{R}} {the exact value of the re-order point} 
#'   \item{\code{SS}} {the safety stock}
#'   \item{\code{z}} {the safety factor}}
#' 
#' @examples 
#' ## ** Example 1 based on Example 5.4 in Axsäter **
#' ## For service level S1 = 0.376, lead-time demand with μ=10 and σ=sqrt(10),
#' ## the corresponding reorder point, R, is sought.
#' ##
#' Reorder_Point_S1(S1=0.376, sigma=sqrt(10), mu=10)
#'
#' ## ** Example 2 **
#' ## Lead-time demand as above but for sevice level S1 of 75%.
#' ## 
#' Reorder_Point_S1(S1=0.75, sigma=sqrt(10), mu=10)
#' 
#' @export
Reorder_Point_S1 <- function(S1, sigma, mu){
  z <- stats::qnorm(S1)
  SS <- z * sigma
  R <- mu + SS
  return( list(R=R, SS=SS, z=z) )
} 