#' @title The classic Newsvendor model
#'
#' @author Athanasios Migdalas, \email{athmig59@gmail.com}
#'
#' @description
#' The newsvendor or newsboy model is a classical sigle period inventory problem with
#' stochastic demand that is assumed to be normally distributed with mean μ and standard deviation σ.
#' It is sought to determine how much to order before the period starts. There are penalties
#' associated with ordering too much (overage cost) and ordering too little (underage cost)
#' compared to the actual demand, revealed after the start of the period.
#'
#' @details
#' There are two alternative ways of specifying the problem:
#' \itemize{
#' \item{Alternative I:} {Specify
#'    \itemize{
#'    \item{Co:} {Overage cost, i.e.,  cost per unit for unsold items (c.f. holding cost)}
#'    \item{Cu:} {Underage cost, i.e., cost per unit for lost sales (c.f. penalty for unsatisfied demand)}
#' }}
#' \item{Alternative II:} {Specify
#'     \itemize{
#'     \item{p:} {Selling price per unit}
#'     \item{c:} {Purchasing price per unit (cost)}
#'     \item{s:} {Salvage price per unit (The relationship, p > c > s, is assumed)}
#'     }
#'     and the function computes the corresponding over- and underage costs from these.
#' }
#' }
#'
#' @param{mu} {demand mean}
#' @param{sigma} {standard deviation of the demand}
#' @param{Co} {overage cost}
#' @param{Cu} {underage cost}
#' @param{p}  {selling price}
#' @param{c}  {purcase price}
#' @param{s}  {salvage price}
#'
#' @return A list containing the following:
#' \itemize{
#'   \item{S}  {the optimal order quantity}
#'   \item{SS}  {the safety stock}
#'   \item{z}  {the safety factor}
#'   \item{S2}  {the service level of type II}
#'   \item{CQ}  {the critical quotient}
#'   \item{CV}  {the coefficient of variation of the demand}
#' }
#' In the case of input of alternative II (see details), the following additional elements are alo returned:
#' \itemize{
#'   \item{EC} {the optimal expected cost}
#'   \item{EP} {the optimal expected profit}
#' }
#'
#' @references
#'  - Axsäter, S. (2015) Inventory Control, 3d Edition, Springer
#'  - Snyder, L. V. and Shen, Z.-J., M. (2011) Fundamentals of Supply Chain Theory, Wiley
#'
#' @examples
#' ## ** Example 1 **
#' ## Data based on Example 5.8 in Axsäter's book
#' ##
#' Newsvendor(300,60,25,50)
#' ##
#' Newsvendor(300,60,Co=25,Cu=50)
#' ##
#' Newsvendor(300,60,p=75,c=25)
#' ##
#' Newsvendor(300,60,p=75,c=25, s=0)
#'
#' ## ** Example 2 **
#' ## Data based on Example 4.4 in Snyder-Shen's book
#' ##
#' Newsvendor(mu=50, sigma=8, Co=0.18, Cu=0.70)
#'
#'## ** Example 3 **
#'## Assume that the newsvendor buys each newspaper for $1 and sells it
#'## for $5. Assume further that the demand is normally distributed with
#'## μ=150 and σ=25. What is the optimal order quantity?
#'##
#'Newsvendor(mu=150, sigma=25, c=1, p=5)
#'##
#'## Suppose that unsold newspapers are salvaged for $0.4 per unit.
#'## Would the optimal order quantity change?
#'##
#'Newsvendor(mu=150, sigma=25, c=1, p=5, s=0.4)
#'
#'
#'@export
Newsvendor <- function(mu,sigma, Co=0, Cu=0, p=0,c=0,s=0){
  if( Co > 0 & Cu > 0 ) {
    CQ <- Cu / ( Cu + Co)
    z <- stats::qnorm(CQ)
    S <- mu + z * sigma
    CV <- sigma/mu
    S2 <- 1 - CV * ( stats::dnorm(z) - ( 1 - stats::pnorm(z) ) * z )
    SS <- z * sigma
    #options(digits=3)
    result <- list(S=S, SS=SS, z=z, S2=S2, CQ=CQ, CV=CV)
    return(result)
  }else if(p > c & c > s){
    Co <- c -s
    Cu <- p - c
    CQ <- Cu / ( Cu + Co)
    c1 <- p - s
    c2 <- p - c
    #c3=c2/c1
    CV <- sigma/mu
    z <- stats::qnorm(CQ)
    #z<- stats::qnorm(p=c3, mean = 0, sd = 1)
    S <- mu + z * sigma
    S2 <- 1 - CV * ( stats::dnorm(z)-( 1 - stats::pnorm(z) ) * z )
    SS <- z * sigma
    EC <- c1 * sigma * stats::dnorm(z)
    EP <- c2 * mu-EC
    #options(digits=3)
    result <- list(S=S, SS=SS, z=z, S2=S2,  CQ=CQ, CV=CV, EC=EC, EP=EP)
    return(result)
  }else {
    return(NaN)
  }
}
