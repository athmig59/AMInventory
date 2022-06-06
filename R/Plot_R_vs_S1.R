#' @title Plots Reorder Point (R) as function of Service Level Type 1 (S1)
#'
#' @author Athanasios Migdalas, \email{athmig59@gmail.com}
#'
#' @description
#' Plots Reorder Point (R), Safety Stock (SS) and Safety Factor (z) as functions
#' of Service Level Type 1 (S1). The purpose is to demonstrate the 
#' nonlinearity of the relations, and particularly the accelation
#' in the increase that occurs as the 100% limit of the service
#' level is approached.
#'
#' @references
#'  - Axsäter, S. (2015) Inventory Control, 3d Edition, Springer
#'
#' @seealso Reorder_Point_S1, Reorder_Point_S2
#' 
#' @return A list with two elements:
#'   \itemize{
#'   \item{\code{Plot}} {contains the plot. If stored, can be desplayed with \code{print()}}
#'   \item{\code{Table}} {contains in four columns, \code{S1}, \code{R}, \code{SS} and \code{z}, the corresponding values for the specified range of service level 1.}
#'   }
#' 
#' @examples 
#' ## ** Example 1 **
#' ## For demand with mean 10 and standard deviation sqrt(10), leaving the lower limit
#' ## of the service level to the default 50%:
#' ##
#' Plot_R_vs_S1(sigma=sqrt(10), mu=10)
#' 
#' ## ** Example 2 **
#' ## As above but specifying a lower limit of 30%
#' ##
#' Plot_R_vs_S1(S1low=0.3, sigma=sqrt(10), mu=10)
#' 
#' @export
Plot_R_vs_S1 <- function(S1low=0.5, sigma, mu){
  
  # Service level type 1
  
  S1 <- seq(S1low,1.0, by=0.05)
  l<- length(S1)
  S1[l] <- 0.99 
  
  # Corresponding safety stock, reorder point and safety factor
  
  SS <- rep(NaN,l)
  R  <- rep(NaN,l)
  z  <- rep(NaN,l)
  
  for(i in 1:l){
    res <- Reorder_Point_S1(S1[i], sigma, mu)
    R[i] <- res$R
    SS[i] <- res$SS
    z[i] <- res$z
  }
  
  df <- data.frame(S1=S1, R=R, SS=SS, z=z)
  
  p <- ggplot2::ggplot(df, aes(x=S1))+
    geom_smooth(formula="y~x", aes(y=R, color="R"))+
    geom_smooth(formula="y~x", aes(y=SS, color="SS"))+
    geom_smooth(formula="y~x", aes(y=z, color="z"))+
    xlab("Service Level Type 1 (S1)") +
    ylab(" R / SS / z")
  
  return(list(Plot=p, Table=df))
}