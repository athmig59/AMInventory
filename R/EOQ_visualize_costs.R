#' @title Visualization of the costs involved in the EOQ model
#'
#' @author Athanasios Migdalas, \email{athmig59@gmail.com}
#'
#' @description
#' A plot of the ordering, holding and total costs is produced for a range
#' of order quantities provided by the user. The optimal order quantity (EOQ)
#' is computed and plotted too.
#'
#' @references
#'  - Axsäter, S. (2015) Inventory Control, 3d Edition, Springer
#'
#' @seealso   EOQsensitivity, DEOQ, IDEOQ, EOQ
#'
#' @param {d} {Yearly demand rate, default value 1500}
#' @param {A} {Fixed ordering cost, default value 10}
#' @param {h} {Yearly holding cost per unit, default value 0.65}
#' @param {Q} {a vector containing range of demands around EOQ. Provide an empty vector c() to let
#'             the automatic generation of a range. Default range is 30:300 by a step 10}
#' @returns A list consisting of two elements:
#' - {Plot} {contains the plot of the costs}
#' - {Table} {is data frame of four columns:\itemize{
#'            \item{Q} {the range of order quantities}
#'            \item{HC} {the corresponding range of holding costs}
#'            \item{OC} {the corresonding range of order costs}
#'            \item{TC} {the corresponding range of total costs}}}
#'  - or NULL, if supplied with inadequate holding cost (<=0)
#'
#' @examples
#' \dontrun{** Example 1 **
#' Run without arguments for the default plot and default output, like this:}
#'
#' EOQ_visualize_costs()
#'
#' \dontrun{** Example 2 **
#' Store the result in a list and view it using print(), like this:}
#'
#' res <- EOQ_visualize_costs()
#' print( res$Plot )
#' print( res$Table )
#'
#' \dontrun{** Example 3 **
#' You can supply your own data, like this:}
#'
#' res <- EOQ_visualize_costs(A=100, d=300, h=20, Q=seq(20,100,by=5))
#'
#' @importFrom ggplot2 ggplot geom_line geom_label aes xlab ylab ggtitle labs geom_vline
#'
#'@export
EOQ_visualize_costs <- function(A=10, d=1500, h=0.65, Q=seq(30,300,by=10)) {

  if(h <= 0){return(NULL)}

  EOQ <- round( sqrt(2*A*d/h) )

  l <- length(Q)

  if( l == 0){
    Q <- seq(10, EOQ, by=EOQ/10)
    Q <- append(Q, seq(EOQ, 2*EOQ, by= EOQ/10))
    Q <- round(Q)
  }

  HC <- function(Q,h){
    return(Q*h/2)
  }

  OC <- function(Q,d,A){
    return(d*A/Q)
  }


  hc <- c()
  oc <- c()
  tc <- c()

  for ( q in Q) {
    hc <- append(hc, HC(q,h))
    oc <- append(oc, OC(q,d,A))
    tc <- append(tc, HC(q,h)+OC(q,d,A))
  }

  df <- data.frame(Q=Q, HC=hc, OC=oc, TC=tc)

  l <- length(Q)
  maxhc <- max(hc)
  minoc <- min(oc)
  mintc <- min(tc)

  #library(ggplot2)

  p <- ggplot2::ggplot(df, aes(x=Q))+
    geom_line(aes(y=HC), color='red')+
    geom_label(x=Q[l],y=maxhc, label='HC')+
    geom_line(aes(y=OC), color='blue')+
    geom_label(x=Q[l],y=minoc, label='OC')+
    geom_line(aes(y=TC), color='darkred')+
    geom_label(x=Q[l],y=mintc*1.2, label='TC')+
    geom_label(x=2+EOQ, y=0, label='EOQ')+
    labs(x='Order Quantity', y='Costs')+
    labs(title='EOQ cost functions') +
    ggtitle("EOQ costs: Order cost (OC) vs Hoding cost (HC) vs Total cost (TC)")

  p <- p + geom_vline(xintercept=EOQ)

  return( list(Plot=p, Table=df) )
}

