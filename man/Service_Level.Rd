% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Service_Level.R
\name{Service_Level}
\alias{Service_Level}
\alias{S1}
\alias{S1b2}
\alias{S2}
\alias{S2b1}
\alias{b1S2}
\alias{b2S1}
\title{Service level functions for single echelon stochastic demand inventory model}
\usage{
S1(sigma, mu, R)

S1b2(h, b2, mu, Q)

S2(sigma, mu, Q, R)

S2b1(h, b1)

b1S2(h, S2)

b2S1(h, mu, Q, S1)
}
\arguments{
\item{{sigma}}{{standard deviation of demand}}

\item{{mu}}{{mean value of demand}}

\item{{R}}{{reorder point}}

\item{{Q}}{{order quantity}}

\item{{h:}}{{holding cost}}

\item{{b1}}{{shortage cost per unit and time unit}}

\item{{b2}}{{shortage cost per unit}}

\item{{S1}}{{service level of type I}}

\item{{S2}}{{service level of type II}}
}
\value{
A single numeric value or NULL, if wrong input
}
\description{
Functions to compute service level of type I (S1) and type II (S2) as
well as shortage
costs (b1 and b2), for given service level,
in the case of normal distribution of the demand
}
\details{
Function S1(): Service level S1
probability of no stockout per order cycle

Function S1b2(): Service level S1, i.e., probability of no stockout
per order cycle in terms of shortage cost per unit, b2
(Approximation)

Function S2(): S2 service level, i.e., "fill rate" fraction, also
S3 service level, i.e., "ready rate" fraction

Function S2b1(): Fill-rate, ready-rate Service Level S2, S3
w.r.t. dhortage cost per unit and time unit, b1

Function b1S2(): Shortage cost per unit per time, b1,
in terms of fill-rate, ready-rate service level S2, S3

Function b2S1(): Shortage cost per unit, b2, in terms of
Service level S1, i.e., probability of no stockout
per order cycle (Approximation)
}
\examples{
S1(sigma, mu, R)
S1b2(h,b2,mu,Q)
S2(sigma,mu,Q,R)
S2b1(h,b1)
b1S2(h,S2)
b2S1(h,mu,Q,S1)
}
\references{
\itemize{
\item Axsäter, S. (2015) Inventory Control, 3d Edition, Springer
}
}
\author{
Athanasios Migdalas, \email{athmig59@gmail.com}
}
