% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PBtRED1.R
\name{PBtRED1}
\alias{PBtRED1}
\title{Partially Balanced Ternary Residual Effect Designs for Prime Number of Treatments}
\usage{
PBtRED1(v)
}
\arguments{
\item{v}{Prime Number of Treatments, v ( >= 5)}
}
\value{
It returns a new class of PBTREDs along with its parameters, Information Matrix (C), Average Variance Factor (AVF), and Canonical Efficiency Factor (CEF) for both treatment and residual effects.
}
\description{
This function generates a class of partially balanced ternary residual effect designs (PBTREDs) for a prime number of treatments (v), where (v >= 5) with p periods, take the value (v+3)/2, and n sequences, take the value v(v-1)/2.
}
\examples{
library(TREDesigns)
PBtRED1(v = 5)
}
