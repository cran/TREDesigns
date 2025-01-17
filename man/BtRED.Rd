% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BtRED.R
\name{BtRED}
\alias{BtRED}
\title{Balanced Ternary Residual Effect Designs for Prime Number of Treatments}
\usage{
BtRED(v)
}
\arguments{
\item{v}{Prime Number of Treatments, v ( >= 5)}
}
\value{
It returns a new class of BTREDs along with its parameters, Information Matrix (C), Average Variance Factor (AVF), and Canonical Efficiency Factor (CEF) for both treatment and residual effects.
}
\description{
This function generates a class of balanced ternary residual effect designs (BTREDs) for a prime number of treatments (v), where (v >= 5) with p periods, determined as (v+3)/2, and n sequences, take the value v(v-1).
}
\examples{
library(TREDesigns)
BtRED(v = 7)
}
