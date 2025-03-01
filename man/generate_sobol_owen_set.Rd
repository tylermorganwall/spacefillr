% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_values.R
\name{generate_sobol_owen_set}
\alias{generate_sobol_owen_set}
\title{Generate Owen-scrambled Sobol Set}
\usage{
generate_sobol_owen_set(n, dim, seed = 0)
}
\arguments{
\item{n}{The number of values (per dimension) to extract.}

\item{dim}{The number of dimensions of the sequence. This has a maximum value of 21201.}

\item{seed}{Default `0`. The random seed.}
}
\value{
An `n` x `dim` matrix with all the calculated values from the set.
}
\description{
Generate a set of values from an Owen-scrambled Sobol set.
}
\examples{
#Generate a 2D sample:
points2d = generate_sobol_owen_set(n=1000, dim = 2)
plot(points2d, xlim=c(0,1),ylim=c(0,1))

#Generate a longer sequence of values from that set
points2d = generate_sobol_owen_set(n=1500, dim = 2)
plot(points2d, xlim=c(0,1),ylim=c(0,1))

#'#Integrate the value of pi by counting the number of randomly generated points that fall
#within the unit circle.
pointset = matrix(generate_sobol_owen_set(10000,dim=2),ncol=2)

pi_estimate = 4*sum(pointset[,1] * pointset[,1] + pointset[,2] * pointset[,2] < 1)/10000
pi_estimate
}
