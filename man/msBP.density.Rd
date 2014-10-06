\name{msBP.pdf}
\alias{msBP.pdf}
\alias{msBP.cdf}
\title{
Random probability density function and random cumulative distribution function  from a msBP
}
\description{
Compute the density and the cumulative distribution functions of a random density drawn from an msBP(a,b) process
}
\usage{
msBP.pdf(weights, n.points)
msBP.cdf(weights, n.points, log)
}
\arguments{
  \item{weights}{
An object of the class binaryTree containing probability weights
}
  \item{n.points}{
Length of the grid over (0,1) in which calculate the value of the random density
}
  \item{log}{
Logical. TRUE for computing the log-cdf
}
}
\value{
Vector of size \code{n.points}
}
\references{
Canale, A. and Dunson, D., B., (2014), Multiscale Bernstein polynomials for densities, arXiv:1410.0827
}
\seealso{
\code{\link{msBP.rtree}}
}
\examples{
prob <-structure(list( T = list(0.15, c(0.05,0.05), c(0.05,0.2,0.1,0.1), 
		c(0,0,0.3,0,0,0,0,0) ), max.s=3), class  = "binaryTree")
density <- msBP.pdf(prob, 100)
probability <- msBP.cdf(prob, 100)
par(mfrow=c(1,2))
plot(density$dens~density$x, ty='l', main = "pdf")
plot(probability$prob~density$x, ty='l', main = "cdf")
}
