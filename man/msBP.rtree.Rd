\name{msBP.rtree}
\alias{msBP.rtree}
\title{
Random msBP tree
}
\description{
Draw a random tree from the msBP process
}
\usage{
msBP.rtree(a, b, max.s = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{
Scalar parameter
}
  \item{b}{
Scalar parameter
}
  \item{max.s}{
Maximum depth of the random trees
}
}
\value{
An object of the class \code{msBPTree}
}
\references{
Canale, A. and Dunson, D. B. (2016), "Multiscale Bernstein polynomials for densities", Statistica Sinica, 26
}
\seealso{
\code{\link{msBP.rsample}}, \code{\link{msBP.compute.prob}} 
}
\examples{
msBP.rtree(2, 2, 4)
}
