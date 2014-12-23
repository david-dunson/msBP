\name{msBP.postCluster}
\alias{msBP.postCluster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Posterior cluster allocation
}
\description{
Perform the posterior multiscale cluster allocation conditionally on a tree of weights. See Algorithm 2 in Canale and Dunson (2014).
}
\usage{
msBP.postCluster(y, weights)
}
\arguments{
  \item{y}{
the sample of individials to be allocated to binary tree structure}
  \item{weights}{
the binary tree of weights (summing to one). An object of the class msBPTree
}
}
\details{
conditionally on the weights contained in \code{weights}, each subject in \code{y} is allocated to a multiscale cluster using Algorithm 2 of Canale and Dunson (2014). It relies on a multiscale modification of the slice sampler of Kalli et al. (2011). 
}
\value{
a matrix with \code{length(y)} row and two columns, denoting the scale and node within the scale, respectively.
}
\references{
Canale, A. and Dunson, D. B. (2014), "Multiscale Bernstein polynomials for densities", arXiv:1410.0827

Kalli, M., Griffin, J., and Walker, S. (2011), "Slice sampling mixture models," Statistics and Computing, 21, 93-105.
}
\seealso{
\code{\link{msBP.Gibbs}}
}
\examples{
set.seed(1)
y <- rbeta(30, 5, 1)
weights <-structure(list( 
	T = list(0, c(0,0.10), c(0.0,0,0.3,0.6)), max.s=2), 
	class  = 'binaryTree')
sh <- msBP.postCluster(y, weights)
clus.size <- msBP.nrvTrees(sh)$n
plot(clus.size)
}
\keyword{multiscale clustering}

