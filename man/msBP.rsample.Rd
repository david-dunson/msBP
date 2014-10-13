\name{msBP.rsample}
\alias{msBP.rsample}

\title{
Random numbers from a random msBP densty 
}
\description{
Random numbers generation from a random density drawn from a msBP process. 
The function implements Algorithm 1 of Canale and Dunson (2014).
}
\usage{
msBP.rsample(n, msBPtree)
}
\arguments{
  \item{n}{
Size of the sample to be generated
}
  \item{msBPtree}{
An object of the class msBPtree
}
}
\value{
A vector containing the random sample
}
\references{
Canale, A. and Dunson, D. B. (2014), "Multiscale Bernstein polynomials for densities", arXiv:1410.0827
}
\seealso{
\code{\link{msBP.rtree}}
}
\examples{
rand.tree <- msBP.rtree(50,2, 4)
rand.samp <- msBP.rsample(50, rand.tree)
hist(rand.samp, prob=TRUE)
prob <- msBP.compute.prob(rand.tree)
density <- msBP.pdf(prob, 100)
points(density$dens~density$x, ty='l', col=4)
}
