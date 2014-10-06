\name{msBP.tree}
\alias{msBP.tree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creating an \code{msBPTree}
}
\description{
Create an object of the class \code{msBPTree}
}
\usage{
msBP.tree(max.s = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{max.s}{
Maximum depth of the binary tree
}
}
\details{
An object of the class  \code{msbpTree} is a list of 5 elements that represent a draw from a msBP(a,b) prior as introduced by Canale and Dunson (2014). The first two elements are the trees of the stopping and descending-to-the-right probabilities, respectively. Both are object of the class \code{binaryTree}. The third and fourth argument are the hyperparameters of the msBP prior, namely \code{a} and \code{b}. The last value is an integer with the maximum depth of both the trees.}
\value{
An object of the class \code{msBPTree} with zero at all nodes and \code{a=b=NULL}.
}
\references{
Canale, A. and Dunson, D., B., (2014), Multiscale Bernstein polynomials for densities, arXiv:1410.0827
}
\author{
Antonio Canale <antonio.canale@unito.it>
}

\seealso{
\code{msBP.rtree}
}
\examples{
tree <- msBP.tree(2)
}
