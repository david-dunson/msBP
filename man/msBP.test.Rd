\name{msBP.test}
\alias{msBP.test}
\title{
Multiscale testing of group differences}
\description{
Performs multiscale hypotesis testing of difference in the distribution of two groups using msBP prior.
}
\usage{
msBP.test(y, a, b, group, priorH0 = 0.5, g0 = "uniform", 
		mcmc, maxScale = 5, plot.it = FALSE, ...)
}
\arguments{
  \item{y}{The pooled sample of observations}
  \item{a,b}{Parameters of the msBP prior}
  \item{group}{Vector of size \code{length(y)} with 0 and 1 denoting the group membership.}
  \item{priorH0}{Prior gues for the probability of H0}
 \item{g0}{prior guess for the density of \code{y}. Currently only "uniform" is supported.}
  \item{mcmc}{a list giving the MCMC parameters. It must include the
          following integers: \code{nb} giving the number of burn-in iterations, \code{nrep} giving
          the total number of iterations (including \code{nb})., and \code{ndisplay} giving
          the multiple of iterations to be displayed on screen while the \code{C++} routine is running (a message will be printed every \code{ndisplay} iterations).} 
  \item{maxScale}{maximum scale of the binary trees.}
  \item{plot.it}{logical. If TRUE a plot of the posterior mean probability of H0 is produced}
  \item{\dots}{additional arguments to be passed.}
}
\value{
a list containing
\item{Ps}{the posterior probabilities of H0 for each scale}
\item{Ps}{the posterior probabilities of the global H0, obtained cumulating \code{Ps}.}
}
\references{
Canale, A. and Dunson, D. B. (2014), "Multiscale Bernstein polynomials for densities", arXiv:1410.0827
}
\examples{
set.seed(1)
y <- runif(100)
g <- c(rep(0,50), rep(2,50))
mcmc <- list(nrep = 5000, nb = 1000, ndisplay = 500)
\dontrun{
test.res <- msBP.test(y, 5, 1, g, mcmc=mcmc, plot.it = TRUE)
}
}
\keyword{ multiscale testing }
