msBP.pdf <-
function(weights, n.points)
{
if(n.points<2) stop("n.points must be an integer greater than 2")
grid <- seq(0.001, 0.999, length=n.points)
w <- tree2vec(weights)
res <- .C("dmsBP_C", as.double(w), as.double(grid), as.integer(n.points), as.integer(weights$max.s), 
	ans = as.double(rep(0, n.points)))
list(x=grid,dens=res$ans)
}

msBP.cdf <-
function(weights, n.points,log=FALSE)
{
if(n.points<2) stop("n.points must be an integer greater than 2")
grid <- seq(0.001, 0.999, length=n.points)
w <- tree2vec(weights)
res <- .C("pmsBP_C", as.double(w), as.double(grid), as.integer(n.points), as.integer(weights$max.s), 
	ans = as.double(rep(0, n.points)), as.integer(log))
list(x=grid,prob=res$ans)
}
