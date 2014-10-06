plot.binaryTree <- function(x, type = "value", precision = 2, ...)
{
	x.coord <- 0.5
	y.coord <- -0.5
	for(i in 1:x$max.s)
	{
		x.coord <- c(x.coord, seq(0,1, by=1/(2^i+1)))
		y.coord <- c(y.coord, rep(-i-0.5, 2^i))
	}
	x.coord <- x.coord[which(x.coord>0 & x.coord<1)]
	x0 <- rep(x.coord[-c((2^x$max+1):2^(x$max.s+1)-1)],each=2)
	x1 <- x.coord[-1]
	y0 <- rep(y.coord[-c((2^x$max+1):2^(x$max.s+1)-1)],each=2)
	y1 <- y.coord[-1]
	plot(0,0, ylim=c(-x$max.s-1,0), xlim=c(0,1), col=0, ylab="scales", xlab="", tcl=0, col.axis=0, ...)
	segments(x0,y0,x1,y1)
	if(type=="value")
	{
		points(x.coord,y.coord, pch=21, cex=3, bg="white", col=0)
		text(x.coord,y.coord, round(tree2vec(x),precision))
	}
	if(type=="color")
	{
		if(max(tree2vec(x))<1) cat("To plot colors on the nodes, the tree must contain integers"); return() 
		points(x.coord,y.coord, pch=21, cex=3, bg=gray.colors(max(tree2vec(x)))[tree2vec(x)], col=1)
	}
	if(type=="size")
	{
		biggest <- max(tree2vec(x))
		scale <- 4/biggest
		points(x.coord,y.coord, pch=21, cex=scale*tree2vec(x), bg="white")
	}
}
#----------------
summary.binaryTree <- function(object, ...)
{
	cat("Binary Tree with ", object$max.s, "scales \n")
	print(object$T)
}
#-----------------

