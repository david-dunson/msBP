msBP.test <- function(y, a, b, group, priorH0 = 0.5, g0 = "uniform", mcmc, maxScale=5, plot.it=FALSE, ...)
{
	priorH1 <- 1 - priorH0
	w <- w1 <- w0 <- msBP.compute.prob(msBP.rtree(a,b,maxScale), root=FALSE)
	Ps <- matrix(0.5, maxScale, mcmc$nrep)
	for(ite in 2:mcmc$nrep)
	{
		if(ite/100 == round(ite/mcmc$ndisplay)) cat(ite)
		#cat(ite)
		Pm <- rep(c(1,Ps[,ite-1]), c(2^(0:(maxScale))))
		W1 <- vec2tree(Pm*tree2vec(w) + (1-Pm)*tree2vec(w1))
		W0 <- vec2tree(Pm*tree2vec(w) + (1-Pm)*tree2vec(w0))
		sh0 <- msBP.postCluster(y[group==0], W0)
		sh1 <- msBP.postCluster(y[group==1], W1)
		nrv0 <- msBP.nrvTrees(sh0, maxS=maxScale)
		nrv1 <- msBP.nrvTrees(sh1, maxS=maxScale)
		n0 <- tree2vec(nrv0$n)
		r0 <- tree2vec(nrv0$r)
		v0 <- tree2vec(nrv0$v)
		n1 <- tree2vec(nrv1$n)
		r1 <- tree2vec(nrv1$r)
		v1 <- tree2vec(nrv1$v)
		n <- n0+n1
		r <- r0+r1
		v <- v0+v1
		S <- rbeta(length(n), 1+n, a+v-n)
		S0 <- rbeta(length(n0), 1+n0, a+v0-n0)
		S1 <- rbeta(length(n1), 1+n1, a+v1-n1)
		R <- rbeta(length(n), b+r, b+v-n-r)
		R0 <- rbeta(length(n), b+r, b+v0-n0-r0)
		R1 <- rbeta(length(n), b+r1, b+v1-n1-r1)
		w <- msBP.compute.prob(structure(list(S = vec2tree(S), R = vec2tree(R)), class  = "msbpTree"), root=FALSE)
		w0 <- msBP.compute.prob(structure(list(S = vec2tree(S0), R = vec2tree(R0)), class  = "msbpTree"), root=FALSE)
		w1 <- msBP.compute.prob(structure(list(S = vec2tree(S1), R = vec2tree(R1)), class  = "msbpTree"), root=FALSE)
		testingprobs <- msBP.nesting(n,r,v,n0,r0,v0,n1,r1,v1,priorH0,a,b,maxScale)
		Ps[,ite] <- testingprobs[1,]
	}
	Ps
	Pp <- apply(Ps,1,mean)
	if(plot.it) plot(cumprod(Pp)~c(1:maxScale), ty='b', xlab="scale", ylim=c(0,1), ylab="P(H_0 | - )", ...)	
	out <- list(Ps=Ps, Pp=Pp)
}
