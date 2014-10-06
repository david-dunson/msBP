/*-------------------------------------------------------------------------------
 Multiscale Bernstein polynomials [msBP]
 wrapper.cpp - Set of wrapper calling C++ functions with bintree struc
 Version 0.1 of September 2013
 2013 - Antonio Canale (antonio.canale@unito.it)
-------------------------------------------------------------------------------*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<R.h>
#include<Rmath.h>
#include "bintree.h"
#include "msBP.h"
#include "auxGibbs.h"
extern "C" {
// compute the weigths given a tree 
void computeprob_C(double *S, double * R, int *maxS, double *a, double *b, double *ans, int *root)
{
	struct bintree *Stree = newtree(1);
	struct bintree *Rtree = newtree(1);
	struct bintree *anstree = newtree(1);
	if(root[0]==0) S[0] = 0;
	array2tree(S, maxS[0], Stree);
	array2tree(R, maxS[0], Rtree);
	anstree = computeprob(Stree, Rtree, a[0], b[0], maxS[0], 1);
	tree2array(anstree, ans, maxS[0], 0);
} 
// compute the density given a tree
void dmsBP_C(double *weights, double * grid, int * ngrid, int *maxS, double * out)
{
	struct bintree *w = newtree(1);
	array2tree(weights, maxS[0], w);
	dmsBP(w, out, grid, ngrid);
}
// compute the cdf given a tree
void pmsBP_C(double *weights, double * grid, int * ngrid, int *maxS, double * out, int *log_p)
{
	struct bintree *w = newtree(1);
	array2tree(weights, maxS[0], w);
	pmsBP(w, out, grid, ngrid, log_p);
}
//random trees
void randtree_C(double *a, double *b, int *maxS, double *ansS, double *ansR)
{
	struct bintree *S = newtree(1);
	struct bintree *R = newtree(1);
	R = rRtree (b[0], maxS[0]);
	S = rStree (a[0], maxS[0]);
	tree2array(S, ansS, maxS[0], 0);
	tree2array(R, ansR, maxS[0], 0);
}
// random sample from a msBP density
void rsample_msBP_C(int *N, double * Rvec, double *Svec, double *a, double *b, int * maxS, double *ans)
{
	int i = 0;	
	struct bintree *Stree = newtree(1);
	struct bintree *Rtree = newtree(1);
	array2tree(Svec, maxS[0], Stree);
	array2tree(Rvec, maxS[0], Rtree);
	for(i=0; i<N[0]; i++)
	{
		ans[i] = rsample_msBP(Rtree,Stree,a[0],b[0]);
	}
}
// marginal beta densities for a given y
void marginalBeta_C(double * out, double *y, int *maxS)
{
     marginalBeta(out, *y, *maxS);
}
//compute the n., r. and v. tree given two vectors (with N rows) of s and h indexes
void allTrees_C(int *s, int *h, int * maxS, int *N, double * nvec, double * rvec, double * vvec)
{
	struct bintree *n = newtree(1);
	struct bintree *r = newtree(1);
	struct bintree *v = newtree(1);
	array2tree(nvec, maxS[0], n);
	array2tree(rvec, maxS[0], r);
	array2tree(vvec, maxS[0], v);
	allTrees(s, h, *maxS, *N, n, r, v);
	tree2array(r, rvec, maxS[0], 0);
	tree2array(n, nvec, maxS[0], 0);
	tree2array(v, vvec, maxS[0], 0);
}
//compute the posterior cluster allocation with slice sampler [Algorithm 2]
void postCluster_C(int *s, int *h, double *y, double *pi, int *maxS, int *N, int *printscreen)
{
	struct bintree *pitree = newtree(1);
	array2tree(pi, maxS[0], pitree);
	postCluster(s, h, y, pitree, *maxS, *N, *printscreen);
} 
}
