/*-------------------------------------------------------------------------------
 Multiscale Bernstein polynomials [msBP]
 gibbs.cpp - Gibbs sampling code for msBP density estimation [Algorithm 3]
 Version 0.3 of September 2014
 2014 - Antonio Canale (antonio.canale@unito.it)
-------------------------------------------------------------------------------*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<R.h>
#include<Rmath.h>
#include "bintree.h"
#include "msBP.h"
#include "auxGibbs.h"
//------------------------------------------------------------------------------
extern "C"{
void msBPgibbs(double *y, double *par, int *sClus, int *hClus, 
	double *Rstart, double *Sstart, double * wstart, int *hyperpriors, 
	int *nrep, int *nb, int *aux, int *printing, double *grid, int *ngrid, double *griddy, int *griddy_length,
	double *postDens, double *postScale, double *postS, double *postR, double *postPi, 
	double *postA, double *postB, int *posts, int *posth, int *type)
{
	int i, j, s, h;
	int N = aux[0];  // 1st element of the auxiliary parameters is the dimension n
	int MAXS = aux[1];  // 2nd element of the auxiliary parameters is the maximum depth allowed
	int MAXVEC = aux[2];  // 3rd element of the auxiliary parameters is length of vector induced by MAXS
	int maxSstart = aux[3];  // 4th element of the auxiliary parameters is scale of starting trees
	int maxS = maxSstart;
	int maxH = 1;
	int flag = 0;
	double a, b;
	a = par[0];
	b = par[1];
	double delta, gamma, lambda, beta;
	beta = par[2];
	gamma = par[3];
	delta = par[4];
	lambda = par[5];
	double nsh, rsh, vsh;
	int every = printing[0];
	int printclustering = printing[1];
	double intoR = 0.5;
	double intoS = 0;
	double log1_S = 0;
	int *ZERO;
	ZERO = ( int* ) R_alloc(1, sizeof(int));
	ZERO[0] = 0;

	//we start from trees of depth 4 (arbitrary choice)
	struct bintree *S = newtree(1); //rStree(a, 4);
	struct bintree *R = newtree(0.5); //rRtree(b, 4);
	struct bintree *w = newtree(1); //computeprob(S, R, a, b, 4, 0);
	array2tree(Sstart, maxSstart, S);
	array2tree(Rstart, maxSstart, R);
	array2tree(wstart, maxSstart, w);

	struct bintree *n = newtree(0);
	struct bintree *r = newtree(0);
	struct bintree *v = newtree(0);
	auxiliaryTrees(sClus, hClus, N, n, r, v);

	GetRNGstate();
	for(i=1; i<nrep[0]; i++)
	{
		if(((i-1) % every) == 0) flag = 1;
		else flag = 0;
		if(flag)  Rprintf("Iteration %i over %i \n",i-1,nrep[0]);
  		R_CheckUserInterrupt();
	/*	printTree(S, maxS);
		printTree(R, maxS);*/
		postCluster(sClus, hClus, y, w, maxS+1, N, printclustering); 
		clearTree(n);
		clearTree(r);
		clearTree(v);
		allTrees(sClus, hClus, maxS, N, n, r, v);
		maxS = 0;
		for(j=0; j<N; j++) maxS = (int) fmax(maxS, sClus[j]);
		maxS = fmin(maxS, MAXS);

		clearTree(S);
		clearTree(R);
		clearTree(w);
		log1_S = 0;
		for(s=0; s<=maxS; s++)
		{
			maxH = (int) pow(2.0, s);		
			for(h=1; h<=maxH; h++)
			{
				vsh = (double) extractNode(v, s, h, 0);
				if(vsh!=0)
				{
				nsh = (double) extractNode(n, s, h, 0);
				rsh = (double) extractNode(r, s, h, 0);
				}
				else
				{
				nsh = 0; rsh = 0;	
				}
				intoS = rbeta(1 + nsh, postA[i-1] + vsh - nsh);
				intoR = rbeta(postB[i-1] + rsh, postB[i-1] + vsh - nsh - rsh);
				if(s==0)
				{
					writeNode(S, 0, s, h);
				}
				else 
				{
					writeNode(S, intoS, s, h);
					log1_S += log(1 - intoS);
				}
				writeNode(R, intoR, s, h);

			}
		}
		if(hyperpriors[0])
		{
			postA[i] = rgamma(beta+pow(2.0,maxS+1)-1, 1/(gamma - log1_S));
			postB[i] = griddy_B(delta, lambda, R, v, n, r, maxS, griddy, griddy_length[0]);
			if(isnan(postA[i]) || (postA[i]<=0)) {
			//Rprintf("NA! or 0 at ite %i\n", i+1);
			postA[i] = rgamma(beta, 1/gamma);
			}
			if(isnan(postB[i]) || (postB[i]<=0)) {
			//Rprintf("NA! or 0 at ite %i\n", i+1);
			postB[i] = rgamma(delta, 1/lambda);
			}
			//Rprintf("Sample a = %f, b = %f (ite %i)\n", postA[i-1], postB[i-1], i);
		}

		//compute probability tree
		w = computeprob(S, R, a, b, maxS, 0);
		//store posterior quantities
		if(i >= nb[0]) 
		{
			if(*type==0)
			{//total weight for each scale
			scaleProb(w, &(postScale[(MAXS+1)*(i-1)]));
			//vec version of S
			tree2array(S, &(postS[(i)*MAXVEC]), MAXS, 0);
			tree2array(R, &(postR[(i)*MAXVEC]), MAXS, 0);
			tree2array(computeprob(S, R, postA[i], postB[i], maxS, 1), &(postPi[(i)*MAXVEC]), MAXS, 0);
			memcpy(&(posts[(i)*N]), &(sClus[0]), N*sizeof(int));
			memcpy(&(posth[(i)*N]), &(hClus[0]), N*sizeof(int));
			dmsBP(computeprob(S, R, postA[i], postB[i], maxS, 1), &postDens[ngrid[0]*(i)], grid, ngrid);
			}
			if(*type==1)
			{
			// note that if type = 1 we are running Alg. 3 for group testing
			// hence we need to save the chains for v, n, and r descending trees 
			// but keep the names postS, postR, and postPi for simplicity
			tree2array(n, &(postS[(i)*MAXVEC]), MAXS, 0);
			tree2array(r, &(postR[(i)*MAXVEC]), MAXS, 0);
			tree2array(v, &(postPi[(i)*MAXVEC]), MAXS, 0);
			}
		}
	}
	PutRNGstate();
}

//------------------------------------------------------------------------------
}
