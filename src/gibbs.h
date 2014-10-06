/*-------------------------------------------------------------------------------
 Multiscale Bernstein polynomial [msBP]
 gibbs.hpp - Header fot Gibbs sampling code 
 Version 2.2 of May 2014
 2013 - Antonio Canale (antonio.canale@unito.it)
-------------------------------------------------------------------------------*/
void msBPgibbs(double *y, double *par, int *sClus, int *hClus, 
	double *Rstart, double *Sstart, double * wstart, int *hyperpriors, 
	int *nrep, int *nb, int *aux, int *printing, double *grid, int *ngrid, double *griddy, int *griddy_length,
	double *postDens, double *postScale, double *postS, double *postR, double *postPi, 
	double *postA, double *postB, int *posts, int *posth, int *type)
