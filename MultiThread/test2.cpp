#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <myompsubs.h>
#include <mychebfuns.h>
#include <HamOp.h>
#include <ChebPoly.h>
#include <ChebPolyOp.h>
#include "arssym.h"
//#include "smatrixa.h"
#include "symsol.h"


// function to approximate
double testfun(double x){
	
	double y, eps = 1e-3, shift = 0.0;

	//y = 1.0/((x-shift)+eps);
	//y = exp(-(x-shift)*(x-shift)/eps);
	//y = exp(-(x-shift)/eps);
	//y = (1.0-x)*(1.0-x)*(1.0-x)*(1.0-x)*(1.0-x)*(1.0-x)*(1.0-x)*(1.0-x);
	//y = eps/(x-shift);
	//y = pow(1.0-x,100.0);
	y = x;

	return(y);
}

int main(int argc, char *argv[]){

	// compute variables
	int N, M, ii, flag;
	int Dims, Nx, Ny, Nz;
	int sd = sizeof(double);
	float time;
	double a, b;
	double xmin, xmax, ymin, ymax, zmin, zmax, specrad;
	double *vecin, *vecout, *coeffs;
	double begin, end;

	// set variables
	a = 0.0;
	b = 1.0;
	Dims = 3;
	Nx = pow(2,4);
	Ny = pow(2,4);
	Nz = pow(2,4);
	xmin = -1.0;
	xmax = 1.0;
	ymin = -1.0;
	ymax = 1.0;
	zmin = -1.0;
	zmax = 1.0;
	specrad = 10.0;
	HamOp H(Dims,Nx,Ny,Nz,xmin,xmax,ymin,ymax,zmin,zmax,&wr2,specrad);
	N = H.getNx()*H.getNy()*H.getNz();
	
	// print info
	printnumthreads();
	printf("N = %d\n",N);

	// allocate memory
	vecin = (double*)malloc(N*sd); if(vecin == NULL){printf("Out of memory\n"); exit(-1);}
	vecout = (double*)malloc(N*sd); if(vecout == NULL){printf("Out of memory\n"); exit(-1);}
	coeffs = (double*)malloc((pow(2,17)+1)*sd); if(coeffs == NULL){printf("Out of memory\n"); exit(-1);}

	// construct polynomial
	begin = omp_get_wtime();
	chebctor(&testfun, a, b, &M, coeffs, &flag);
	ChebPoly P(M,a,b,coeffs);
	end = omp_get_wtime();
	printf("time to construct poly = %f\n",end-begin);

	ChebPolyOp C(&P,&H);
	H.print();
	printf("M = %d\n",P.getM());

	// fill arrays
	begin = omp_get_wtime();
	init(N, 1.0, vecin);
	init(N, 0.0, vecout);
	end = omp_get_wtime();
	printf("time to init arrays = %f\n",end-begin);

	// call chebmult
	begin = omp_get_wtime();
	C.multChebPolyOp(vecin,vecout);
	end = omp_get_wtime();
	printf("time to mult = %f\n",end-begin);


	// Defining what we need: the four eigenvectors of A with smallest magnitude.
	// A.MultMv is the function that performs the product w <- A.v.
	ARSymStdEig<double, ChebPolyOp> dprob(N, 4L, &C, &ChebPolyOp::multChebPolyOp, "LM");

	// Finding eigenvalues and eigenvectors.
	dprob.FindEigenvectors();

	// Printing solution.
	Solution(C, dprob);

	H.setspecrad(H.getspecrad()*dprob.Eigenvalue(3));
	H.print();
	ChebPolyOp D(&P,&H);

	// Defining what we need: the four eigenvectors of A with smallest magnitude.
	// A.MultMv is the function that performs the product w <- A.v.
	ARSymStdEig<double, ChebPolyOp> dprob2(N, 4L, &D, &ChebPolyOp::multChebPolyOp, "LM");

	// Finding eigenvalues and eigenvectors.
	dprob2.FindEigenvectors();

	// Printing solution.
	Solution(D, dprob2);


	// print arrays
/*	printf("\nvecin, vecout\n");
	for(ii=0;ii<N;ii++){
		printf("%+f, %+f\n",vecin[ii],vecout[ii]);
	}
*/
	// free memory
	free(vecin);
	free(vecout);

}

