#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_cblas.h>
#include <myompsubs.h>
#include <myfortompsubs.h>
#include <mychebfuns.h>
#include <HamOp.h>
#include <ChebPoly.h>
#include <ChebPolyOp.h>

// function to approximate
double testfun(double x){
	
	double y, eps = 1e-4, shift = 0.0;

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
	double a, b, scl, nrm, dt, res;
	double xmin, xmax, ymin, ymax, zmin, zmax, specrad;
	double *vecin, *vecout, *coeffs;
	double begin, end;

	// set variables
	a = 0.0;
	b = 1.0;
	Dims = 1;
	Nx = pow(2,7);
	Ny = pow(2,0);
	Nz = pow(2,0);
	xmin = 0.0;
	xmax = 1.0;
	ymin = 0.0;
	ymax = 1.0;
	zmin = 0.0;
	zmax = 1.0;
	specrad = 1.0;
	HamOp H(Dims,Nx,Ny,Nz,xmin,xmax,ymin,ymax,zmin,zmax,&zero,specrad);
	N = H.getNx()*H.getNy()*H.getNz();

	// print info
	printnumthreads();
	printf("N = %d\n",N);

	// compute specrad
	begin = omp_get_wtime();
	//H.computespecrad();
	end = omp_get_wtime();
	H.print();
	printf("time to compute specrad = %f\n",end-begin);
	H.setspecrad(1.0);

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

	//ChebPolyOp C(&P,&H);
	ChebPolyOp C(P,H);
	C.print();
	
	// fill arrays
	begin = omp_get_wtime();
	init(N, 1.0/sqrt((double)N), vecin);
	init(N, 0.0, vecout);
	end = omp_get_wtime();
	printf("time to init arrays = %f\n",end-begin);

	begin = omp_get_wtime();
	for(ii=0;ii<1000;ii++){
		printf("\niter = %d\n",ii+1);
		// call chebmult
		//C.multChebPolyOp(vecin,vecout);
		H.multHamOp(vecin,vecout);

		dt = dot(N,vecout,vecin);
		printf("dot = %e\n",dt);

		scaleadd(N,-dt,vecin,vecout);
		res = norm(N,vecin);
		printf("res = %e\n",fabs(res/dt));

		nrm = norm(N,vecout);
		printf("norm = %e\n",nrm);
		scale(N,1.0/nrm,vecout);
		copy(N,vecout,vecin);
	}
	end = omp_get_wtime();
	printf("time to mult = %f\n",end-begin);


/*	// print arrays
	printf("\nvecin, vecout\n");
	for(ii=0;ii<N;ii++){
		printf("%+e, %+e\n",vecin[ii],vecout[ii]);
	}
*/
	// free memory
	free(vecin);
	free(vecout);

}

