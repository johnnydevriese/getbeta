#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <myfortompsubs.h>
#include <HamOp.h>
#include "arssym.h"

// zero potential
double zero(double x, double y, double z){

	return 0.0;
}

// one potential
double one(double x, double y, double z){

	return 1.0;
}

// weighted r^2 potential
double wr2(double x, double y, double z){

	return 2.0*(x*x + 4.0*y*y + 9.0*z*z);
}

// function to print a Hamop struct
void HamOp::print() const{
	printf("\nHamOp:\n");
	printf(" Dims = %d\n",Dims);
	printf(" Nx = %d\n",Nx);
	printf(" Ny = %d\n",Ny);
	printf(" Nz = %d\n",Nz);
	printf(" xmin = %+f\n",xmin);
	printf(" xmax = %+f\n",xmax);
	printf(" ymin = %+f\n",ymin);
	printf(" ymax = %+f\n",ymax);
	printf(" zmin = %+f\n",zmin);
	printf(" zmax = %+f\n",zmax);
	printf(" specrad = %+f\n\n",specrad);
}

void HamOp::computespecrad(){

	int N = Nx*Ny*Nz;

	ARSymStdEig<double, HamOp> dprob(N, 2, this, &HamOp::multHamOp, "LM", (int)fmin((double)40,(double)N));

	// Finding eigenvalues and eigenvectors.
	dprob.FindEigenvalues();

	setspecrad(fabs(specrad*dprob.Eigenvalue(1)));

}

void HamOp::multHamOp(double *vecin, double *vecout){

	// thread variables
	int nthds, tid;

	// compute variables
	int N, stride, start, stop, ii;
	int xpos, ypos, zpos;
	double hx, hy, hz;

	// Fork a team of threads giving them their own copies of variables
	#pragma omp parallel private(nthds, tid, N, stride, start, stop, ii, hx, hy, hz, xpos, ypos, zpos) shared(vecin, vecout)
	{
		// compute thread variables
		nthds = omp_get_num_threads();
		tid = omp_get_thread_num();

		// compute problem size 
		N = Nx*Ny*Nz;

		// compute stride
		stride = ceil((long double)N/nthds);

		// compute start and stop
		start = tid*stride;
		stop = (int)fminl((long double)(tid+1)*stride,(long double)N);

		// print info
		//printf("id, stride, start, stop = %d, %10d, %10d, %10d\n",tid,stride,start,stop);

		// compute directional differences
		hx = (double)(xmax - xmin)/(Nx + 1);
		hy = (double)(ymax - ymin)/(Ny + 1);
		hz = (double)(zmax - zmin)/(Nz + 1);

		// loop for multiplication
		for(ii=start;ii<stop;ii++){
			// compute xpos, ypos, zpos
			xpos = ii%Nx;
			ypos = ((ii-xpos)/Nx)%Ny;
			zpos = ((ii-xpos)/Nx-ypos)/Ny;

			// negative x laplacian
			if(xpos == 0){
				vecout[ii] = (2.0*vecin[ii] - vecin[ii+1])/hx/hx/specrad;
			}
			else if((xpos+1) == Nx){
				vecout[ii] = (2.0*vecin[ii] - vecin[ii-1])/hx/hx/specrad;
			}
			else if(ii < N){
				vecout[ii] = (2.0*vecin[ii] - vecin[ii+1] - vecin[ii-1])/hx/hx/specrad;
			}

			// negative y laplacian
			if(Dims == 2){
				if(ypos == 0){
					vecout[ii] = vecout[ii] + (2.0*vecin[ii] - vecin[ii+Nx])/hy/hy/specrad;
				}
				else if((ypos+1) == Ny){
					vecout[ii] = vecout[ii] + (2.0*vecin[ii] - vecin[ii-Nx])/hy/hy/specrad;
				}
				else if(ii < N){
					vecout[ii] = vecout[ii] + (2.0*vecin[ii] - vecin[ii+Nx] - vecin[ii-Nx])/hy/hy/specrad;
				}
			}

			// negative z laplacian
			if(Dims == 3){
				if(zpos == 0){
					vecout[ii] = vecout[ii] + (2.0*vecin[ii] - vecin[ii+Nx*Ny])/hz/hz/specrad;
				}
				else if((zpos+1) == Nz){
					vecout[ii] = vecout[ii] + (2.0*vecin[ii] - vecin[ii-Nx*Ny])/hz/hz/specrad;
				}
				else if(ii < N){
					vecout[ii] = vecout[ii] + (2.0*vecin[ii] - vecin[ii+Nx*Ny] - vecin[ii-Nx*Ny])/hz/hz/specrad;
				}
			}

			// potential
			if(Dims == 1){
				vecout[ii] = vecout[ii] + vecin[ii]*(*potential)(xmin+(double)(xpos+1)*hx,0.0,0.0)/specrad;
			}
			else if(Dims == 2){
				vecout[ii] = vecout[ii] + vecin[ii]*(*potential)(xmin+(double)(xpos+1)*hx,ymin+(double)(ypos+1)*hy,0.0)/specrad;
			}
			else if(Dims == 3){
				vecout[ii] = vecout[ii] + vecin[ii]*(*potential)(xmin+(double)(xpos+1)*hx,ymin+(double)(ypos+1)*hy,zmin+(double)(zpos+1)*hz)/specrad;
			}
		}

	}  // All threads join master thread and disband

}
