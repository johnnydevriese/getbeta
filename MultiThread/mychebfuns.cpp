#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fftw3.h>
#include <math.h>
#include <unistd.h>
#include <omp.h>

#define PI 3.1415926535897932

void dct(int N, double *in, double *out){

	// compute variables
	int ii;
	fftw_plan my_plan;
	fftw_init_threads();

	// define plan
	fftw_plan_with_nthreads(omp_get_max_threads());
	my_plan = fftw_plan_r2r_1d(N, in, out, FFTW_REDFT00, FFTW_ESTIMATE);

	//execute plan
	fftw_execute(my_plan);

	// scale output
	for(ii=0; ii < N; ii++){
		if(ii == 0 || ii == N-1){
			out[ii] = out[ii]/(double)(N-1)/2.0;
		}
		else{
			out[ii] = out[ii]/(double)(N-1);
		}
	}

	// destroy plan
	fftw_destroy_plan(my_plan);
	fftw_cleanup_threads();

}

void chebctor(double (*fun)(double), double a, double b, int *M, double *coeffs, int *flag){

	// compute variables
	int sz = 17;
	int ii, jj, N = pow(2,sz)+1, stride;
	double *x,*temp, maxnum, tol = 2.2204e-16;
	double pt;

	// allocate memory
	temp = (double*)malloc(sizeof(double)*N);
	x = (double*)malloc(sizeof(double)*N);

	// initialize 
	*flag = 1;

	// compute shifted chebyshev pts
	for(ii=0; ii < N; ii++){
		pt = cos((long double)(ii)*PI/(long double)(N-1));
		x[ii] = ((b-a)*pt+(b+a))/2.0;
	}

	// compute y = f(x) and maxnum
	maxnum = 0.0;
	for(ii=0; ii < N; ii++){
		x[ii] = (*fun)(x[ii]);
		if(fabs(x[ii]) > maxnum){
			maxnum = fabs(x[ii]);
		}
	}

	// exit if overflow
	if(maxnum == INFINITY){
		printf("Function causes overflow!\n");
		*flag = 2;
		return;
	}
	tol = tol*maxnum;

	// compute coefficients adaptively
	for(ii=3; ii<sz+1; ii++){
		// set current degree
		*M = pow(2,ii)+1;

		// set stride
		stride = pow(2,sz-ii);

		// fill temp
		for(jj=0;jj<*M;jj++){
			temp[jj] = x[jj*stride];
		}

		// compute coeffs
		dct(*M, temp, coeffs);

		// check for convergence
		for(jj=0;jj<*M;jj++){
			if(fabs(coeffs[*M-1-jj]) >= tol){
				if(jj > 7){
					*M = *M-1-jj;
					*flag = 0;

					// store coeffs in correct order
					for(jj=0;jj<*M+1;jj++){
						temp[jj] = coeffs[*M-jj];
					}
					for(jj=0;jj<*M+1;jj++){
						coeffs[jj] = temp[jj];
					}
					for(jj=*M+1;jj<N;jj++){
						coeffs[jj] = 0.0f;
					}

					// free memory
					free(temp);
					free(x);

					return;
				}
				break;
			}
		}
	}

	// free memory
	free(temp);
	free(x);

}
