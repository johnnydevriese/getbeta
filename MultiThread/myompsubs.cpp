#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_cblas.h>
#include <HamOp.h>

void printnumthreads(void){

	// thread variables
	int tid, nthds;

	/* Fork a team of threads giving them their own copies of variables */
	#pragma omp parallel private(tid, nthds)
	{
		
		// compute thread variables
		tid = omp_get_thread_num();
		if(tid == 0){
			nthds = omp_get_num_threads();
			printf("\nNumber of Threads = %d\n",nthds);
		}

	}  /* All threads join master thread and disband */

}

void lap(int N, double scl, double *vecin, double *vecout){

	// thread variables
	int nthds, tid;

	// compute variables
	int stride, start, stop, ii;

	/* Fork a team of threads giving them their own copies of variables */
	#pragma omp parallel private(nthds, tid, stride, start, stop, ii) shared(N, scl, vecin, vecout)
	{
		// compute thread variables
		nthds = omp_get_num_threads();
		tid = omp_get_thread_num();

		// compute stride
		stride = ceil((long double)N/nthds);

		// compute start and stop
		start = tid*stride;
		stop = (int)fminl((long double)(tid+1)*stride,(long double)N);

		// print info
		//printf("id, stride, start, stop = %d, %10d, %10d, %10d\n",tid,stride,start,stop);

		// multiply vecin by lap and write to vecout
		for(ii=start;ii<stop;ii++){
			if(ii == 0){
				vecout[ii] = (2.0*vecin[ii] - vecin[ii+1])/scl;
			}
			else if(ii == N-1){
				vecout[ii] = (2.0*vecin[ii] - vecin[ii-1])/scl;
			}
			else if(ii < N){
				vecout[ii] = (2.0*vecin[ii] - vecin[ii+1] - vecin[ii-1])/scl;
			}
		}

	}  /* All threads join master thread and disband */

}

void init(int N, double scl, double *vec){

	// thread variables
	int nthds, tid;

	// compute variables
	int stride, start, stop, ii;

	/* Fork a team of threads giving them their own copies of variables */
	#pragma omp parallel private(nthds, tid, stride, start, stop, ii) shared(N, scl, vec)
	{
		// compute thread variables
		nthds = omp_get_num_threads();
		tid = omp_get_thread_num();

		// compute stride
		stride = ceil((long double)N/nthds);

		// compute start and stop
		start = tid*stride;
		stop = (int)fminl((long double)(tid+1)*stride,(long double)N);

		// print info
		//printf("id, stride, start, stop = %d, %10d, %10d, %10d\n",tid,stride,start,stop);

		// initialize vec
		for(ii=start;ii<stop;ii++){
			vec[ii] = scl;
		}

	}  /* All threads join master thread and disband */

}

double norm(int N, double *vec){

	// thread variables
	int nthds, tid;

	// compute variables
	int m, stride, start, stop;
	double nrm;

	/* Fork a team of threads giving them their own copies of variables */
	#pragma omp parallel private(nthds, tid) shared(m)
	{
		// compute thread variables
		nthds = omp_get_num_threads();
		tid = omp_get_thread_num();
		
		if(tid == 0){
			m = nthds;
		}

	}

	//printf("m = %d\n",m);

	double pnrms[m]; 

	/* Fork a team of threads giving them their own copies of variables */
	#pragma omp parallel private(nthds, tid, stride, start, stop) shared(N, vec, pnrms)
	{
		// compute thread variables
		nthds = omp_get_num_threads();
		tid = omp_get_thread_num();

		// compute stride
		stride = ceil((long double)N/nthds);

		// compute start and stop
		start = tid*stride;
		stop = (int)fminl((long double)(tid+1)*stride,(long double)N);

		pnrms[tid] = cblas_dnrm2(stop-start,&vec[start],1);
		//printf("pnrms[%d] = %+e\n",tid,pnrms[tid]);
	} 

	nrm = cblas_dnrm2(m,&pnrms[0],1);
	//printf("nrm = %+e\n",nrm);

	return nrm;
}

double dot(int N, double *vec1, double *vec2){

	// thread variables
	int nthds, tid;

	// compute variables
	int m, stride, start, stop;
	double dot;

	/* Fork a team of threads giving them their own copies of variables */
	#pragma omp parallel private(nthds, tid) shared(m)
	{
		// compute thread variables
		nthds = omp_get_num_threads();
		tid = omp_get_thread_num();
		
		if(tid == 0){
			m = nthds;
		}

	}

	//printf("m = %d\n",m);

	double pnrms[m]; 

	/* Fork a team of threads giving them their own copies of variables */
	#pragma omp parallel private(nthds, tid, stride, start, stop) shared(N, vec1, vec2, pnrms)
	{
		// compute thread variables
		nthds = omp_get_num_threads();
		tid = omp_get_thread_num();

		// compute stride
		stride = ceil((long double)N/nthds);

		// compute start and stop
		start = tid*stride;
		stop = (int)fminl((long double)(tid+1)*stride,(long double)N);

		pnrms[tid] = cblas_ddot(stop-start,&vec1[start],1,&vec2[start],1);
		//printf("pnrms[%d] = %+e\n",tid,pnrms[tid]);
	} 

	dot = cblas_dasum(m,&pnrms[0],1);
	//printf("nrm = %+e\n",nrm);

	return dot;
}

void scale(int N, double scl, double *vec){

	// thread variables
	int nthds, tid;

	// compute variables
	int stride, start, stop, ii;

	/* Fork a team of threads giving them their own copies of variables */
	#pragma omp parallel private(nthds, tid, stride, start, stop, ii) shared(N, scl, vec)
	{
		// compute thread variables
		nthds = omp_get_num_threads();
		tid = omp_get_thread_num();

		// compute stride
		stride = ceil((long double)N/nthds);

		// compute start and stop
		start = tid*stride;
		stop = (int)fminl((long double)(tid+1)*stride,(long double)N);

		// print info
		//printf("id, stride, start, stop = %d, %10d, %10d, %10d\n",tid,stride,start,stop);

		// initialize vec
		for(ii=start;ii<stop;ii++){
			vec[ii] = scl*vec[ii];
		}

	}  /* All threads join master thread and disband */

}

void scaleadd(int N, double scl, double *vec1, double *vec2){

	// thread variables
	int nthds, tid;

	// compute variables
	int stride, start, stop, ii;

	/* Fork a team of threads giving them their own copies of variables */
	#pragma omp parallel private(nthds, tid, stride, start, stop, ii) shared(N, scl, vec1, vec2)
	{
		// compute thread variables
		nthds = omp_get_num_threads();
		tid = omp_get_thread_num();

		// compute stride
		stride = ceil((long double)N/nthds);

		// compute start and stop
		start = tid*stride;
		stop = (int)fminl((long double)(tid+1)*stride,(long double)N);

		// print info
		//printf("id, stride, start, stop = %d, %10d, %10d, %10d\n",tid,stride,start,stop);

		// initialize vec
		for(ii=start;ii<stop;ii++){
			vec1[ii] = scl*vec1[ii] + vec2[ii];
		}

	}  /* All threads join master thread and disband */

}

void copy(int N, double *vecsrc, double *vectar){

	// thread variables
	int nthds, tid;

	// compute variables
	int stride, start, stop, ii;

	/* Fork a team of threads giving them their own copies of variables */
	#pragma omp parallel private(nthds, tid, stride, start, stop, ii) shared(N, vecsrc, vectar)
	{
		// compute thread variables
		nthds = omp_get_num_threads();
		tid = omp_get_thread_num();

		// compute stride
		stride = ceil((long double)N/nthds);

		// compute start and stop
		start = tid*stride;
		stop = (int)fminl((long double)(tid+1)*stride,(long double)N);

		// print info
		//printf("id, stride, start, stop = %d, %10d, %10d, %10d\n",tid,stride,start,stop);

		// initialize vec
		for(ii=start;ii<stop;ii++){
			vectar[ii] = vecsrc[ii];
		}

	}  /* All threads join master thread and disband */

}

void clenshaw(int N, double coeff, double *vecin, double alpha, double *temp2, double beta, double *vecout, double scl, double *temp1) {

	// thread variables
	int nthds, tid;

	// compute variables
	int stride, start, stop, ii;
	double temp;

	/* Fork a team of threads giving them their own copies of variables */
	#pragma omp parallel private(nthds, tid, stride, start, stop, ii, temp) shared(N, coeff, vecin, alpha, temp2, beta, vecout, scl, temp1)
	{
		// compute thread variables
		nthds = omp_get_num_threads();
		tid = omp_get_thread_num();

		// compute stride
		stride = ceil((long double)N/nthds);

		// compute start and stop
		start = tid*stride;
		stop = (int)fminl((long double)(tid+1)*stride,(long double)N);

		// print info
		//printf("id, stride, start, stop = %d, %10d, %10d, %10d\n",tid,stride,start,stop);

		// initialize vec
		for(ii=start;ii<stop;ii++){
			temp = vecout[ii];
			vecout[ii] = coeff*vecin[ii] + alpha*temp2[ii] + beta*vecout[ii] + scl*temp1[ii];
			temp1[ii] = temp;
		}

	}  /* All threads join master thread and disband */
}
