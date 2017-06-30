#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <myompsubs.h>
#include <myfortompsubs.h>
#include <HamOp.h>
#include <ChebPoly.h>
#include <ChebPolyOp.h>

void ChebPolyOp::print() const{
	printf("\nChebPolyOp\n");
	HamOp::print();
	ChebPoly::shortprint();
}

void ChebPolyOp::multChebPolyOp(double *vecin, double *vecout) {

	// compute variables
	int ii, N;
	int sd = sizeof(double);
	double beta, alpha;
	double *temp1, *temp2;

	// N
	N = HamOp::getNx()*HamOp::getNy()*HamOp::getNz();

	// alpha and beta
	alpha = 2.0/(ChebPoly::getb()-ChebPoly::geta());
	beta = -(ChebPoly::getb()+ChebPoly::geta())/(ChebPoly::getb()-ChebPoly::geta());

	// allocate memory
	temp1 = (double*)malloc(N*sd); if(temp1 == NULL){printf("Out of memory\n"); exit(-1);}
	temp2 = (double*)malloc(N*sd); if(temp2 == NULL){printf("Out of memory\n"); exit(-1);}

	// initialize vecin
	init(N, 0.0, temp1);
	init(N, 0.0, temp2);

	// loop to evaluate p(A)v
	if(ChebPoly::getM() == 0){
		clenshaw(N, ChebPoly::getcoeffs()[0], vecin, 0.0, temp2, 0.0, vecout, 0.0, temp1);
	}
	else{
		clenshaw(N, ChebPoly::getcoeffs()[0], vecin, 0.0, temp2, 0.0, vecout, 0.0, temp1);
		init(N, 0.0, temp1);

		for(ii = 0; ii < (ChebPoly::getM()-1); ii++){
			HamOp::multHamOp(vecout, temp2);
			clenshaw(N, ChebPoly::getcoeffs()[ii+1], vecin, 2.0*alpha, temp2, 2.0*beta, vecout, -1.0, temp1);
		}

		HamOp::multHamOp(vecout, temp2);
		clenshaw(N, ChebPoly::getcoeffs()[ChebPoly::getM()], vecin, alpha, temp2, beta, vecout, -1.0, temp1);
	}

	// free memory
	free(temp1);
	free(temp2);
}
