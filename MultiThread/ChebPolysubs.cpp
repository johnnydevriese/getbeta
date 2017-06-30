#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ChebPoly.h>


// function to print a ChebPoly struct
void ChebPoly::print() const{
	int ii;

	printf("\nChebPoly:\n");
	printf(" M = %d\n",M);
	printf(" a = %+f\n",a);
	printf(" b = %+f\n",b);
	for(ii=0;ii<M+1;ii++){
		printf(" coeffs[%d] = %+e\n",ii,coeffs[ii]);
	}
	printf("\n");
}

// function to print a ChebPoly struct
void ChebPoly::shortprint() const{

	printf("\nChebPoly:\n");
	printf(" M = %d\n",M);
	printf(" a = %+f\n",a);
	printf(" b = %+f\n\n",b);
}
