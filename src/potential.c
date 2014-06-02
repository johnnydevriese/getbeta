#include <math.h>

//creating the potential operator 

//flag should be 1 here because this should be added to the value of Y now because lapmult will have already acted on Y. So, you do not want to overwrite Y. 

//computing h, the size of the interval, locally in this function. 

//beta is computed in HamOp function and is an arbitrary scalar. 

void potential(double beta, int n, double a, double b, double*X, double*Y,int flag){

	int jj;
	double h = (a - b)/(n+1.0);

	if(flag == 0){	
	//jj+1 because it picks the first interval + 1
	for(jj=0;jj<n;jj++){
		Y[jj] = (pow(( b + (jj+1.0)*h),2) * X[jj])*beta; 
	}
	
	}
	else{  //if else do this
	for(jj=0;jj<n;jj++){
		Y[jj] += (pow(( b +(jj+1.0)*h),2) * X[jj])*beta;
	}

	}
 

}

