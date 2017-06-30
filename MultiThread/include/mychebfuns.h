#ifndef mychebfuns_h__included
#define mychebfuns_h__included
void dct(int N, double *in, double *out);
void chebctor(double (*fun)(double), double a, double b, int *M, double *coeffs, int *flag);
#endif 


