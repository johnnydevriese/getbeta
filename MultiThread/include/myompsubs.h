#ifndef myompsubs_h__included
#define myompsubs_h__included
void printnumthreads(void);
void lap(int N, double scl, double *vecin, double *vecout);
void init(int N, double scl, double *vec);
void clenshaw(int N, double coeff, double *vecin, double alpha, double *temp2, double beta, double *vecout, double scl, double *temp1);
void scale(int N, double scl, double *vec);
void copy(int N, double *vecsrc, double *vectar);
double norm(int N, double *vec);
double dot(int N, double *vec1, double *vec2);
void scaleadd(int N, double scl, double *vec1, double *vec2);
#endif 


