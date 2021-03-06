#include <myompsubs.h>

#ifndef ChebPoly_h__included
#define ChebPoly_h__included
class ChebPoly{
	private:
		// private variables
		int M;
		double a;
		double b;
		double *coeffs;

/*		// set info functions
		void setM(int num) {M = num;}
		void setM(double num) {a = num;}
		void setM(double num) {b = num;}
		void setM(double *num) {coeffs = num;}
*/
	public:
		// public variables

		// getinfo functions
		int getM() {return M;}
		double geta() {return a;}
		double getb() {return b;}
		double *getcoeffs() {return coeffs;}

		// printer
		void print() const;
		void shortprint() const;

		// constructors
		//ChebPoly(const ChebPoly &other) : M((other).M), a((other).a), b((other).b), coeffs((other).coeffs) {} 
		ChebPoly(int i1=0, double d1=-1.0, double d2=1.0) : M(i1), a(d1), b(d2){
			coeffs = (double*)malloc((pow(2,17)+1)*sizeof(double)); if(coeffs == NULL){printf("Out of memory for coeffs in ChebPoly!\n"); exit(-1);}
			init(pow(2,17)+1, 0.0, coeffs);
		}
		ChebPoly(int i1, double d1, double d2, double *d3) : M(i1), a(d1), b(d2), coeffs(d3){}

		// destructor
		~ChebPoly(){}
};
#endif
