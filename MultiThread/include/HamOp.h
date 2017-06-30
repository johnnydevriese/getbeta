#ifndef HamOp_h__included
#define HamOp_h__included
double zero(double x, double y, double z);
double one(double x, double y, double z);
double wr2(double x, double y, double z);

class HamOp{
	private:
		// private variables
		int Dims;
		int Nx;
		int Ny;
		int Nz;
		double xmin;
		double xmax;
		double ymin;
		double ymax;
		double zmin;
		double zmax;
		double (*potential)(double x, double y, double z);

	public:
		// public variables
		double specrad;

		// getinfo functions
		int getDims() {return Dims;}
		int getNx() {return Nx;}
		int getNy() {return Ny;}
		int getNz() {return Nz;}
		double getxmin() {return xmin;}
		double getxmax() {return xmax;}
		double getymin() {return ymin;}
		double getymax() {return ymax;}
		double getzmin() {return zmin;}
		double getzmax() {return zmax;}
		double getspecrad() {return specrad;}

		// set info functions
		void setspecrad(double num) {specrad = num;}
	
		// compute specrad using arpack
		void computespecrad();

		// matrix vector multiplier
		void multHamOp(double *vecin, double *vecout);

		// printer
		void print() const;

		// constructor
		HamOp(int i1=1, int i2=1, int i3=1, int i4=1, double d1=0.0, double d2=1.0, double d3=0.0, double d4=1.0, 
			double d5=0.0, double d6=1.0, double (*d7)(double , double , double )=&zero, double d8=1.0) :
			Dims(i1), Nx(i2), Ny(i3), Nz(i4), xmin(d1), xmax(d2), ymin(d3), ymax(d4), zmin(d5), zmax(d6), potential(d7), specrad(d8) {}

		// destructor
		~HamOp(){}
};
#endif
