#include <HamOp.h>
#include <ChebPoly.h>

#ifndef ChebPolyOp_h__included
#define ChebPolyOp_h__included
class ChebPolyOp : public ChebPoly, public HamOp
{

	public:
		// constructor
		ChebPolyOp(const ChebPoly &P, const HamOp &H) : ChebPoly(P), HamOp(H){}

		// matrix vector multiplier
		void multChebPolyOp(double *vecin, double *vecout);

		// printer
		void print() const;
};
#endif
