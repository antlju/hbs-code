#include "includes.h"

void lapl(const Bundle &B, Pencil &P, const Real xfac, const Real yfac, const Real zfac)
{
	for (size_t i=0;i<B.nx_;i++)
	{
		P(i,0,0) = del2x(B,xfac,i,0)+del2y(B,yfac,i,0)+del2z(B,zfac,i,0);
	}
}
