#include "includes.h"

void sgrad(const Bundle &B,Pencil &P, const Real xfac, const Real yfac, const Real zfac)
{
	for (size_t i=0;i<B.nx_;i++)
	{
		P(i,0,0) = delx(B,xfac,i,0);
		P(i,0,1) = dely(B,yfac,i,0);
		P(i,0,2) = delz(B,zfac,i,0);
	}
}
