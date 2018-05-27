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

/// Takes a vector bundle and computes P = (B dot grad)B, a vector Pencil.
void udotgradu(const Bundle &B, Pencil &P, const Real xfac, const Real yfac, const Real zfac)
{
        assert((B.nx_ == P.nx_) && (B.nvars_ == P.nvars_) );
        for (size_t i=0;i<B.nx_;i++)
        {
                for (size_t vi=0;vi<P.nvars_;vi++)
                {
                        P(i,0,vi) = B(i,0,0)*delx(B,xfac,i,vi)
                                +B(i,0,1)*dely(B,yfac,i,vi)
                                +B(i,0,2)*delz(B,zfac,i,vi);
                }
        }
}
