#include "includes.h"

/// Divergence of a vector field p = div(u).
/// p is a scalar pencil and u is a 3-vector bundle.
void div(const Bundle &u, Pencil &p, const Real xfac, const Real yfac, const Real zfac)
{
	
        for (size_t i = 0;i<u.nx_;i++)
        {
                p(i,0,0) = delx(u,xfac,i,0)+dely(u,yfac,i,1)+delz(u,zfac,i,2);
        }
}
