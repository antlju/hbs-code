#include "typedef.h"
#include "diffops.h"
/// First derivatives
Real delx(const Bundle &B, const Real dxfactor, const Int i, const Int vi)
{
        return dxfactor * fd4d1(B(i-2,0,vi),B(i-1,0,vi),
                                B(i+1,0,vi),B(i+2,0,vi));
}

Real dely(const Bundle &B, const Real dyfactor, const Int i, const Int vi)
{
        return dyfactor * fd4d1(B(i,7,vi),B(i,3,vi), //refer to qjkmap.h for q->j vals.
                                B(i,1,vi),B(i,5,vi));
}

Real delz(const Bundle &B, const Real dzfactor, const Int i, const Int vi)
{
        return dzfactor * fd4d1(B(i,6,vi),B(i,2,vi), //refer to qjkmap.h for q->k vals.
                                B(i,4,vi),B(i,8,vi));
}



/// Second derivatives
Real del2x(const Bundle &B, const Real dxfactor, const Int i, const Int vi)
{
        return dxfactor * fd4d2(B(i-2,0,vi),B(i-1,0,vi),
				B(i,0,vi),
                                B(i+1,0,vi),B(i+2,0,vi));
}

Real del2y(const Bundle &B, const Real dyfactor, const Int i, const Int vi)
{
        return dyfactor * fd4d2(B(i,7,vi),B(i,3,vi),
				B(i,0,vi),
                                B(i,1,vi),B(i,5,vi));
}

Real del2z(const Bundle &B, const Real dzfactor, const Int i, const Int vi)
{
        return dzfactor * fd4d2(B(i,6,vi),B(i,2,vi),
				B(i,0,vi),
                                B(i,4,vi),B(i,8,vi));
}

