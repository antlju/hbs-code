#pragma once
#include "typedef.h"

/// Finite differences
Real fd4d1(const Real m2h, const Real m1h, const Real p1h, const Real p2h);
Real fd4d2(const Real m2h, const Real m1h, const Real mid, const Real p1h, const Real p2h);

/// Bundle partial derivatives
Real delx(const Bundle &B, const Real dxfactor, const Int i, const Int vi);
Real dely(const Bundle &B, const Real dyfactor, const Int i, const Int vi);
Real delz(const Bundle &B, const Real dzfactor, const Int i, const Int vi);

/// Curl operator on bundle
void curl(const Bundle &u, Pencil &omega, const Real xfac, const Real yfac, const Real zfac);
