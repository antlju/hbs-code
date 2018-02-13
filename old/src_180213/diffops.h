#pragma once
#include "common.h"

/// DERIVATIVES
Real der1(const Real m2h, const Real m1h, const Real p1h, const Real p2h);
Real der2(const Real m2h, const Real m1h, const Real mid, const Real p1h, const Real p2h);

/// GRADIENT
void grad(const sBundle &f, vPencil &g);

/// CURL
Real pointPartial(const vBundle &B,const Int i,const Real dfactor,const Int mu, const Int vi);
void bundlecurl(const vBundle &u, sPencil &cu, Int vi, const Real dfactor1, const Real dfactor2);
void curl(const Field4 &u, Field4 &cu, const Real xfactor, const Real yfactor, const Real zfactor);
