#ifndef _DIFFOPS_H_
#define _DIFFOPS_H_

Real der1(Real m2h, Real m1h, Real p1h, Real p2h);
Real der2(Real m2h, Real m1h, Real mid, Real p1h, Real p2h);
void curl(Rfield omega, const Rfield U, const Real xfactor, const Real yfactor, const Real zfactor);
void lapl(Rfield d2f, const Rfield f, const Real xfactor, const Real yfactor, const Real zfactor);
void div(Rfield g, const Rfield u, const Real xfactor, const Real yfactor, const Real zfactor);
void grad(Rfield u, const Rfield f, const Real xfactor, const Real yfactor, const Real zfactor);

#endif
