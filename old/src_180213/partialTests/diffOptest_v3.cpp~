#include "common.h"
#include "includes.h"

/// This is a test of the differential operators.
/// Let's start with getting a simple x-derivative working on a bundle.
/// And this seems to work fine!

sBundle sp2sb(const sPencil &p)
{
        sBundle ret(p.nx_);
        for (size_t i=0;i<p.nx_;i++)
        {
                ret(i,0) = p(i);
        }
        return ret;
}

sPencil dxpen(const sBundle &b,Real xfactor)
{
        sPencil ret(b.nx_);
        for (size_t i=0;i<b.nx_;i++)
        {
                ret(i) = xfactor*der1(b(i-2,0),b(i-1,0),b(i+1,0),b(i+2,0));
        }
        return ret;
}

Int main()
{
        Int Nx = NX;
        Int Ny = NY; 
        Int Nz = NZ;
        
        Real L0 = 0.0;
        Real L1 = 2*M_PI;
        Real dx = (L1-L0)/Nx;
        
        sPencil x(Nx); 
        linspace(x,L0,L1,dx); /// x-grid
        x.printpencil();

        /// define function f(x) = lambda*x
        Real lambda = 3.0;
        sPencil fpen = x*lambda;

        fpen.printpencil();
        sBundle fbun = sp2sb(fpen);

        fbun.printpencil(0,0);

        sPencil dfpen = dxpen(fbun,1.0/dx);

        dfpen.printpencil();
        
        return 0;
}
