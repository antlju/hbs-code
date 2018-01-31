#include "common.h"

void linspace(spencil &u, Real start, Real end, Real dx)
{
        Real size = u.nx_;
        //Real dx = (end-start)/size;
        for (size_t i=0;i<size;i++)
        {
                u(i) = start+dx*i;
                //std::cout << u(i) << std::endl;
        }
        //std::cout << std::endl;
}

void printspencil(spencil &u,const std::string title="")
{
        std::cout << "------------ spencil print: " << title << std::endl;
        for (size_t i=0;i<u.nx_;i++)
                std::cout << u(i) << std::endl;
        
        std::cout << std::endl;
        //std::cout << "----------------------------" << std::endl;
}
Int main()
{
        Int Nx = NX;

        Real L0 = 0.0;
        Real L1 = 2*M_PI;
        Real dx = (L1-L0)/Nx;
        
        spencil x(Nx);
        spencil a(Nx);
        a.fillUp(3.0);
        
        linspace(x,L0,L1,dx);

        /// some nice tests of arithmetics
        /*
        printspencil(x, "x");
        x = x*6.0;
        printspencil(x, "6.0*x");

        printspencil(a, "a");
        a = a-x;
        printspencil(a, "a-x");
        */
        
        return 0;
}
