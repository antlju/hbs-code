#include "common.h"

void linspace(sPencil &u, Real start, Real end, Real dx)
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

void printvpencil(vPencil &u, const std::string title="")
{
        std::cout << "------------ spencil print: " << title << std::endl;
        for (size_t i=0;i<u.nx_;i++)
                std::cout << u(i,0,0) << " " << u(i,0,1) << " " << u(i,0,2) << std::endl;
        
        std::cout << std::endl;
        //std::cout << "----------------------------" << std::endl;
}

void printvbundle(vBundle &u, Int q, const std::string title="")
{
        std::cout << "------------ sbundle print: " << title << std::endl;
        for (size_t i=0;i<u.nx_;i++)
                std::cout << u(i,q,0) << " " << u(i,q,1) << " " << u(i,q,2) << " " << u.nvars_ << std::endl;
        
        std::cout << std::endl;
}


Int main()
{
        Int Nx = NX;

        Real L0 = 0.0;
        Real L1 = 2*M_PI;
        Real dx = (L1-L0)/Nx;
        
        sPencil x(Nx);
        vPencil a(Nx);
        a.fillPencil(3.0);
        
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

        vBundle B(Nx);

        //B.fillPencil(3.0, 0);
        //B.fillPencil(2.0, 1);
        //B.fillPencil(1.0, 2);

        //printsbundle(B,0);

        for (size_t qi=0;qi<5;qi++)
        {
                B.fillPencil(1.0+qi,qi,0);
                B.fillPencil(1.0+qi,qi,1);
                B.fillPencil(1.0+qi,qi,2);
        }
        
        for (size_t j =0;j<5;j++)
                printvbundle(B,j);

        B = B-(B*2.0);

        for (size_t j =0;j<5;j++)
                printvbundle(B,j);
        
        return 0;
}
