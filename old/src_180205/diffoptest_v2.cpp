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

Int main()
{
        Int Nx = NX;
        //Int Ny = NY;
        //Int Nz = NZ;
        
        Real L0 = 0.0;
        Real L1 = 2*M_PI;
        Real dx = (L1-L0)/Nx;
        
        sPencil x(Nx); //x-grid
        x.printpencil(0,0);
        linspace(x,L0,L1,dx);
        x.printpencil(0,0);

        vPencil g(Nx);
        sBundle f(Nx);
        //f.printpencil(0,0);

        createqtojkmap();
        

        return 0;
}
