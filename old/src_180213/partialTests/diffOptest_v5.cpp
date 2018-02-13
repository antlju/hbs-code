#include "common.h"
#include "includes.h"


/// WORKING PARTIALS IN ALL DIRECTIONS!

/// This is a test of the differential operators.
/// x-derivative on a bundle is working (diffOptest_v2.cpp), PBC are working (diffOptest_v3.cpp)
/// x-derivative is working ff->bundle->pencil->dff framwork (diffOptest_v4.cpp)

/// Let's now make a general partial derivative on a bundle that gives back a pencil.
/// Let's also base this on a periodic scalar function of space u(x,y,z);


void u_setup(Field4 &u, const sPencil &x,const sPencil &y,const sPencil &z)
{
        Real l1 = 1.0,l2=1.0,l3=1.0;
        
        for (size_t i = 0;i<u.nx_;i++)
        {
                for (size_t j = 0;j<u.ny_;j++)
                {
                        for (size_t k =0;k<u.nz_;k++)
                        {
                                u(i,j,k) = sin(l1*x(i))+sin(l2*y(j))+sin(l3*z(k));
                        }
                }
        }
}

void vBpartial(const vBundle &B,sPencil &P,const Real dfactor,Int mu, Int vi)
{
        assert(mu>=0 && mu<3);
        if (mu == 0) //x-derivative
        {
                for (size_t i=0;i<B.nx_;i++)
                {
                        P(i,0) = dfactor * der1(B(i-2,0,vi),B(i-1,0,vi),
                                                B(i+1,0,vi),B(i+2,0,vi));
                }
        }
        else if (mu == 1) //y-derivative
        {
                for (size_t i=0;i<B.nx_;i++)
                {
                        P(i,0) = dfactor * der1(B(i,7,vi),B(i,3,vi), //refer to qjkmap.h for q->j vals.
                                                B(i,1,vi),B(i,5,vi));
                }
        }
        else //z-derivative
        {
                for (size_t i=0;i<B.nx_;i++)
                {
                        P(i,0) = dfactor * der1(B(i,6,vi),B(i,2,vi), //refer to qjkmap.h for q->k vals.
                                                B(i,4,vi),B(i,8,vi));
                }
        }
}

void partial(const Field4 &ff, Field4 &dff, const Real dfactor, Int mu, Int vi=0 )
{
        vBundle Dbundle((Int)ff.nx_);
        sPencil Dpencil((Int)ff.nx_);

        //ff2sbundle(ff,Dbundle,0,0);
        
        for (size_t j=0;j<ff.ny_;j++)
        {
                for (size_t k=0;k<ff.nz_;k++)
                {
                        ff2vbundle(ff,Dbundle,j,k,vi);
                        vBpartial(Dbundle,Dpencil,dfactor,mu,vi);
                        spencil2ff(Dpencil,dff,j,k,vi);
                }
        }

}

Int main()
{
        Int Nx = NX;
        Int Ny = NY; 
        Int Nz = NZ;
        
        Real L0 = 0.0;
        Real L1 = 2*M_PI;
        Real dx = (L1-L0)/(Nx);
        
        sPencil x(Nx); 
        linspace(x,L0,L1,dx); /// x-grid
        //x.printpencil();

        Field4 u(Nx,Ny,Nz,1);
        u_setup(u,x,x,x);
        apply_pbc(u);

        Field4 du(Nx,Ny,Nz,1);
        partial(u,du,1.0/dx,2,0); //partial_z

        printfield(du);
        /*
        for (size_t i=0;i<du.nx_;i++)
        {
                std::cout << du(0,i,0) << "\t \t" << cos(x(i)) << "\t \t "
                          << fabs(cos(x(i))-du(0,i,0)) << std::endl;
        }
        */
        for (size_t i=0;i<du.nx_;i++)
        {
                std::cout <<  cos(x(i)) << "\t";
        }
        std::cout << std::endl;
        
        return 0;
}
