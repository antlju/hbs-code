#include "common.h"
#include "includes.h"

/// This is a test of the differential operators.
/// x-derivative on a bundle is working (diffOptest_v2.cpp), PBC are working (diffOptest_v3.cpp)
/// Let's now make a general partial derivative on a bundle that gives back a pencil.
/// Let's also base this on a periodic scalar function of space u(x,y,z);

/// THIS VERSION HAS A WORKING X-DERIVATIVE USING FF->BUNDLE-> DBUNDLE->PENCIL->DFF
/// FRAMEWORK!

void u_setup(Field4 &u, const sPencil &x,const sPencil &y,const sPencil &z)
{
        Real l1 = 1.0,l2=1.0,l3=1.0;
        
        for (size_t i = 0;i<u.nx_;i++)
        {
                for (size_t j = 0;j<u.ny_;j++)
                {
                        for (size_t k =0;k<u.nz_;k++)
                        {
                                u(i,j,k) = sin(l1*x(i));//+sin(l2*y(j))+sin(l3*z(k));
                        }
                }
        }
}

void Bpartial_x(const sBundle &B,sPencil &P,const Real xfactor)
{
        for (size_t i=0;i<B.nx_;i++)
        {
                P(i,0) = xfactor * der1(B(i-2,0),B(i-1,0),
                                        B(i+1,0),B(i+2,0));
        }
        
}

void partial_x(const Field4 &ff, Field4 &dff, const Real xfactor, Int vi=0 )
{
        sBundle Dbundle((Int)ff.nx_);
        sPencil Dpencil((Int)ff.nx_);

        //ff2sbundle(ff,Dbundle,0,0);
        
        for (size_t j=0;j<ff.ny_;j++)
        {
                for (size_t k=0;k<ff.nz_;k++)
                {
                        ff2sbundle(ff,Dbundle,j,k);
                        Bpartial_x(Dbundle,Dpencil,xfactor);
                        spencil2ff(Dpencil,dff,j,k);
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
        partial_x(u,du,1.0/dx);

        for (size_t i=0;i<du.nx_;i++)
        {
                std::cout << du(i,0,0) << "\t \t" << cos(x(i)) << "\t \t "
                          << fabs(cos(x(i))-du(i,0,0)) << std::endl;
        }
        
        
        return 0;
}

/// Old functions
/*
void partial_x(const Field4 &ff, Field4 &dff, const Real xfactor, Int vi=0 )
{
                
        for (size_t i=0;i<ff.nx_;i++)
        {
                for (size_t j=0;j<ff.ny_;j++)
                {
                        for (size_t k=0;k<ff.nz_;k++)
                        {

                                dff(i,j,k,vi) = xfactor * der1(ff(i-2,j,k,vi),ff(i-1,j,k,vi),
                                                               ff(i+1,j,k,vi),ff(i+2,j,k,vi));

                                                     
                        }
                }

        }
}
*/

/*
void apply_pbc(Field4 &u)
{
        size_t Ng = u.ng_,Nvars=u.nvar_;
        size_t Nx=u.nx_,Ny=u.ny_,Nz=u.nz_;
        
        for (size_t vi=0;vi<Nvars;vi++)
        {
                
                for (size_t i=0;i<Nx;i++)
                {
                        for (size_t j=0;j<Ny;j++)
                        {
                                for (size_t k=0;k<Nz;k++)
                                {

                                        //set pbc along x
                                        u(-Ng,j,k,vi) = u(Nx-2,j,k,vi);
                                        u(-Ng+1,j,k,vi) = u(Nx-1,j,k,vi);
                                        u(Nx,j,k,vi) = u(0,j,k,vi);
                                        u(Nx+1,j,k,vi) = u(1,j,k,vi);

                                        //set pbc along y
                                        u(i,-Ng,k,vi) = u(i,Ny-2,k,vi);
                                        u(i,-Ng+1,k,vi) = u(i,Ny-1,k,vi);
                                        u(i,Ny,k,vi) = u(i,0,k,vi);
                                        u(i,Ny+1,k,vi) = u(i,1,k,vi);

                                        //set pbc along z
                                        u(i,j,-Ng,vi) = u(i,j,Nz-2,vi);
                                        u(i,j,-Ng+1,vi) = u(i,j,Nz-1,vi);
                                        u(i,j,Nz,vi) = u(i,j,0,vi);
                                        u(i,j,Nz+1,vi) = u(i,j,1,vi);

                                }
                        }
                }
               
        }
}
*/
