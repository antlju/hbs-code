#include "common.h"
#include "includes.h"

/// This is a test of the differential operators.
/// x-derivative on a bundle is working (diffOptest_v2.cpp)
/// Let's now make a general partial derivative on a bundle that gives back a pencil.
/// Let's also base this on a periodic scalar function of space u(x,y,z);

/// THIS VERSION I HAVE TRIED THE PBC AND IT WORKS!
///

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

void bpartial_x(const sBundle &u, sPencil &du, const Real xfactor)
{

        for (size_t i=0;i<du.nx_;i++)
        {
                du(i,0) = xfactor*der1(u(i-2,0,0),u(i-1,0,0),u(i+1,0,0),u(i+2,0,0));
        }
        
}

Field4 partial_x(const Field4 &u, const Real xfactor)
{
        sBundle dbundle(u.nx_);
        sPencil dpencil(u.nx_); /// The result of differentiation.
        Field4 du(u.nx_,u.ny_,u.nz_,u.nvar_);
        
        for (size_t j=0;j<u.ny_;j++)
        {
                for (size_t k=0;k<u.nz_;k++)
                {
                        ff2sbundle(u,dbundle,j,k);
                        bpartial_x(dbundle,dpencil,xfactor);
                        spencil2ff(dpencil,du,j,k);
                }
        }

        return du;
}

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
        x.printpencil();

        Field4 u(Nx,Ny,Nz,1);
        u_setup(u,x,x,x);
        for (size_t k=0;k<u.nz_+2*u.ng_;k++)
                std::cout << u(k-u.ng_,0,0) << std::endl;

        std::cout << " --------------- " << std::endl;
        
        apply_pbc(u);

        for (size_t k=0;k<u.nz_+2*u.ng_;k++)
                std::cout << u(k-u.ng_,0,0) << std::endl;
        
        std::cout << " --------------- " << std::endl;
        
        
        /*
        Field4 du = partial_x(u,1.0/dx);
        
        for (size_t i = 0;i<du.nx_;i++)
        {
                std::cout << fabs(du(i,0,0)-cos(x(i))) << std::endl;
        }
        */
        
        return 0;
}
