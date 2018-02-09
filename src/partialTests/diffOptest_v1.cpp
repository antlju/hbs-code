#include "common.h"
#include "includes.h"


/// THIS VERSION (V1) IS NOT DOING AT ALL WHAT IT SHOULD
/// GETTING WEIRD VALUES IN DERIVATIVE PENCILS.
/// F.EX IF FIELD IS EXP(X) THEN PENCIL CORRESPONDING TO d_x(exp(x)) IS
/// NOT EXP(X) AGAIN! THIS DISCREPANCY IS NOT ONLY CLOSE TO "EDGES".



void expfieldsetup(Field4 &u,sPencil &x, sPencil &y, sPencil &z)
{
        Real l1 = 1.0,l2=1.0,l3=1.0;
        
        for (size_t i = 0;i<u.nx_;i++)
        {
                for (size_t j = 0;j<u.ny_;j++)
                {
                        for (size_t k =0;k<u.nz_;k++)
                        {
                                //u(i,j,k) = exp(l1*x(i))+exp(l2*y(j))+exp(l3*z(k));
                                u(i,j,k) = exp(l1*x(i));
                        }
                }
        }
}

void bpartial(const sBundle &u, sPencil &du, Int mu)
{
        assert(mu == 0 || mu == 1 || mu == 2);
        if (mu == 0)
        {
                for (size_t i=0;i<du.nx_;i++)
                {
                        du(i,0) = der1(u(i-2,0,0),u(i-1,0,0),u(i+1,0,0),u(i+2,0,0));
                }
        }
        
}

void ffpartial(const Field4 &u, Field4 &du, const Int mu, const Int vi=0)
{
        sBundle dbundle(u.nx_);
        sPencil dpencil(u.nx_); /// The result of differentiation.
        
        for (size_t j=0;j<u.ny_;j++)
        {
                for (size_t k=0;k<u.nz_;k++)
                {

                        ff2sbundle(u,dbundle,j,k);
                        bpartial(dbundle,dpencil,mu);
                        dpencil.printpencil();
                        spencil2ff(dpencil,du,j,k);
                }
        }
}

/// This is a test of the differential operators.
/// 
///
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

        /// Construct a scalar field.
        Field4 sff(Nx,Ny,Nz,1);

        /// Let's try calculating derivatives on a scalar function
        /// u(x,y,z) = sum_k=0^2 exp^(lambda_k x_k).

        expfieldsetup(sff,x,x,x);

        printfield(sff);

        Field4 dsff(Nx,Ny,Nz,1);
        ffpartial(sff,dsff,0,0);

        printfield(pwiseffdiv(dsff,sff));
        
        return 0;
}
