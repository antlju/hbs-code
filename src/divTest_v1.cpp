/// Library 
#include <iostream>
#include <cmath>
#include <cassert>

/// Headers
#include "inputparams.h"
#include "vartypedef.h"
#include "typedef.h"
#include "fmesh.h"
#include "pmesh.h"
#include "utility.h"
#include "diffops.h"


/// Test program for computing the curl of the ABC-field.
/// - version 1.

/// Takes two vector fields and produces a scalar function.
void dotprod(const Mesh &A, const Mesh &B, Mesh &C)
{

        for (size_t i=0;i<A.nx_;i++)
        {
                for (size_t j=0;j<A.ny_;j++)
                {
                        for (size_t k=0;k<A.nz_;k++)
                        {
                                
                                C(i,j,k) = A(i,j,k,0)*B(i,j,k,0)
                                        +A(i,j,k,1)*B(i,j,k,1)
                                        +A(i,j,k,2)*B(i,j,k,2);
                        }
                }
        }

}

void setup(Mesh &u, const Pencil &x, const Pencil &y, const Pencil &z)
{
        //Real A=1.0,B=1.0,C=1.0;
        Real L=1.0;
        for (size_t i=0;i<u.nx_;i++)
        {
                for (size_t j=0;j<u.ny_;j++)
                {
                        for (size_t k=0;k<u.nz_;k++)
                        {
                                u(i,j,k,0) = 3*x(i);//sin(x(i));
                                u(i,j,k,1) = 2*y(j);//sin(y(j));
                                u(i,j,k,2) = 1*z(k);//sin(z(k));
                                //u(i,j,k,0) = x[i];
                                //u(i,j,k,1) = y[i];
                                //u(i,j,k,2) = z[i];
                        }
                }
        }

}

/// Divergence of a vector field p = div(u).
/// p is a scalar pencil and u is a 3-vector bundle.
void div(const Bundle &u, Pencil &p, const Real xfac, const Real yfac, const Real zfac)
{
        for (size_t i = 0;i<u.nx_;i++)
        {
                p(i,0,0) = delx(u,xfac,i,0)+dely(u,yfac,i,1)+delz(u,zfac,i,2);
        }
}

void apply_div(const Mesh &u, Mesh &divu, const Real xfac, const Real yfac, const Real zfac)
{
        Bundle Dbundle(u.nx_,3); /// 3-vector bundle
        Pencil Dpencil(u.nx_,1); /// Scalar pencil
        for (size_t j=0;j<u.ny_;j++)
        {
                for (size_t k=0;k<u.nz_;k++)
                {
                        ff2bundle(u,Dbundle,j,k);
                        div(Dbundle,Dpencil,xfac,yfac,zfac);
                        pencil2ff(Dpencil,divu,j,k);
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
        Real dx = (L1-L0)/Nx;
        
        Pencil x(Nx);
        x.printpencil();
        linspace(x,L0,L1,dx);
        x.printpencil();

        Mesh u(Nx,Ny,Nz,3);
        setup(u,x,x,x);
        //apply_pbc(u);
        
        Mesh g(Nx,Ny,Nz,1);

        apply_div(u,g,1.0/dx,1.0/dx,1.0/dx);

        for (Int i = 0;i<g.nx_;i++)
                std::cout << g(i,0,0) << "\t" << g(0,i,0) << "\t" << g(0,0,i) << std::endl;
        /*
        for (size_t vii=0;vii<g.nvar_;vii++)
        {
                std::cout << "------- component: " << vii << " -------" << std::endl;
                std::cout << "---------------------------------------" << std::endl;
                for (size_t ii=0;ii<g.nx_;ii++)
                {
                        for (size_t jj=0;jj<g.ny_;jj++)
                        {
                                for (size_t kk=0;kk<g.nz_;kk++)
                                {
                                        //std::cout << ii << " " << jj << " " << kk << " " << vii << " " << u(ii,jj,kk,vii) << std::endl;
                                        std::cout << g(ii,jj,kk,vii)-(cos(x(jj))+cos(x(jj))+cos(x(kk))) << "\t";
                                }
                                std::cout << "\n";
                        }
                        std::cout << "---------------------------------------" << std::endl;
                        std::cout << std::endl;
                }
                        
        }
        
        */
        
        return 0;
}
