#include "common.h"
#include "includes.h"

/// Test program for computing the curl of the ABC-field.
/// - version 1.

/// Takes two vector fields and produces a scalar function.
void dotprod(const Field4 &A, const Field4 &B, Field4 &C)
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

void setup(Field4 &u, const sPencil &x, const sPencil &y, const sPencil &z)
{
        Real A=1.0,B=1.0,C=1.0;
        Real L=1.0;
        for (size_t i=0;i<u.nx_;i++)
        {
                for (size_t j=0;j<u.ny_;j++)
                {
                        for (size_t k=0;k<u.nz_;k++)
                        {
                                u(i,j,k,0) = A*sin(L*z(k))+C*cos(L*y(j));
                                u(i,j,k,1) = B*sin(L*x(i))+A*cos(L*z(k));
                                u(i,j,k,2) = C*sin(L*y(j))+B*cos(L*x(i));
                                //u(i,j,k,0) = x[i];
                                //u(i,j,k,1) = y[i];
                                //u(i,j,k,2) = z[i];
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
        //x.printpencil();

        Field4 u(Nx,Ny,Nz,3);
        setup(u,x,x,x);
        apply_pbc(u);
        //printfield(u);
        
        Field4 cu(Nx,Ny,Nz,3);
        curl(u,cu,1.0/dx,1.0/dx,1.0/dx); // Curl is NOT WORKING!
        // GIVES SOME WEIRD THINGS!

        Field4 udotu(Nx,Ny,Nz,1);
        Field4 udotcu(Nx,Ny,Nz,1);
        dotprod(u,u,udotu); //CORRECT!
        dotprod(u,cu,udotcu);


        /// PRINTING FOR DEBUGGING! Ideally we want to this to print udotcu/udotu = k everywhere.
        for (size_t vii=0;vii<udotu.nvar_;vii++)
        {
                std::cout << "------- component: " << vii << " -------" << std::endl;
                std::cout << "---------------------------------------" << std::endl;
                for (size_t ii=0;ii<udotu.nx_;ii++)
                {
                        for (size_t jj=0;jj<udotu.ny_;jj++)
                        {
                                for (size_t kk=0;kk<udotu.nz_;kk++)
                                {
                                        std::cout << cu(ii,jj,kk,0) << "\t";
                                }
                                std::cout << "\n";
                        }
                        std::cout << "---------------------------------------" << std::endl;
                        std::cout << std::endl;
                }
                        
        }
        
        
        return 0;
}
