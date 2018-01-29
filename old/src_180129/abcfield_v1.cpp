#include "common.h"
#include "utils.h"

/// Utility methods
void printfield(const Rfield &u,const Int useNg = 0) //Prints contents to command line (excluding ghost points)
{
        if (useNg !=0)
        {
                for (size_t vii=0;vii<u.nvar_;vii++)
                {
                        std::cout << "------- component: " << vii << " -------" << std::endl;
                        std::cout << "---------------------------------------" << std::endl;
                        for (size_t ii=0;ii<u.nx_+u.ng_;ii++)
                        {
                                for (size_t jj=0;jj<u.ny_+u.ng_;jj++)
                                {
                                        for (size_t kk=0;kk<u.nz_+u.ng_;kk++)
                                        {
                                                std::cout << ii << " " << jj << " " << kk << " " << vii << " " << u(ii-u.ng_,jj-u.ng_,kk-u.ng_,vii) << std::endl;
                                                //std::cout << u(ii,jj,kk,vii) << "\t";
                                        }
                                        std::cout << "\n";
                                }
                                std::cout << "---------------------------------------" << std::endl;
                        }
                        
                }
        }
        else
        { 
                for (size_t vii=0;vii<u.nvar_;vii++)
                {
                        std::cout << "------- component: " << vii << " -------" << std::endl;
                        std::cout << "---------------------------------------" << std::endl;
                        for (size_t ii=0;ii<u.nx_;ii++)
                        {
                                for (size_t jj=0;jj<u.ny_;jj++)
                                {
                                        for (size_t kk=0;kk<u.nz_;kk++)
                                        {
                                                std::cout << ii << " " << jj << " " << kk << " " << vii << " " << u(ii,jj,kk,vii) << std::endl;
                                                //std::cout << u(ii,jj,kk,vii) << "\t";
                                        }
                                        std::cout << "\n";
                                }
                                std::cout << "---------------------------------------" << std::endl;
                        }
                        
                }
        }
}
                
Real *linspace(Int size, Real start, Real end)
{

        Real dx = (end-start)/size;
        Real *axis = new Real[size];
        for (Int i=0;i<size;i++)
        {
                axis[i] = start+dx*i;
                std::cout << axis[i] << std::endl;
        }
        std::cout << std::endl;
        return axis;

}

void ABC_setup(Rfield &u, const Real kk, const Real A, const Real B, const Real C, const Real *x, const Real *y, const Real *z)
{
        for (size_t i=0;i<u.nx_;i++)
        {
                for (size_t j=0;j<u.ny_;j++)
                {
                        for (size_t k=0;k<u.nz_;k++)
                        {
                                u(i,j,k,0) = A*sin(kk*z[k])+C*cos(kk*y[j]);
                                u(i,j,k,1) = B*sin(kk*x[i])+A*cos(kk*z[k]);
                                u(i,j,k,2) = C*sin(kk*y[j])+B*cos(kk*x[i]);
                                //u(i,j,k,0) = x[i];
                                //u(i,j,k,1) = y[i];
                                //u(i,j,k,2) = z[i];
                        }
                }
        }

}

Int main()
{
        Int Nx=NX,Ny=NY,Nz=NZ;

        /// ABC-field parameters.
        Real kk=1.0,A=1.0,B=1.0,C=1.0;

        /// Construct space grids
        Real L0 = 0.0;
        Real L1 = 2*M_PI;
        Real dx = (L1-L0)/Nx;

        Real *x = linspace(Nx,L0,L1);
        Real *z = linspace(Nz,L0,L1);
        Real *y = linspace(Ny,L0,L1);
        
/// Introduce the vectorfield object with zeroed entries.
        Rfield uABC(Nx,Ny,Nz,3);

        ABC_setup(uABC,kk,A,B,C,x,y,z);
        //printfield(uABC);
        apply_pbc(uABC);

        return 0;
}
