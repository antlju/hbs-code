#include "common.h"

void printfield(const Field4 &u) //Prints contents to command line (excluding ghost points)
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
                                        //std::cout << ii << " " << jj << " " << kk << " " << vii << " " << u(ii,jj,kk,vii) << std::endl;
                                        std::cout << u(ii,jj,kk,vii) << "\t";
                                }
                                std::cout << "\n";
                        }
                        std::cout << "---------------------------------------" << std::endl;
                }
                        
        }
}

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

void ABC_setup(Field4 &u, sPencil &x, sPencil &y, sPencil &z)
{
        Real kk=1.0,A=1.0,B=1.0,C=1.0;
        for (size_t i=0;i<u.nx_;i++)
        {
                for (size_t j=0;j<u.ny_;j++)
                {
                        for (size_t k=0;k<u.nz_;k++)
                        {
                                u(i,j,k,0) = A*sin(kk*z(k,0,0))+C*cos(kk*y(j,0,0));
                                u(i,j,k,1) = B*sin(kk*x(i,0,0))+A*cos(kk*z(k,0,0));
                                u(i,j,k,2) = C*sin(kk*y(j,0,0))+B*cos(kk*x(i,0,0));
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
        
        size_t Hsize = 4*NGHOSTS+1;
        
        Real L0 = 0.0;
        Real L1 = 2*M_PI;
        Real dx = (L1-L0)/Nx;
        
        sPencil x(Nx); //x-grid
        x.printpencil(0,0);
        linspace(x,L0,L1,dx);
        x.printpencil(0,0);

        vPencil g(Nx);
        sBundle f(Nx);
        for (size_t q=0;q<Hsize;q++)
        {

                        f.fillPencil(q,q,0);
                        f.printpencil(q,0,"");
                
        }

        Field4 ABC(Nx,Ny,Nz,3);
        ABC_setup(ABC,x,x,x);
        //printfield(ABC); //seems to be setup good.
        return 0;
}