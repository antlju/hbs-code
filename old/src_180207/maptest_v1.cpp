#include "common.h"
#include "qjkmap.h"

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

void ff2sbundle(const Field4 &ff, sBundle &B, Int j, Int k)
{
        for (Int q=0;q<(Int)B.stencilSize_;q++)
        {
                Intpair jk = qtojk(q);
                for (Int i=0;i<(Int)B.nx_;i++)
                {
                        B(i,q,0) = ff(i,j+jk.first,k+jk.second,0);
                }
                //B.printpencil(q,vi);
                
        }
} // end ff2sbundle()

void ff2vbundle(const Field4 &ff, vBundle &B, Int j, Int k, Int extvi=0)
{
        if (extvi==0)
        {
                for (Int q=0;q<(Int)B.stencilSize_;q++)
                {
                        Intpair jk = qtojk(q);
                        for (Int vi=0;vi<(Int)B.nvars_;vi++)
                        {
                                for (Int i=0;i<(Int)B.nx_;i++)
                                {
                                        B(i,q,vi) = ff(i,j+jk.first,k+jk.second,vi);
                        
                                }
                                //B.printpencil(q,vi);
                        }
                }
        }
        else
        {
                for (Int q=0;q<(Int)B.stencilSize_;q++)
                {
                        Intpair jk = qtojk(q);
                        for (Int i=0;i<(Int)B.nx_;i++)
                        {
                                B(i,q,extvi) = ff(i,j+jk.first,k+jk.second,extvi);
                        
                        }
                        
                }
        }
} // end ff2vbundle()



                
Int main()
{


        Int Nx = NX;
        Int Ny = NY;
        Int Nz = NZ;
        
        Real L0 = 0.0;
        Real L1 = 2*M_PI;
        Real dx = (L1-L0)/Nx;
        
        sPencil x(Nx); 
        linspace(x,L0,L1,dx); //x-grid

        Field4 field(Nx,Ny,Nz,3);
        ABC_setup(field,x,x,x);
        printfield(field);
        
        sBundle g(Nx);
        vBundle vg(Nx);
        //ff2sbundle(field,g,(Int)Ny/2,(Int)Nz/2);
        //ff2sbundle(field,g,2,2);
        ff2vbundle(field,vg,2,2);
        
        return 0;

}

/*
  std::cout << "j" << "\t" << "k" << std::endl;
  for (Int q = 0;q<4*NGHOSTS+1;q++)
  {
  Intpair jk = qtojk(q);

  std::cout << jk.first << "\t" << jk.second << std::endl;
  }
*/
