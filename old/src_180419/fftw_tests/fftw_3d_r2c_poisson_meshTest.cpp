
#include "includes.h"

#include <fftw3.h>

#define SIZE 4
#define NXX SIZE
#define NYY SIZE
#define NZZ SIZE

typedef fftwMesh<Real> fftwM;

void printIn(fftwM &A,Int ifft)
{
        Real scalefac;
        if (ifft == 1)
                scalefac = A.nx_*A.ny_*A.nz_;
        else
                scalefac = 1.0;
        std::cout << "--------------------------------------" << std::endl;
	for (size_t i=0;i<A.nx_;i++)
	{
		for (size_t j=0;j<A.ny_;j++)
		{ 
                        for (size_t k=0;k<A.nz_;k++)
                        {

                                if (fabs(A.in(i,j,k)) < 1e-14)
                                        std::cout << 0.0 << "\t";
                                else
                                std::cout << A.in(i,j,k)/scalefac << "\t";
                        }
                        std::cout << std::endl;

		}
                std::cout << std::endl;
	}
        std::cout << "--------------------------------------" << std::endl;
}

void printOut(fftwM &A)
{
        std::cout << "--------------------------------------" << std::endl;
	for (size_t i=0;i<A.nx_;i++)
	{
		for (size_t j=0;j<A.ny_;j++)
		{ 
                        for (size_t k=0;k<A.nzh_;k++)
                        {
                                Real outRe = A.outRe(i,j,k);
                                Real outIm = A.outIm(i,j,k);
                                if (fabs(outRe) < 1e-14)
                                        outRe = 0.0;
                                if (fabs(outIm) < 1e-14)
                                        outIm = 0.0;

                                std::cout << "(" << outRe << "," << outIm << ")" << "\t";
                        }
                        std::cout << std::endl;

		}
                std::cout << std::endl;
	}
        std::cout << "--------------------------------------" << std::endl;
}

void initialise(fftwM &A, Pencil &x, Pencil &y, Pencil &z,Real q1,Real q2,Real q3)
{
	for (size_t i=0;i<A.nx_;i++)
	{
		for (size_t j=0;j<A.ny_;j++)
		{ 
                        for (size_t k=0;k<A.nz_;k++)
                        {
                                
                                A.in(i,j,k) = sin(q1*x(i))*sin(q2*y(j))*sin(q3*z(k));
                        }

		}

	}
}

void PssnFreqDiv(fftwM &A, Real xlen)
{
        Int II,JJ;
        Int Nx=A.nx_,Ny=A.ny_,Nzh=A.nzh_;
        Real k1,k2,k3,Pi=M_PI;
        for (Int i=0;i<Nx;i++)
        {
                if (2*i<Nx)
                        II = i;
                else
                        II = Nx-i;
                k1 = 2*Pi*II/xlen;
                
                for (Int j=0;j<Ny;j++)
                {
                        if (2*j<Ny)
                                JJ = j;
                        else
                                JJ = Ny-j;
                        k2 = 2*Pi*JJ/xlen;
                        
                        for (Int k=0;k<Nzh;k++)
                        {
                                k3 = 2*Pi*k/xlen;
                                Real fac = -1.0*(pow(k1,2)+pow(k2,2)+pow(k3,2));
                                if (fabs(fac) < 1e-14)
                                {
                                        A.outRe(i,j,k) = 0.0;
                                        A.outIm(i,j,k) = 0.0;
                                }
                                else
                                {
                                        A.outRe(i,j,k) = A.outRe(i,j,k)/fac;
                                        A.outIm(i,j,k) = A.outIm(i,j,k)/fac;
                                }
                                
                        }
                }
        }
}

Int main()
{
  	Real Pi = M_PI;
	// Grid size
	Int Nx=NXX,Ny=NYY,Nz=NZZ;


	// Declare FFTW things.
        fftwM fftwMem(Nx,Ny,Nz);

	fftw_plan fwrd = fftw_plan_dft_r2c_3d(Nx,Ny,Nz,fftwMem.inPtr(),fftwMem.outPtr(),FFTW_MEASURE);
	fftw_plan bwrd = fftw_plan_dft_c2r_3d(Nx,Ny,Nz,fftwMem.outPtr(),fftwMem.inPtr(),FFTW_MEASURE);
	
	// Parameters
	Real q1=1.0;
	Real q2=q1;
        Real q3=q2;

	Real L0=0.0;
	Real L1=2*Pi;
	Real xlen = (L1-L0);
	Real dx=xlen/Nx;

	// Setup grid and initialise
	Pencil x(Nx);
        linspace(x,L0,L1,dx);
        initialise(fftwMem,x,x,x,q1,q2,q3);
        printIn(fftwMem,0);

        fftw_execute(fwrd);
        printOut(fftwMem);
        
        PssnFreqDiv(fftwMem,xlen);
        fftw_execute(bwrd);
        
        printIn(fftwMem,1);
        
	return 0;
}
