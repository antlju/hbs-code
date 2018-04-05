/// This program performs a 3d real-to-complex in place FFT,
/// and then solves the Poisson equation in Fourierspace and then backtransforms.
/// Seem to work up to size 256 at least!


#include "includes.h"

#include <fftw3.h>

#define SIZE 8
#define NXX SIZE
#define NYY SIZE
#define NZZ SIZE

void initialise(Real *in, Pencil &x, Pencil &y, Pencil &z,Real q1,Real q2,Real q3)
{
	for (size_t i=0;i<x.nx_;i++)
	{
		for (size_t j=0;j<y.nx_;j++)
		{ 
                        for (size_t k=0;k<z.nx_;k++)
                        {
			//in[j+(x.nx_+2)*i] = cos(q1*x(i));
                                in[k+(z.nx_+2)*(j+y.nx_*i)] = sin(q1*x(i))*sin(q2*y(j))*sin(q3*z(k));
                        }
		}
	}
}

Int main()
{
  	Real Pi = M_PI;
	// Grid size
	Int Size = SIZE;
	Int Nx=NXX,Ny=NYY,Nz=NZZ,Nzh=(Nz/2+1);


	// Declare FFTW things.
	fftw_complex *mem;
	fftw_complex *out;
	Real *in;
	mem = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * Nx * Ny * Nzh);
	out = mem;
	in = mem[0];

	fftw_plan fwrd = fftw_plan_dft_r2c_3d(Nx,Ny,Nz,in,out,FFTW_MEASURE);
	fftw_plan bwrd = fftw_plan_dft_c2r_3d(Nx,Ny,Nz,out,in,FFTW_MEASURE);
	
	// Parameters
	Real q1=2.0;
	Real q2=q1;
        Real q3=q2;

	Real L0=0.0;
	Real L1=2*Pi;
	Real xlen = (L1-L0);
	Real dx=xlen/Nx;

	// Setup grid and initialise
	Pencil x(Nx);
        linspace(x,L0,L1,dx);
        initialise(in,x,x,x,q1,q2,q3);
        
        for (Int i=0;i<Nx;i++)
	{
		for (Int j=0;j<Ny;j++)
		{
                        for (Int k=0;k<Nz;k++)
                        {
                                Real inn = in[k+(Nz+2)*(j+Ny*i)];
                                if (fabs(inn) < 1e-14)
                                        inn = 0.0;
                                std::cout << inn/ << "\t";

                        }
                        std::cout << std::endl;
		}
                std::cout << "---------------------------------------" << std::endl;
	}
        std::cout << std::endl;
        std::cout << std::endl;
        
        fftw_execute(fwrd);

        Int II,JJ;
        Real k1,k2,k3;
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
                                        out[k+Nzh*(j+Ny*i)][0] = 0.0;
                                        out[k+Nzh*(j+Ny*i)][1] = 0.0;
                                }
                                else
                                {
                                        out[k+Nzh*(j+Ny*i)][0] /= fac;
                                        out[k+Nzh*(j+Ny*i)][1] /= fac;
                                }
                                
                        }
                }
        }
        
        fftw_execute(bwrd);
        
        for (Int i=0;i<Nx;i++)
	{
		for (Int j=0;j<Ny;j++)
		{
                        for (Int k=0;k<Nz;k++)
                        {
                                Real analytic = -sin(q1*x(i))*sin(q2*x(j))*sin(q3*x(k))/(pow(q1,2)+pow(q2,2)+pow(q3,2));
                                Real ifft = in[k+(Nz+2)*(j+Ny*i)]/(Nx*Ny*Nz);
                                Real diff = fabs(ifft-analytic);
                                if (diff < 1e-14)
                                        diff = 0.0;

                                if (fabs(ifft) < 1e-14)
                                        ifft = 0.0;
                                
                                std::cout << ifft << "\t";
                        }
                        std::cout << std::endl;
		}
                std::cout << "---------------------------------------" << std::endl;
	}
	return 0;
}
