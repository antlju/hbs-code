/// This program performs a 2d real-to-complex in place FFT,
/// and then solves the Poisson equation in Fourierspace and then backtransforms.
/// This seems to work for arbitrary parameters q1,q2 and up to at least SIZE 512!!!!

#include "includes.h"

#include <fftw3.h>

#define SIZE 64
#define NXX SIZE
#define NYY SIZE

void initialise(Real *in, Pencil &x, Pencil &y,Real q1,Real q2)
{
	for (size_t i=0;i<x.nx_;i++)
	{
		for (size_t j=0;j<y.nx_;j++)
		{
			//in[j+(x.nx_+2)*i] = cos(q1*x(i));
			in[j+(x.nx_+2)*i] = sin(q1*x(i))*sin(q2*y(j));
		}
	}
}

Int main()
{
  	Real Pi = M_PI;
	// Grid size
	Int Size = SIZE;
	Int Nx=NXX,Ny=NYY,Nyh=(Ny/2+1);


	// Declare FFTW things.
	fftw_complex *mem;
	fftw_complex *out;
	Real *in;
	mem = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * Nx * Nyh);
	out = mem;
	in = mem[0];

	fftw_plan fwrd = fftw_plan_dft_r2c_2d(Nx,Ny,in,out,FFTW_MEASURE);
	fftw_plan bwrd = fftw_plan_dft_c2r_2d(Nx,Ny,out,in,FFTW_MEASURE);
	
	// Parameters
	Real q1=2.0;
	Real q2=3.0;

	Real L0=0.0;
	Real L1=2*Pi;
	Real xlen = (L1-L0);
	Real dx=xlen/Nx;

	// Setup grid and initialise
	Pencil x(Nx);
	linspace(x,L0,L1,dx);
	initialise(in,x,x,q1,q2);

	// Check input
	/*
	std::cout<< "--------------------------------------------" << std::endl;
	for (Int i=0;i<Nx;i++)
	{
		for (Int j=0;j<Ny;j++)
		{
			std::cout << in[j+(Ny+2)*i] << "\t";
		}
		std::cout << std::endl;
	}
	std::cout<< "--------------------------------------------" << std::endl;
	*/
	
	// Execute forward transform
	fftw_execute(fwrd);

	// Check transform before modification
	/*
	std::cout<< "--------------------------------------------" << std::endl;
	for (Int i=0;i<Nx;i++)
	{
		for (Int j=0;j<Nyh;j++)
		{
			std::cout << "(" << out[j+Nyh*i][0] << "," << out[j+Nyh*i][1] << ")" << "\t";
		}
		std::cout << std::endl;
	}
	std::cout<< "--------------------------------------------" << std::endl;
	*/
	
	/// Modify the transformed array to solve poisson in Fourierspace.
	Int II;
	for (Int i=0;i<Nx;i++)
	{
		if (i<Nx/2)
			II = i;
		else
			II = Nx-i;

		Real k1 = 2*Pi*II/xlen;
		for (Int j=0;j<Nyh;j++)
		{
			Real k2 = 2*Pi*j/xlen;

			Real fac = -1.0*(pow(k1,2)+pow(k2,2));
			if (fac == 0)
			{
				out[j+Nyh*i][0] = 0.0;
				out[j+Nyh*i][1] = 0.0;
			}
			else
			{
				out[j+Nyh*i][0] /= fac;
				out[j+Nyh*i][1] /= fac;	
			}
			

			
		}
	}

	// Check transform after modification
	/*
	std::cout<< "--------------------------------------------" << std::endl;
	for (Int i=0;i<Nx;i++)
	{
		for (Int j=0;j<Nyh;j++)
		{
			std::cout << "(" << out[j+Nyh*i][0] << "," << out[j+Nyh*i][1] << ")" << "\t";
		}
		std::cout << std::endl;
	}
	std::cout<< "--------------------------------------------" << std::endl;
	*/
	
	// Execute backward transform
	fftw_execute(bwrd);

	// Check inverse transform
	
	//std::cout<< "--------------------------------------------" << std::endl;
	for (Int i=0;i<Nx;i++)
	{
		for (Int j=0;j<Ny;j++)
		{
			Real ifft = in[j+(Ny+2)*i]/(Nx*Ny);
			Real ainvfac = pow(q1,2)+pow(q2,2);
			Real analytic = -sin(q1*x(i))*sin(q2*x(j))/ainvfac;

			if (fabs(ifft) < 1e-14)
				ifft = 0.0;
			std::cout << ifft << "\t";
			
			Real diff = fabs(ifft-analytic);
			if (diff < 1e-14)
				diff = 0.0;
			//std::cout << diff << "\t";
			
		}
		std::cout << std::endl;
	}
	//std::cout<< "--------------------------------------------" << std::endl;
	

	
	return 0;
}
