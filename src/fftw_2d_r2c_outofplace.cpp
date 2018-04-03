/// This program performs a 2d real-to-complex out of place FFT.
/// This works for at least FFT->iFFT up to SIZE 512!

#include "includes.h"

#include <fftw3.h>

#define SIZE 512
#define NXX SIZE
#define NYY SIZE

void initialise(Real *in, Pencil &x, Pencil &y,Real q1,Real q2)
{
	for (size_t i=0;i<x.nx_;i++)
	{
		for (size_t j=0;j<y.nx_;j++)
		{
			//in[j+(x.nx_+2)*i] = cos(q1*x(i));
			in[j+(x.nx_)*i] = sin(q1*x(i))*cos(q2*y(j));
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
	fftw_complex *out;
	Real *in;
	out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * Nx * Nyh);
	in = (Real*)fftw_malloc(sizeof(Real) * Nx * Ny);

	fftw_plan fwrd = fftw_plan_dft_r2c_2d(Nx,Ny,in,out,FFTW_MEASURE);
	fftw_plan bwrd = fftw_plan_dft_c2r_2d(Nx,Ny,out,in,FFTW_MEASURE);
	
	// Parameters
	Real q1=1.0;
	Real q2=q1;

	Real L0=0.0;
	Real L1=2*Pi;
	Real dx=(L1-L0)/Nx;

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
			std::cout << in[j+Ny*i] << "\t";
		}
		std::cout << std::endl;
	}
	std::cout<< "--------------------------------------------" << std::endl;
	*/
	
	// Execute forward transform
	fftw_execute(fwrd);

	//Check transform
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
	/*
	std::cout<< "--------------------------------------------" << std::endl;
	for (Int i=0;i<Nx;i++)
	{
		for (Int j=0;j<Ny;j++)
		{
			std::cout << in[j+Ny*i]/(Nx*Ny) << "\t";
		}
		std::cout << std::endl;
	}
	std::cout<< "--------------------------------------------" << std::endl;
	*/

	std::cout<< "--------------------------------------------" << std::endl;
	for (Int i=0;i<Nx;i++)
	{
		for (Int j=0;j<Ny;j++)
		{
			std::cout << abs((in[j+Ny*i]/(Nx*Ny))-sin(q1*x(i))*cos(q2*x(j))) << "\t";
		}
		std::cout << std::endl;
	}
	std::cout<< "--------------------------------------------" << std::endl;
	
	return 0;
}
