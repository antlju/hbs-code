/// This program performs a 2d real-to-complex out of place FFT.
/// It then solves the 2D poisson equation in fourierspace and backtransforms.
/// This seems to work up to SIZE 512!!!!

#include "includes.h"

#include <fftw3.h>

#define SIZE 8
#define NXX SIZE
#define NYY SIZE

void initialise(Real *in, Pencil &x, Pencil &y,Real q1,Real q2)
{
	for (size_t i=0;i<x.nx_;i++)
	{
		for (size_t j=0;j<y.nx_;j++)
		{
			//in[j+(x.nx_+2)*i] = cos(q1*x(i));
			in[j+(x.nx_)*i] = sin(q1*x(i))*sin(q2*x(j));
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
	Real q1=2.0;
	Real q2=q1;

	Real L0=0.0;
	Real L1=2*Pi;
	Real xlen = (L1-L0);
	Real dx=xlen/Nx;

	// Setup grid and initialise
	Pencil x(Nx);
	linspace(x,L0,L1,dx);
	initialise(in,x,x,q1,q2);

	// Check input
	
	std::cout<< "IN AFTER INIT --------------------------------------------" << std::endl;
	for (Int i=0;i<Nx;i++)
	{
		for (Int j=0;j<Ny;j++)
		{
			Real inin = in[j+Ny*i];
			if (fabs(inin) < 1e-14)
				inin = 0.0;
			std::cout << inin << "\t";
		}
		std::cout << std::endl;
	}
	std::cout<< "--------------------------------------------" << std::endl;
	
	
	// Execute forward transform
	fftw_execute(fwrd);

	//Check FFT before modification
	/*
	std::cout<< "OUT BEFORE MOD --------------------------------------------" << std::endl;
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

	//Check FFT after modification
	/*
	std::cout<< "OUT AFTER MOD --------------------------------------------" << std::endl;
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
	std::cout<< "IN AFTER IFFT--------------------------------------------" << std::endl;
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

	// Check difference of result from analytic
	std::cout<< "DIFF AFTER IFFT--------------------------------------------" << std::endl;
	for (Int i=0;i<Nx;i++)
	{
		for (Int j=0;j<Ny;j++)
		{
			Real ifft = in[j+Ny*i]/(Nx*Ny);
			Real invanafactor = pow(q1,2)+pow(q2,2);
			Real analytic = -sin(q1*x(i))*sin(q2*x(j))/invanafactor;
			Real diff = fabs(ifft-analytic);
			if (diff < 1e-14)
				diff = 0.0;
			std::cout << diff << "\t";
			if (fabs(ifft) < 1e-14)
				ifft = 0.0;
			//std::cout << invanafactor*ifft << "\t";
		}
		std::cout << std::endl;
	}
	std::cout<< "--------------------------------------------" << std::endl;

	return 0;
}
