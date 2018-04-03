#include "includes.h"

#include <fftw3.h>

#define SIZE 2
#define NXX SIZE
#define NYY SIZE

void initialise(Real *in, Pencil &x, Pencil &y,Real q1,Real q2)
{
	for (size_t i=0;i<x.nx_;i++)
	{
		for (size_t j=0;j<y.nx_;j++)
		{
			//in[j+(x.nx_+2)*i] = cos(q1*x(i));
			in[j+(x.nx_)*i] = sin(q1*x(i));
		}
	}
}

Int main()
{
	Real Pi = M_PI;
	// Grid size
	Int Size = SIZE;
	Int Nx=NXX,Ny=NYY;

	// Declare FFTW things. Memory and plans.
	fftw_complex *mem,*mem2;
	double *mem1,*mem3;
	//mem = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*(Ny/2+1));
	mem1 = (double*)fftw_malloc(sizeof(double)*Nx*Ny);
	mem2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*(Ny/2+1));
	mem3 = (double*)fftw_malloc(sizeof(double)*Nx*Ny);
	
	Real *in,*in2;
	fftw_complex *out;
	in = mem1;//mem[0];
	out = mem2;
	in2 = mem3;

	fftw_plan fwrd = fftw_plan_dft_r2c_2d(Nx,Ny,in,out,FFTW_MEASURE);
	fftw_plan bwrd = fftw_plan_dft_c2r_2d(Nx,Ny,out,in2,FFTW_MEASURE);

	// Parameters
	Real q1=1.0;
	Real q2=q1;

	Real L0=0.0;
	Real L1=2*Pi;
	Real dx=(L1-L0)/Nx;

	Pencil x(Nx);
	linspace(x,L0,L1,dx);

	
	initialise(in,x,x,q1,q2);

	for (size_t i=0;i<x.nx_;i++)
	{
		
		for (size_t j=0;j<x.nx_;j++)
		{
			//std::cout << in[j+(x.nx_+2)*i] << " ";
			std::cout << in[j+(x.nx_)*i] << " ";
		}
		std::cout << std::endl;
	}
	
	fftw_execute(fwrd);

	// Fourier space operation

	
	Int II,JJ;
	for (Int i=0;i<Nx;i++)
	{
		Real k1 = 0.0;
		if (i<Nx/2)
			II = i;
		else
			II = Nx-i;
		
		k1 = 2*Pi*II/(L1-L0);
		for (Int j=0;j<Ny/2+1;j++)
		{
			Real k2 = 0.0;
			JJ = j;
			k2 = 2*Pi*JJ/(L1-L0);

			Real Reold = out[j+(Ny)*i][0];
			Real Imold = out[j+(Ny)*i][1];
			std::cout << "(" << Reold << "," << Imold << ")" << "\t";


			out[j+(Ny/2+1)*i][1] = Reold;
			out[j+(Ny/2+1)*i][0] = -Imold;

			std::cout << "(" << out[j+(Ny/2+1)*i][0] << "," << out[j+(Ny/2+1)*i][1] << ")" << "\t";
		}
		std::cout << std::endl;
	}

	
	
	fftw_execute(bwrd);

	std::cout << "-------------------------------" << std::endl;

	for (size_t i=0;i<x.nx_;i++)
	{
		for (size_t j=0;j<x.nx_;j++)
		{
			//std::cout << in[j+(x.nx_+2)*i]/(Nx*Ny) << " ";
			std::cout << in2[j+(x.nx_)*i]/(Nx*Ny) << " ";
		}
		std::cout << std::endl;
	}

	fftw_destroy_plan(fwrd); fftw_destroy_plan(bwrd); fftw_free(mem1); fftw_free(mem2); fftw_free(mem3);
}
