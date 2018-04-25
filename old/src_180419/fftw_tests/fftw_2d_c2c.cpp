/// ---------------------------------------
///
/// In place complex-to-complex FFT -> IFFT
///
/// ---------------------------------------

#include "includes.h"

#include <fftw3.h>

#define SIZE 2
#define NXX SIZE
#define NYY SIZE

void initialise(fftw_complex *in, Pencil &x, Pencil &y,Real q1,Real q2)
{
	Int Nx = x.nx_,Ny=y.nx_;
	for (Int i=0;i<Nx;i++)
	{
		for (Int j=0;j<Ny;j++)
		{
			Int idx = j+Ny*i;
			in[idx][0] = cos(x(i));
			in[idx][1] = 0.0;
		}
	}
}

void printC(const fftw_complex *in, Int Nx, Int Ny)
{
	std::cout << "------------------ print ------------------" << std::endl;
	for (Int i=0;i<Nx;i++)
	{
		for (Int j=0;j<Ny;j++)
		{
			Int idx = j+Ny*i;
			std::cout << idx <<": "<< "("<<in[idx][0]<<","<<in[idx][1]<<")"<<std::endl;
		}
	}
	std::cout << "-------------------------------------------" << std::endl;
}

Int main()
{
	Real Pi = M_PI;
	// Grid size
	//Int Size = SIZE;
	Int Nx=NXX,Ny=NYY;//Nyhf=Ny/2+1;
	
	// Declare FFTW things. Memory and plans.
	fftw_complex *mem;
	mem = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);

	
	fftw_complex *in = mem;
	fftw_complex *out = mem;

	
	fftw_plan fwrd = fftw_plan_dft_2d(Nx,Ny,
					  in,out,
					  FFTW_FORWARD,FFTW_MEASURE);
	fftw_plan bwrd = fftw_plan_dft_2d(Nx,Ny,
					  out,in,
					  FFTW_BACKWARD,FFTW_MEASURE);
	
	// Parameters
	Real q1=1.0;
	Real q2=q1;

	Real L0=0.0;
	Real L1=2*Pi;
	Real dx=(L1-L0)/Nx;

	Pencil x(Nx);
	linspace(x,L0,L1,dx);
	//x.printpencil();
	//printC(in,Nx,Ny);
	
	initialise(in,x,x,q1,q2);

	//printC(in,Nx,Ny);

	fftw_destroy_plan(fwrd); fftw_destroy_plan(bwrd);
	
	return 0;
}
