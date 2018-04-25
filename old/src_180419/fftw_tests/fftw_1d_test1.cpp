#include "includes.h"

#include <fftw3.h>

void initialise(Real *in, Pencil &x)
{
	size_t Nx = x.nx_;
	Real q = 1.0;
	for (size_t i=0;i<Nx;i++)
	{

                //in[i] = sin(q*x(i));
                in[i] = 0;
	}
        in[Nx/2] = 1.0;
        
}
	       
Int main()
{
	Int Nx=16+1;
        Int pNx = 2*(Nx/2+1);
	
	//fftwMesh<Real> U(Nx,Ny,Nz);
	fftw_complex *mem_;
	mem_ = (fftw_complex*)
                fftw_malloc(
			sizeof(fftw_complex) * 2*(Nx/2+1)
                        );
               
	
//fftw_complex *out = mem_;
        Real *in = mem_[0];
        fftw_complex *out = mem_;
        fftw_plan plan = fftw_plan_dft_r2c_1d(Nx,
                                              in,out,
                                              FFTW_MEASURE);

        Real L0=-2*M_PI;
        Real L1=2*M_PI;
        Real dx = (L1-L0)/(Nx-1);
        Pencil x(Nx);
        linspace(x,L0,L1,dx);
	
        initialise(in,x);
        for (size_t i=0;i<Nx;i++)
	{

                std::cout << in[i] << std::endl;
	}
        std::cout << std::endl;


        fftw_execute(plan);


        for (size_t i=0;i<pNx/2;i++)
	{

                std::cout << out[i][0] << std::endl;
	}










        fftw_destroy_plan(plan);
        fftw_free(mem_);

//std::cout << "Hej testar fftw!" << std::endl;
        return 0;
}
