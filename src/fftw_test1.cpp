#include "includes.h"

#include <fftw3.h>

void initialise(Real *in, Pencil &x,Pencil &y,Pencil &z)
{
	size_t Nx = x.nx_,Ny = y.nx_,Nz = z.nx_;
	Real q = 1.0;
	for (size_t i=0;i<Nx;i++)
	{
		for (size_t j=0;j<Ny;j++)
		{
			for (size_t k=0;k<Nz;k++)
			{
				Real r2 = pow(x(i),2)+pow(y(j),2)+pow(z(k),2);
				in[k+Ny*(j+Nx*i)] = sin(q*sqrt(r2));
			}
		}
	}
}
	       
Int main()
{
	Int Nx=NX+1,Ny=NY+1,Nz=NZ+1;
	
	//fftwMesh<Real> U(Nx,Ny,Nz);
	fftw_complex *mem_;
	mem_ = (fftw_complex*)
			fftw_malloc(
			sizeof(fftw_complex) * (Nx*Ny*((size_t)(2*(Nz/2+1))))
				);
	
	//fftw_complex *out = mem_;
	Real *in = mem_[0];
	fftw_complex *out = mem_;
	fftw_plan plan = fftw_plan_dft_r2c_3d(Nx,Ny,Nz,
				       in,out,
				       FFTW_MEASURE);

	Real L0=-2*M_PI;
        Real L1=2*M_PI;
        Real dx = (L1-L0)/(Nx-1);
	Pencil x(Nx);
        linspace(x,L0,L1,dx);
	
	initialise(in,x,x,x);


	fftw_execute(plan);
	
	fftw_destroy_plan(plan);
	fftw_free(mem_);

        size_t pNz = 2*(Nz/2+1);
        
        for (size_t i=0;i<Nx;i++)
	{
		for (size_t j=0;j<Ny;j++)
		{
			for (size_t k=0;k<(Nz/2+1);k++)
			{
                                std::cout << "[ " << out[k+Ny*(j+Nx*i)][0] << ", " << out[k+Ny*(j+Nx*i)][1] << "i ]";
			}
                        std::cout << std::endl;
		}
                std::cout << std::endl;
                std::cout << "-----------------------";
                std::cout << std::endl;
	}

	std::cout << "Hej testar fftw!" << std::endl;
        return 0;
}
