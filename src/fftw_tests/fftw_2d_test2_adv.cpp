#include "includes.h"

#include <fftw3.h>

void initialise(Real *in, Real q1, Real q2,Pencil &x,Pencil &y)
{
	size_t Nx = x.nx_,Ny = y.nx_;
	for (size_t i=0;i<Nx;i++)
	{
		for (size_t j=0;j<Ny;j++)
		{
                        Real val = cos(q1*x(i))*cos(q2*y(j));
                        in[j+Nx*i] = val;
		}
	}
}
	       
Int main()
{
        Real q1=1.0,q2=1.0;
        Real pi = M_PI;
        Int Size = 4;
	Int Nx=Size+1,Ny=Size+1; 
	fftw_complex *mem_;
	mem_ = (fftw_complex*)
			fftw_malloc(
                                sizeof(fftw_complex) * (Nx*(2*(Ny/2+1)))
				);
	
	//fftw_complex *out = mem_;
	Real *in = mem_[0];
	fftw_complex *out = mem_;
        Int dims[2] = {Nx,Ny};
        Int inembed1[2] = {Nx,2*(Ny/2+1)}; Int onembed1[2] = {Nx,Ny/2+1};
        Int inembed2[2] = {Nx,Ny};
        Int *nullPtr = NULL;
	fftw_plan forward = fftw_plan_many_dft_r2c(2,dims,1,
                                                   in,inembed1,1,1,
                                                   out,onembed1,1,1,
                                                   FFTW_MEASURE);
        fftw_plan backward = fftw_plan_many_dft_c2r(2,dims,1,
                                                   out,onembed1,1,1,
                                                   in,inembed1,1,1,
                                                   FFTW_MEASURE);
        

	Real L0=-2*M_PI;
        Real L1=2*M_PI;
        Real dx = (L1-L0)/(Nx);
	Pencil x(Nx);
        linspace(x,L0,L1,dx);
        Real dxx = x(1)-x(0);
        
	
	initialise(in,q1,q2,x,x);
        

	fftw_execute(forward);
        fftw_execute(backward);
 
        std::cout.precision(4);
        //Printing result of inverse transform compared to analytical
        for (size_t i=0;i<Nx;i++)
	{
		for (size_t j=0;j<2*(Ny/2+1);j++)
		{
                        //Real ana_psi = -cos(q1*x(i))*cos(q2*x(j))/(pow(q1,2)+pow(q2,2));
                        Real ana_A = cos(q1*x(i))*cos(q2*x(j));
                        Real diff = fabs(in[j+Nx*i]/(Nx*Ny)-ana_A);
                        if (diff < 1e-14)
                                diff = 0.0;
                        
                        std::cout << std::left << diff << "\t";
		}
                std::cout << std::endl;
	}
        
	fftw_destroy_plan(forward);
        fftw_destroy_plan(backward);
	fftw_free(mem_);








	std::cout << "Hej testar 2D fftw!" << std::endl;
        return 0;
}
