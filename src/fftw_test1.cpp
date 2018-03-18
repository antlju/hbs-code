#include "includes.h"

#include <fftw3.h>

template<class T>
class fftwMesh {
public:
        /// Internal storage.
        //std::vector<T> mem_;
	fftw_complex *mem_;

        /// "Real" x-size
        size_t nx_;

        /// "Real" y-size
        size_t ny_;

        /// "Real" z-size
        size_t nz_;

	/// Logical size of the "padding" (last) dimension
	size_t pdimsize_;

	/// 3D constructor
	fftwMesh(size_t Nx, size_t Ny, size_t Nz) :
		nx_(Nx), ny_(Ny), nz_(Nz), pdimsize_((size_t)(2*(Nz/2+1)))
	{
		mem_ = (fftw_complex*)
			fftw_malloc(
			sizeof(fftw_complex) * (Nx*Ny*((size_t)(2*(Nz/2+1))))
				);
		
		//mem_.resize( Nx*Ny*((size_t)(2*(Nz/2+1))) ); //Allocate memory
		// Note how last dimension is padded for in-place Real-to-Complex transform
                //std::fill(mem_.begin(),mem_.end(), T() ); //Fill with zeros of type T.
	}

	T* dataPtr()
	{
		return &mem_;
	}
	
};

void initialise(Real *in, Pencil &x,Pencil &y,Pencil &z)
{
	size_t Nx = x.nx_,Ny = y.nx_,Nz = z.nx_;
	
	for (size_t i=0;i<Nx;i++)
	{
		for (size_t j=0;j<Ny;j++)
		{
			for (size_t k=0;k<Nz;k++)
			{
				Real r2 = pow(x(i),2)+pow(y(j),2)+pow(z(k),2);
				in[i+Ny*(j+Nz*k)] = sin(sqrt(r2));
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

	std::cout << "Hej testar fftw!" << std::endl;
        return 0;
}
