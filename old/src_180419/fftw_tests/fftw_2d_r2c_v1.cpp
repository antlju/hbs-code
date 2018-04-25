#include "includes.h"

#include <fftw3.h>

/// This program solves the 2D Poisson equation via the in-place complex-to-complex
/// Fourier Transform and its inverse, using the FFTW library.

void InitialCondition(Real *in, Real q1, Real q2, Pencil x, Pencil y)
{
        //We will use real initial condition so the imaginary part in[][1] = 0.0
        size_t Nx = x.nx_,Ny = y.nx_;
	for (size_t i=0;i<Nx;i++)
	{
		for (size_t j=0;j<Ny;j++)
		{
                        Real val = sin(q1*x(i))*sin(q2*y(j));
                        in[j+(Nx+2)*i] = val; //Real part.
		}
	}
}

/// ------------- MAIN ------------------
Int main()
{
        // Define grid size
        Int Size=32;
        Int Nx=Size,Ny=Size;
        
        // Make space in memory,
        fftw_complex *mem_;
        mem_ = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx * (Ny+2));

        Real *in_;
        in_ = mem_[0];

        fftw_complex *out_;
        out_ = mem_;
        
        // Make fftw plans for forward and backward transforms.
        fftw_plan fplan = fftw_plan_dft_r2c_2d(Nx,Ny,in_,out_,FFTW_MEASURE);
        fftw_plan bplan = fftw_plan_dft_c2r_2d(Nx,Ny,out_,in_,FFTW_MEASURE);

        // Construct grid dimensions and coordinate values
        // We will have symmetric square grid so only need one linspace.
        Real pi = M_PI;
        Real L0 = 0;
        Real L1 = 2*pi;
        Real dx = (L1-L0)/Nx;
        Pencil x(Nx);
        linspace(x,L0,L1,dx);
        
        // Set parameters for initial condition and apply it.
        Real q1=1.0,q2=1.0;
        InitialCondition(in_,q1,q2,x,x);

        // Execute forward transform
        fftw_execute(fplan);

        //Solve FourierSpace Pssn eq
        Int k1,k2;
        for (Int i=0;i<Nx;i++)
        {
                k1 = i-Nx*(Int)(i/(Nx/2+1));
                for (Int j=0;j<Ny/2;j++)
                {
                        k2 = j;
                        Real fac = pow((2*pi*k1/(L1-L0)),2)+pow((2*pi*k2/(L1-L0)),2);
                        if (fac == 0)
                        {
                                out_[j+(Nx)*i][0] = 0.0; //Re
                                out_[j+(Nx)*i][1] = 0.0; //Im
                        }
                        else
                        {
                                out_[j+(Nx)*i][0] = (-1.0/fac)*out_[j+(Nx)*i][0]; //Re
                                out_[j+(Nx)*i][1] = (-1.0/fac)*out_[j+(Nx)*i][1]; //Im
                        }
                }
        }

        // Now we take the inverse transform to find our solution to the Pssn eq.
        // Note that this result is scaled by total grid size, i.e Nx*Ny !
        fftw_execute(bplan);

        //std::cout.precision(4);
        //Printing result of inverse transform
        
        for (Int i=0;i<Nx;i++)
	{
		for (Int j=0;j<Ny;j++)
		{
                        Real repart = (Real)in_[j+(Nx+2)*i]/(Nx*Ny);
                        Real output = repart;
                        
                        std::cout << output << ",";
		}
                std::cout << std::endl;
	}
        
        
        // Remember to destroy plans and free memory
        fftw_destroy_plan(fplan); fftw_destroy_plan(bplan); fftw_free(mem_);


        return 0;
        
        // End main().

}
