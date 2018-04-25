#include "includes.h"

#include <fftw3.h>

/// This program solves the 2D Poisson equation via the in-place complex-to-complex
/// Fourier Transform and its inverse, using the FFTW library.

void InitialCondition(fftw_complex *in, Real q1, Real q2, Pencil x, Pencil y)
{
        //We will use real initial condition so the imaginary part in[][1] = 0.0
        size_t Nx = x.nx_,Ny = y.nx_;
	for (size_t i=0;i<Nx;i++)
	{
		for (size_t j=0;j<Ny;j++)
		{
                        Real val = sin(q1*x(i))*sin(q2*y(j));
                        in[j+Nx*i][0] = val; //Real part.
                        in[j+Nx*i][1] = 0.0; //Im part.
		}
	}
}

/// ------------- MAIN ------------------
Int main()
{
        // Define grid size
        Int Size=64;
        Int Nx=Size,Ny=Size;
        
        // Make space in memory,
        fftw_complex *mem_;
        mem_ = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx * Ny);

        // Make fftw plans for forward and backward transforms.
        fftw_plan fplan = fftw_plan_dft_2d(Nx,Ny,mem_,mem_,FFTW_FORWARD,FFTW_MEASURE);
        fftw_plan bplan = fftw_plan_dft_2d(Nx,Ny,mem_,mem_,FFTW_BACKWARD,FFTW_MEASURE);

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
        InitialCondition(mem_,q1,q2,x,x);

        // Execute forward transform
        fftw_execute(fplan);

        // To solve the FTransformed Pssn eq we need to solve for f = -1/(k1^2+k2^2)*g.
        // This amounts to multiplying our transform by -1/(k1^2+k2^2).
        // We find k1_i = 2*pi*i/(L1-L0) and k2_j = 2*pi*j/(L1-L0).
        
        for (Int i=0;i<Nx;i++)
        {
                for (Int j=0;j<Ny;j++)
                {
                        
                        Int II = i;
                        Int JJ = j;
                        
                        if (2*i<Nx)
                                II = i;
                        else
                                II = Nx-i;

                        if (2*j<Ny)
                                JJ = j;
                        else
                                JJ = Ny-j;
                        
                        
                        Real fac = pow((2*pi*II/(L1-L0)),2)+pow((2*pi*JJ/(L1-L0)),2);

                        if (fac == 0)
                        {
                                mem_[j+Nx*i][0] = 0.0; //Re
                                mem_[j+Nx*i][1] = 0.0; //Im
                        }
                        else
                        {
                                mem_[j+Nx*i][0] = (-1.0/fac)*mem_[j+Nx*i][0]; //Re
                                mem_[j+Nx*i][1] = (-1.0/fac)*mem_[j+Nx*i][1]; //Im
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
                        Real repart = (Real)mem_[j+(Nx)*i][0]/(Nx*Ny);
                        Real impart = (Real)mem_[j+(Nx)*i][1]/(Nx*Ny);
                        Real out = repart;
                        
                        std::cout << out << ",";
		}
                std::cout << std::endl;
	}
        


        // Remember to destroy plans and free memory
        fftw_destroy_plan(fplan); fftw_destroy_plan(bplan); fftw_free(mem_);


        return 0;
} // End main().
