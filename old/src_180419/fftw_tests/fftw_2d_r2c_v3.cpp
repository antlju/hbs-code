#include "includes.h"

#include <fftw3.h>

/// This program solves the 2D Poisson equation via the in-place complex-to-complex
/// Fourier Transform and its inverse, using the FFTW library.

void InitialCondition(Real *in, Real q1, Real q2, const Pencil x, const Pencil y)
{
        //We will use real initial condition so the imaginary part in[][1] = 0.0
        size_t Nx = x.nx_,Ny = y.nx_;
	for (size_t i=0;i<Nx;i++)
	{
		for (size_t j=0;j<Ny;j++)
		{
                        Real val = sin(q2*y(j))*sin(q1*x(i));
                        in[j+(Nx+2)*i] = val; //Real part.
		}
	}
}

void ksetup(fftw_complex *outArr, const Int Nx, const Int Ny, const Real L1,const Real L0)
{
        Real k1,k2;
        Int II,JJ;
        Real pi = M_PI;
        for (Int i=0;i<Nx;i++)
        {
                //k1 = (Int)i-Nx*(i/(Nx/2+1));
                if (i<Nx/2)
                        II = i;
                else
                        II = Nx-i;
                
                k1 = 2*pi*II/(L1-L0);
                for (Int j=0;j<Ny/2+1;j++)
                {
                        JJ = j;
                        k2 = 2*pi*JJ/(L1-L0);
                        //partial_x v => -i*k_x*vhat in Fourier space.
                        //What we have from forward transform is vhat in the "out" array.
                        //So we need to set this to -i*k_x*vhat.

                        //Re(out) = -k_x*Im(out_old)
                        Real Imoutold = outArr[j+Nx*i][1];
                        Real Reoutold = outArr[j+Nx*i][0];
                        outArr[j+Nx*i][0] = k1*Imoutold;
                        outArr[j+Nx*i][1] = -1.0*k1*Reoutold;
                }
                
        }
}

/// ------------- MAIN ------------------
Int main()
{
        // Define grid size
        Int Size=16;
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
        Real q1=1.0;
        Real q2=q1;
        InitialCondition(in_,q1,q2,x,x);

        // Execute forward transform
        fftw_execute(fplan);

        //Apply partial x deriv in Fourier space.
        ksetup(out_,Nx,Ny,L1,L0);
        

        // Now we take the inverse transform to find our solution to the Pssn eq.
        // Note that this result is scaled by total grid size, i.e Nx*Ny !
        fftw_execute(bplan);

        //std::cout.precision(4);
        //Printing result of inverse transform
        
        for (Int i=0;i<Nx;i++)
	{
                for (Int j=0;j<Ny;j++)
                {
                        Real output = in_[j+(Nx+2)*i]/((Size/4)*Nx*Ny);
                        std::cout << output << " ";
                }
                std::cout << std::endl;
	}
        
        
        // Remember to destroy plans and free memory
        fftw_destroy_plan(fplan); fftw_destroy_plan(bplan); fftw_free(mem_);


        return 0;
        
        // End main().

}
