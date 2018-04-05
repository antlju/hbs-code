#include "includes.h"

#include <fftw3.h>

#define SIZE 32
#define NXX SIZE
#define NYY SIZE
#define NZZ SIZE
/// This program computes the partial 1st derivative w.r.t coordinate i of
/// a 3D scalar function, by Fourier Transform.
/// The point is to learn how the in-place real-to-complex FFTW implementation works.

Int idx(Int i,Int j,Int k)
{
        return k+(NYY+2)*(j+NXX*i);
}

void setIn(Real *in,Pencil &x,Pencil &y,Pencil &z,Real q1)
{
        for (Int i=0;i<NXX;i++)
        {
                for (Int j=0;j<NYY;j++)
                {
                        for (Int k=0;k<NZZ;k++)
                        {
                                in[idx(i,j,k)] = sin(x(i));
                        }
               
                }
               
        }
}

void printOut(const fftw_complex *out,Int Nx,Int Ny,Int Nz)
{
        std::cout << "Print start --------------------------------" << std::endl;
        for (Int i=0;i<Nx;i++)
        {
                for (Int j=0;j<Ny;j++)
                {
                        for (Int k=0;k<Nz/2+1;k++)
                        {
                                std::cout << "(" <<out[k+Ny*(j+Nx*i)][0] << ","
                                          << out[k+Ny*(j+Nx*i)][1] << ")" << " ";
                        }
                        std::cout << std::endl;
                }
                std::cout << std::endl;
        }
        std::cout << "Print end --------------------------------" << std::endl;
}


void printIn(const Real *in,Int Nx,Int Ny,Int Nz,Int isIFFT)
{
        std::cout << "Print start --------------------------------" << std::endl;
        for (Int i=0;i<Nx;i++)
        {
                for (Int j=0;j<Ny;j++)
                {
                        for (Int k=0;k<Nz+2;k++)
                        {
                                if (isIFFT == 1)
                                        std::cout << in[idx(i,j,k)]/(Nx*Ny*Nz) << " ";
                                else
                                        std::cout << in[idx(i,j,k)] << " ";
                        }
                        std::cout << std::endl;
                }
                std::cout << std::endl;
        }
        std::cout << "Print end --------------------------------" << std::endl;
}

void ksetup(fftw_complex *outArr,Int Nx,Int Ny,Int Nz,Real L1,Real L0)
{
        Real k1,k2,k3;
        Int II,JJ,KK;
        Real pi = M_PI;
        for (Int i=0;i<Nx;i++)
        {
                //k1 = (Int)i-Nx*(i/(Nx/2+1));
                if (i<Nx/2)
                        II = i;
                else
                        II = Nx-i;
                
                k1 = 2*pi*II/(L1-L0);
                for (Int j=0;j<Ny;j++)
                {
                        //k2 = (Int)j-Ny*(j/(Ny/2+1));
                        if (j<Ny/2)
                                JJ = j;
                        else
                                JJ = Ny-j;
                        
                        k2 = 2*pi*JJ/(L1-L0);
                        for (Int k=0;k<Nz/2+1;k++)
                        {
                                KK = k;
                                k3 = 2*pi*KK/(L1-L0);
                                //partial_x v => -i*k_x*vhat in Fourier space.
                                //What we have from forward transform is vhat in the "out" array.
                                //So we need to set this to -i*k_x*vhat.

                                //Re(out) = -k_x*Im(out_old)
                                Real Imoutold = outArr[k+Ny*(j+Nx*i)][1];
                                Real Reoutold = outArr[k+Ny*(j+Nx*i)][0];
                                outArr[k+Ny*(j+Nx*i)][0] = k1*Imoutold;
                                outArr[k+Ny*(j+Nx*i)][1] = -k1*Reoutold;
                        }
                }
        }
}

Int main()
{
        // Define grid size
        Int Size=SIZE;
        Int Nx=Size,Ny=Size,Nz=Size;

        // Allocate memory
        fftw_complex *mem_;
        mem_ = (fftw_complex*) fftw_malloc(
                sizeof(fftw_complex)*(Nx*Ny*(Nz+2))
                );

        Real *in;
        fftw_complex *out;
        in = mem_[0];
        out = mem_;

        // Make fftw plans for forward and backward transforms.
        fftw_plan fplan = fftw_plan_dft_r2c_3d(Nx,Ny,Nz,in,out,FFTW_MEASURE);
        fftw_plan bplan = fftw_plan_dft_c2r_3d(Nx,Ny,Nz,out,in,FFTW_MEASURE);
        
        //printIn(in,Nx,Ny,Nz);

        // Construct grid dimensions and coordinate values
        // We will have symmetric square grid so only need one linspace.
        Real pi = M_PI;
        Real L0 = 0;
        Real L1 = 2*pi;
        Real dx = (L1-L0)/Nx;
        Pencil x(Nx);

        linspace(x,L0,L1,dx);
        //x.printpencil();

        // Set parameters for initial condition and apply it.
        Real q1 = 1.0;
        setIn(in,x,x,x,q1);

        //printIn(in,Nx,Ny,Nz);        

        // Execute forward transform
        fftw_execute(fplan);

        printOut(out,Nx,Ny,Nz);
        
        // Set-up k = (k1,k2,k3).
        // This is really where the derivative operation happens in Fourier space.
        ksetup(out,Nx,Ny,Nz,L1,L0);

        printOut(out,Nx,Ny,Nz);
        
        // Execute backward transform
        fftw_execute(bplan);

        //Print result
        for (Int i=0;i<Nx;i++)
                std::cout << in[idx(i,0,0)]/((SIZE/4)*Nx*Ny*Nz) <<"\t" << cos(x(i)) << std::endl;
        //printIn(in,Nx,Ny,Nz,1);
        // Remember to destroy plans and free memory
        fftw_destroy_plan(fplan); fftw_destroy_plan(bplan); fftw_free(mem_);
        
        return 0;
        
}
// End main().
