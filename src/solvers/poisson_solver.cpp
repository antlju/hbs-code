#include "includes.h"

/// Frequency division step.
void PssnFreqDiv(fftwMesh &A, const Real xlen)
{
        Int II,JJ;
        Int Nx=A.nx_,Ny=A.ny_,Nzh=A.nzh_;
        Real k1,k2,k3,Pi=M_PI;
        for (Int i=0;i<Nx;i++)
        {
                if (2*i<Nx)
                        II = i;
                else
                        II = Nx-i;
                k1 = 2*Pi*II/xlen;
                
                for (Int j=0;j<Ny;j++)
                {
                        if (2*j<Ny)
                                JJ = j;
                        else
                                JJ = Ny-j;
                        k2 = 2*Pi*JJ/xlen;
                        
                        for (Int k=0;k<Nzh;k++)
                        {
                                k3 = 2*Pi*k/xlen;
                                Real fac = -1.0*(pow(k1,2)+pow(k2,2)+pow(k3,2));
                                if (fabs(fac) < 1e-14)
                                {
                                        A.freqdomRe(i,j,k) = 0.0;
                                        A.freqdomIm(i,j,k) = 0.0;
                                }
                                else
                                {
                                        A.freqdomRe(i,j,k) = A.freqdomRe(i,j,k)/fac;
                                        A.freqdomIm(i,j,k) = A.freqdomIm(i,j,k)/fac;
                                }
                                
                        }
                }
        }
}

/// Solver for 3D Poisson equation using FFTW in-place real-to-complex Fourier transforms.
void PoissonSolve3D(Mesh &inMesh, fftwMesh &fftwMem, const Int vi, const Real xlen)
{
        //Set dimensions.
        Int Nx=inMesh.nx_,Ny=inMesh.ny_,Nz=inMesh.nz_;
        
        // Declare FFTW things.

	fftw_plan fwrd = fftw_plan_dft_r2c_3d(Nx,Ny,Nz,fftwMem.inPtr(),fftwMem.outPtr(),FFTW_MEASURE);
	fftw_plan bwrd = fftw_plan_dft_c2r_3d(Nx,Ny,Nz,fftwMem.outPtr(),fftwMem.inPtr(),FFTW_MEASURE);

        /// Copy from finite difference mesh to FFTW mesh.
        mesh2fftw(inMesh,fftwMem,vi);

        /// Forward transform
        fftw_execute(fwrd);

        /// Perform Pssn frequency division step
        PssnFreqDiv(fftwMem,xlen);

        /// Backward transform
        fftw_execute(bwrd);

        /// Copy from FFTW mesh to finite difference Mesh
        fftw2mesh(fftwMem, inMesh, vi);
}
