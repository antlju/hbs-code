#include "includes.h"

// Copy should only be done to and from the "space domain".
void mesh2fftw(const Mesh &input, fftwMesh &out, const Int vi=0)
{
        // Make sure that dimensions match.
        assert(input.nx_ == out.nx_);
        assert(input.ny_ == out.ny_);
        assert(input.nz_ == out.nz_);

        //Copy loop
        for (size_t i=0;i<input.nx_;i++)
        {
                for (size_t j=0;j<input.ny_;j++)
                {
                        for (size_t k=0;k<input.nz_;k++)
                        {
                                out.spacedom(i,j,k) = input(i,j,k,vi);
                        }
                }
        }
                
}

// Copy should only be done to and from the "space domain".
void fftw2mesh(const fftwMesh &input, Mesh &out, const Int vi=0)
{
        // Make sure that dimensions match.
        assert(input.nx_ == out.nx_);
        assert(input.ny_ == out.ny_);
        assert(input.nz_ == out.nz_);

        //Copy loop
        for (size_t i=0;i<input.nx_;i++)
        {
                for (size_t j=0;j<input.ny_;j++)
                {
                        for (size_t k=0;k<input.nz_;k++)
                        {
                                out(i,j,k,vi) = input.spacedom(i,j,k);
                        }
                }
        }
                
}


Int main()
{
        Int Nx=NX,Ny=NY,Nz=NZ;

        Real fillvalA = 1.0;
        Real fillvalB = 3.0;
        
        // Test copy from mesh to fftw
        fftwMesh FFTW_A(Nx,Ny,Nz);
        Mesh FMESH_A(Nx,Ny,Nz,1);
        FMESH_A.fill(fillvalA);
        FFTW_A.printspacedom();
        mesh2fftw(FMESH_A,FFTW_A);
        FFTW_A.printspacedom();

        // Test copy from fftw to mesh
        fftwMesh FFTW_B(Nx,Ny,Nz);
        Mesh FMESH_B(Nx,Ny,Nz,1);
        FFTW_B.fill(fillvalB);
        FMESH_B.print();
        fftw2mesh(FFTW_B,FMESH_B);
        FMESH_B.print();
        
        
}
